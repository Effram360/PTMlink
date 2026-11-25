#!/usr/bin/env python3
"""ptmlink -- a thin Python reader of the PTMlink release, mirroring the R package function for function.

It reads the same per-gene shards the website and the R package read (a local data directory, or
the release address with checksum verification), applies the same default view in the same order,
sorts on the same one key (evidence), reads the same reading rules from the same reference file, and
recomputes nothing. The repository's parity test runs both readers on one fixture and holds them to
the same rows, the same withheld counts, the same edges and the same rule verdicts.

Standard library only for the shards; pyarrow (optional) for query() over site_served.parquet.

    from ptmlink import PTMlink
    db = PTMlink()        # the website's data; PTMlink(data_dir="path/to/data") reads a local copy instead
    db.status(["EGFR", "HER2"])
    rows = db.gene("EGFR", "glyco")["rows"]
    db.pair("TFRC_N727", "COIL_T303")["cohorts"]
    db.explain("TFRC_N727", "COIL_T303")
    db.query(where=[("q_value", "<", 1e-10), ("k_cohorts", ">=", 5)], limit=100)

Command line (a readable summary on stdout; --json prints the full result):

    python3 ptmlink.py status EGFR HER2
    python3 ptmlink.py gene EGFR --side glyco --include-cis
    python3 ptmlink.py pair TFRC_N727 COIL_T303
    python3 ptmlink.py pair "PTMG:P02786:N727|PTMP:P38432:T303"        # the same pair, by its pair_id
    python3 ptmlink.py --release-dir path/to/release query --where q_value '<' 1e-10 --where k_cohorts '>=' 5 --limit 20

Optional: --data-dir path/to/data reads a local copy of the website's data instead of downloading it.
"""
import argparse
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import sys
import tempfile
import textwrap
import time
import urllib.error
import urllib.parse
import urllib.request
from statistics import median

DEFAULT_FILTERS = {
    "include_cis": False, "min_k": 2, "hide_composition": False, "anchored_only": False, "significant_only": True,
    "max_i2": None, "adjusted_excludes_zero": False, "shadow_clean_only": False, "hide_self_coupled": False, "min_loco_pctile": None,
}
WITHHELD_KEYS = ["cis", "k", "composition", "anchored", "significant", "heterogeneous", "loco", "adjusted", "shadow", "self"]
NETWORK_COLUMNS = ["gene", "partner", "role", "rank", "drawn", "n_site_pairs", "n_own_sites", "n_partner_sites",
                   "n_site_pairs_positive", "n_site_pairs_negative", "sign_mixed", "evidence_best", "evidence_median", "pctile_at_best",
                   "rho_at_best", "k_at_best", "n_at_best", "I2_at_best", "q_at_best", "q_anchored_at_best",
                   "rho_prot", "shadow_clean_at_best", "self_coupling_at_best", "self_glyco_at_best", "self_phospho_at_best",
                   "n_site_pairs_anchored", "anch_any", "anch_dose", "anch_gene_editing", "anch_enzymatic",
                   "own_site_at_best", "best_glyco_site", "best_phospho_site"]
NOTE_SYMMETRY = "the score is symmetric and licenses no direction"


# Some hosts refuse urllib's default User-Agent (the live site's firewall returns 403), so every download names the reader.
USER_AGENT = "ptmlink-python/1.0.0"
# With no source given, the reader uses the public website's data, as the R package does.
DEFAULT_DATA_URL = "https://ptmlink.biomarkercenter.org/data"


def _site_id_map():
    """site key -> stable identifier for every reported site, from the map the R package ships
    (ptmlink/inst/extdata/site_ids.tsv.gz); a shard alone knows only its own gene's sites."""
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(here, "..", "ptmlink", "inst", "extdata", "site_ids.tsv.gz"), os.path.join(here, "site_ids.tsv.gz")):
        if os.path.exists(cand):
            with gzip.open(cand, "rt", encoding="utf-8") as fh:
                next(fh, None)
                return dict(line.rstrip("\n").split("\t", 1) for line in fh if line.strip())
    return {}


# THREE KINDS OF FAILURE, KEPT APART, as in the R package. A file the data source does not list (DataAbsent) is a
# fact about the release: a gene with no reported pair has no data file. A host that cannot be reached
# (DataUnreachable) or a download that does not match the checksum list (ChecksumMismatch) is a fact about the
# connection, and is never reported as the absence of a gene or a pair. Callers catch DataAbsent alone.
class PTMlinkError(RuntimeError):
    """An error the reader explains in one sentence; the command line and the MCP server print it as it is."""


class DataAbsent(PTMlinkError):
    pass


class DataUnreachable(PTMlinkError):
    pass


class ChecksumMismatch(PTMlinkError):
    pass


# ONE WORDING for a pair the release does not report, in R, Python and the MCP server
ABSENT_PAIR = "not among the reported pairs of this release; the complete scored pairs are released with the article"


def _host(url):
    h = urllib.parse.urlparse(url).netloc
    return h or url


def _unreachable(base, detail, answered=False):
    if answered:   # the host answered with an HTTP status: the address, not the connection, is the first thing to check
        return DataUnreachable("%s answered %s. If the address is wrong, set it with --data-url (data_url= in Python): it is the folder that "
                               "holds SHA256SUMS.json. Otherwise try again in a few minutes, or use a local copy of the data with --data-dir "
                               "(data_dir= in Python)." % (_host(base), detail))
    return DataUnreachable("could not reach %s%s. Check the internet connection and try again, or use a local copy of the data "
                           "with --data-dir (data_dir= in Python)." % (_host(base), " (%s)" % detail if detail else ""))


def _download(url, dest, base=None):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=60) as r, open(dest, "wb") as fh:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
    except urllib.error.HTTPError as e:
        raise _unreachable(base or url, "%s %s" % (e.code, e.reason), answered=True) from None
    except urllib.error.URLError as e:
        raise _unreachable(base or url, str(getattr(e, "reason", e))) from None
    except OSError as e:   # timeouts, resets
        raise _unreachable(base or url, str(e)) from None


class AmbiguousName(ValueError):
    """An alias that more than one symbol carries."""


def _not_symbol(q):
    if re.match(r"^.+_[NSTY][0-9]+$", str(q).strip().upper()):
        raise ValueError("%s is a site key, not a gene; the site command (site() in Python, ptm_site in MCP) lists its partners" % q)
    raise ValueError("%s is not a recognized symbol, alias or UniProt accession; status %s says how symbols resolve" % (q, q))


def _is_na(v):
    return v is None or (isinstance(v, float) and math.isnan(v))


GENE_EXAMPLE = 'a gene symbol, alias or UniProt accession such as "EGFR"'
SITE_EXAMPLE = 'a site key such as "EGFR_N603" or "EGFR_T693"'


def _names(x, arg, example, one=False):
    """A name argument as a list of non-empty strings; anything else is refused with the argument's name and an
    example, before any file is read (the R package's .ptm_chr)."""
    xs = [x] if isinstance(x, str) else (list(x) if isinstance(x, (list, tuple)) else None)
    if xs is None:
        raise ValueError("%s is %s, not text: give %s" % (arg, "missing" if x is None else type(x).__name__, example))
    if not xs:
        raise ValueError("%s is empty: give %s" % (arg, example))
    bad = [v for v in xs if not isinstance(v, str)]
    if bad:
        raise ValueError("%s holds %s, not text: give %s" % (arg, "a missing entry" if bad[0] is None else type(bad[0]).__name__, example))
    xs = [v.strip() for v in xs]
    if any(not v for v in xs):
        raise ValueError("%s holds a missing or empty entry: give %s" % (arg, example))
    if one and len(xs) != 1:
        raise ValueError("%s takes one value here, got %d" % (arg, len(xs)))
    return xs


def check_filters(f):
    """The view's settings, checked: a partial dict is completed from the defaults, an unknown name or a value of the
    wrong kind is refused (the R package's .ptm_filters)."""
    out = dict(DEFAULT_FILTERS)
    if f is None:
        return out
    if not isinstance(f, dict):
        raise ValueError('filters is a set of named settings, such as {"min_k": 1}')
    bad = [k for k in f if k not in DEFAULT_FILTERS]
    if bad:
        raise ValueError("unknown filter %s; the filters are %s" % (", ".join('"%s"' % k for k in bad), ", ".join(DEFAULT_FILTERS)))
    for k, v in f.items():
        if isinstance(DEFAULT_FILTERS[k], bool):
            if not isinstance(v, bool):
                raise ValueError("%s is true or false" % k)
        elif k == "min_k":
            if isinstance(v, bool) or not isinstance(v, (int, float)) or v != int(v) or v < 1:
                raise ValueError("min_k is a whole number of cohorts, 1 or more")
            v = int(v)
        else:
            if _is_na(v):
                v = None
            elif isinstance(v, bool) or not isinstance(v, (int, float)) or not 0 <= v <= 1:
                raise ValueError("%s is a number between 0 and 1, or null for no limit" % k)
        out[k] = v
    return out


def no_pair_note(st):
    """What was measured for a gene that has no reported pair, in one line."""
    return "%s has no reported pair. Glycosites: %s. Phosphosites: %s." % (st.get("gene"), st.get("status_glyco"), st.get("status_phospho"))


def _reference_path():
    p = os.environ.get("PTMLINK_REFERENCE")
    if p and os.path.exists(p):
        return p
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(here, "..", "ptmlink", "inst", "reference.json"), os.path.join(here, "reference.json")):
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError("reference.json not found; set PTMLINK_REFERENCE to ptmlink/inst/reference.json")


# The label a reader sees. The build record keeps the internal build label; the twelve-cohort build v1.2 is the first
# public release, v1.0.0, as on the website. A build record that states `public_release` wins.
PUBLIC_LABEL = {"v1.2": "v1.0.0"}


def _release_date():
    """The public release date, as the package's citation file states it."""
    try:
        p = os.path.join(os.path.dirname(_reference_path()), "extdata", "ptmlink_citation.cff")
        with open(p, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("date-released:"):
                    return line.split(":", 1)[1].strip().strip('"')
    except OSError:
        pass
    return None


def public_label(release, public=None):
    if public:
        return public
    return PUBLIC_LABEL.get(release, release)


# The article's display names, as on the website: GC for the gastric cohort, ccRCC, rare kidney
# cancers, and NAT for normal adjacent tissue. "cohort" keeps the analysis key the release files and every argument use.
COHORT_LABEL = {"CCRCC": "ccRCC", "RARE_KIDNEY": "Rare kidney cancers", "STAD": "GC"}
TUMOUR_TYPE = {"STAD": "gastric cancer", "RARE_KIDNEY": "non-clear-cell renal cancer"}


def cohort_label(code):
    return COHORT_LABEL.get(code, code)


def display_cohort(c):
    d = dict(c, label=c.get("label") or cohort_label(c.get("cohort")))
    if "tumour_type" in c:
        base = TUMOUR_TYPE.get(c.get("cohort")) or re.sub(r"\s*\(rare kidney histologies\)", "",
                                                        re.sub(r",\s*normal adjacent tissue( \(NAT\))?$", "", str(c.get("tumour_type") or "")))
        d["tumour_type"] = base + ", normal adjacent tissue (NAT)" if d.get("role") == "normal" else base
    return d


class PTMlink:
    def __init__(self, data_dir=None, data_url=None, release_dir=None, cache_dir=None, verify_release=True):
        if data_dir and not os.path.isdir(data_dir):
            raise FileNotFoundError("no such data directory: %s" % data_dir)
        self.data_dir = os.path.abspath(data_dir) if data_dir else None
        self.data_url = (data_url or os.environ.get("PTMLINK_DATA_URL", "") or ("" if data_dir else DEFAULT_DATA_URL)).rstrip("/")
        self.cache_dir = cache_dir or os.environ.get("PTMLINK_CACHE_DIR") or os.path.join(os.path.expanduser("~"), ".cache", "ptmlink", "data")
        if self.data_url:
            self.cache_dir = os.path.join(self.cache_dir, hashlib.sha256(self.data_url.encode("utf-8")).hexdigest())
        self._remote_sums = None
        self.release_dir = None
        self.release_sha = None
        self._cache = {}
        self._ids = {}
        self._ids_release_read = False
        self._ids_package_read = False
        self._reference = None
        self.messages = []
        self.echo_messages = True
        if release_dir:
            self.set_release_dir(release_dir, verify=verify_release)

    # ---- data sources -------------------------------------------------------------------------
    def set_release_dir(self, path, verify=True):
        if not os.path.isdir(path):
            raise FileNotFoundError("no such release directory: %s" % path)
        pq = os.path.join(path, "site_served.parquet")
        if verify:
            if not os.path.exists(pq):
                raise FileNotFoundError("no site_served.parquet in " + path)
            sums = os.path.join(path, "SHA256SUMS.json")
            if not os.path.exists(sums):
                raise FileNotFoundError("no SHA256SUMS.json in " + path + "; pass verify=False to use the directory unverified")
            want = json.load(open(sums)).get("site_served.parquet")
            got = _sha256(pq)
            if want is None or got != want:
                raise ValueError("site_served.parquet in %s does not match SHA256SUMS.json" % path)
            self.release_sha = got
        self.release_dir = os.path.abspath(path)
        self._ids_release_read = False
        return self.release_dir

    def _manifest(self):
        # Refresh once per reader instance, including corrected releases at the same URL.
        if self._remote_sums is None:
            os.makedirs(self.cache_dir, exist_ok=True)
            fd, tmp = tempfile.mkstemp(dir=self.cache_dir)
            os.close(fd)
            cached = os.path.join(self.cache_dir, "SHA256SUMS.json")
            try:
                try:
                    _download(self.data_url + "/SHA256SUMS.json", tmp, self.data_url)
                except DataUnreachable:
                    # offline with a cache: the checksum list saved by the last session still verifies every cached file
                    old = None
                    if os.path.exists(cached):
                        try:
                            with open(cached, encoding="utf-8") as fh:
                                old = json.load(fh)
                        except ValueError:
                            old = None
                    if not isinstance(old, dict) or not old:
                        raise
                    self._say("could not reach %s; using the checksum list cached on %s (files not cached yet cannot be read)"
                              % (_host(self.data_url), time.strftime("%Y-%m-%d", time.localtime(os.path.getmtime(cached)))))
                    self._remote_sums = old
                    return old
                with open(tmp, encoding="utf-8") as fh:
                    sums = json.load(fh)
                if not isinstance(sums, dict) or not sums:
                    raise RuntimeError("invalid checksum manifest at " + self.data_url)
                os.replace(tmp, os.path.join(self.cache_dir, "SHA256SUMS.json"))
                self._remote_sums = sums
            finally:
                if os.path.exists(tmp):
                    os.remove(tmp)
        return self._remote_sums

    def _say(self, text):
        """A message, as R's message(): to stderr, and kept in self.messages for a caller (the MCP server) to return."""
        self.messages.append(text)
        if self.echo_messages:
            sys.stderr.write("ptmlink: %s\n" % text)

    def _has(self, rel):
        """Does the data source carry this file? Asks the checksum list (remote) or the directory (local)."""
        if self.data_dir:
            return os.path.exists(os.path.join(self.data_dir, rel))
        if not self.data_url:
            return False
        return rel in self._manifest()

    def _path(self, rel):
        # a name built from what the caller typed (a gene, a site) never leaves the data source
        if rel.startswith(("/", "\\")) or "\\" in rel or ".." in rel.split("/") or "\x00" in rel:
            raise DataAbsent("%s is not in this data source" % rel)
        if self.data_dir:
            p = os.path.join(self.data_dir, rel)
            if not os.path.exists(p):
                raise DataAbsent("%s is not in this data source" % rel)
            return p
        if not self.data_url:
            raise RuntimeError("no data source: give data_dir or data_url, or set PTMLINK_DATA_URL")
        sums = self._manifest()
        local = os.path.join(self.cache_dir, rel)
        if rel == "SHA256SUMS.json":
            return local
        want = sums.get(rel)
        if want is None:
            raise DataAbsent("%s is not in this data source" % rel)
        os.makedirs(os.path.dirname(local), exist_ok=True)
        if not os.path.exists(local) or _sha256(local) != want:
            fd, tmp = tempfile.mkstemp(dir=os.path.dirname(local))
            os.close(fd)
            try:
                _download(self.data_url + "/" + rel, tmp, self.data_url)
                try:
                    self._verify(tmp, rel)
                except ChecksumMismatch:
                    # the site may have been updated since this reader read its checksum list (a long-running server):
                    # read the list again, once, and accept the file only if the new list names it
                    self._remote_sums = None
                    if self._manifest().get(rel) == want:
                        raise
                    self._cache.clear()
                    _download(self.data_url + "/" + rel, tmp, self.data_url)
                    self._verify(tmp, rel)
                os.replace(tmp, local)
            finally:
                if os.path.exists(tmp):
                    os.remove(tmp)
        return local

    def _verify(self, local, rel):
        want = self._manifest().get(rel)
        if want is None:
            os.remove(local)
            raise DataAbsent("%s is not in this data source" % rel)
        if _sha256(local) != want:
            os.remove(local)
            raise ChecksumMismatch("checksum mismatch for %s -- download removed. The site may be in the middle of an update; "
                                   "try again in a few minutes." % rel)

    def _json(self, rel):
        if rel in self._cache:
            return self._cache[rel]
        p = self._path(rel)
        opener = gzip.open if p.endswith(".gz") else open
        with opener(p, "rt", encoding="utf-8") as fh:
            v = json.load(fh)
        self._cache[rel] = v
        if rel.startswith("genes/"):
            self._learn_ids(v)
        return v

    @staticmethod
    def _inference(meta=None):
        """Release metadata, never the release name, says the method: student_normal_bh (Student tails, q_value BH, the
        only correction) or the legacy Fisher-normal/BH of earlier internal builds written before the marker. Any
        other marker reads as legacy, as the release build itself reads it."""
        meta = meta or {}
        meta = meta.get("inference") or meta
        modern = meta.get("inference_revision") == "student_normal_bh"
        out = {"inference_revision": "student_normal_bh" if modern else "legacy_fisher_bh",
               "p_method": "signed_student_normal_weighted" if modern else "fisher_normal",
               "q_method": "BH", "q_family": "all_scored_site_pairs",
               "anchored_q_family": "perturbation_anchored_rows_within_each_layer",
               "served_rule": "q_value<0.05 OR anch_any"}
        out.update({k: meta[k] for k in ("q_family", "anchored_q_family", "glycoform_q_family", "cohort_q_family")
                    if meta.get(k) is not None})
        return out

    def _shard(self, gene):
        sh = self._json("genes/%s.json.gz" % gene)
        sh["inference"] = self._inference(sh)
        return sh

    def _shard_or_none(self, gene):
        try:
            return self._shard(gene)
        except DataAbsent:
            return None

    def _structure(self, gene):
        """A gene's AlphaFold record, found as the website finds it: through structure/bundles.json (gene -> bundle file,
        each bundle an object keyed by gene) when the source carries the index, else structure/<GENE>.json.gz. None
        when the source has no model for the gene; a network or checksum failure propagates."""
        if self._has("structure/bundles.json"):
            path = (self._json("structure/bundles.json").get("genes") or {}).get(gene)
            return self._json(path).get(gene) if path else None
        try:
            return self._json("structure/%s.json.gz" % gene)
        except DataAbsent:
            return None

    # ---- identifiers --------------------------------------------------------------------------
    def _learn_ids(self, sh):
        for tab in (sh.get("glyco_sites") or [], sh.get("phospho_sites") or []):
            for rec in tab:
                if rec.get("site") and rec.get("row_id"):
                    self._ids[rec["site"]] = rec["row_id"]

    def _ids_from_release(self):
        if not self.release_dir or self._ids_release_read:
            return
        for f in ("features_glyco.tsv", "features_phospho.tsv"):
            p = os.path.join(self.release_dir, f)
            if os.path.exists(p):
                with open(p, encoding="utf-8") as fh:
                    head = fh.readline().rstrip("\n").split("\t")
                    if "site" in head and "row_id" in head:
                        si, ri = head.index("site"), head.index("row_id")
                        for line in fh:
                            cells = line.rstrip("\n").split("\t")
                            if len(cells) > max(si, ri) and cells[si] and cells[ri] and cells[ri] != "NA":
                                self._ids[cells[si]] = cells[ri]
        self._ids_release_read = True

    def _add_ids(self, rows):
        if not self._ids_package_read:
            for k, v in _site_id_map().items():
                self._ids.setdefault(k, v)
            self._ids_package_read = True
        self._ids_from_release()
        for r in rows:
            if _is_na(r.get("glyco_row_id")):
                r["glyco_row_id"] = self._ids.get(r["glyco_site"])
            if _is_na(r.get("phospho_row_id")):
                r["phospho_row_id"] = self._ids.get(r["phospho_site"])
            r["pair_id"] = None if _is_na(r["glyco_row_id"]) or _is_na(r["phospho_row_id"]) else r["glyco_row_id"] + "|" + r["phospho_row_id"]
        return rows

    # ---- small readers ------------------------------------------------------------------------
    def release(self):
        b = self._json("build.json")
        return dict(b, inference=self._inference(b), build_release=b.get("release"),
                    release=public_label(b.get("release"), b.get("public_release")), released=_release_date())

    def cohorts(self):
        """One record per cohort with its role: "tumour" (a tumor cohort, pooled) or "normal" (view-only, never pooled
        or counted); a release whose file has no role holds tumor cohorts only."""
        return [display_cohort(dict(c, role=c.get("role") or "tumour")) for c in self._json("cohorts.json")]

    def _kcov(self):
        """cohort -> k_cov of the tumor cohorts the dictionary states; none stated gives every cohort k = 1."""
        try:
            d = self._json("cohorts.json")
        except Exception:
            return {}
        return {c["cohort"]: float(c["k_cov"]) for c in d if c.get("k_cov") is not None and c.get("role") != "normal"}

    def glycans(self):
        return self._json("glycans.json.gz")

    @staticmethod
    def default_filters():
        return dict(DEFAULT_FILTERS)

    def _index(self):
        return self._json("index.json.gz")

    def _resolve(self, gene):
        idx = self._index()
        up = gene.strip().upper()
        for g in idx["genes"]:
            if g.get("gene") and g["gene"].upper() == up:
                return {"gene": g["gene"], "via": "current symbol"}
        hits = [a for a in idx["aliases"] if a.get("query") and a.get("resolves_to") and a["query"].upper() == up]
        # an alias listed twice for one symbol (two letter cases, or an alias that is also a previous symbol) names that symbol
        opts = sorted(set(h["resolves_to"] for h in hits))
        if len(opts) == 1:
            return {"gene": opts[0], "via": "%s of %s" % (re.sub(r" of$", "", hits[0]["kind"]), opts[0])}
        if len(opts) > 1:
            raise AmbiguousName("%s is ambiguous: it is an alias of %s; give one of these symbols" % (gene, " and ".join([", ".join(opts[:-1]), opts[-1]])))
        return {"gene": None, "via": "unknown"}

    def status(self, gene):
        genes = _names(gene, "gene", GENE_EXAMPLE)
        out = []
        for q in genes:
            try:
                r = self._resolve(q)
            except AmbiguousName as e:
                why = str(e).split(": ", 1)[1]
                out.append({"query": q, "gene": None, "resolved_via": "ambiguous: " + why, "status_glyco": "ambiguous alias", "status_phospho": "ambiguous alias"})
                continue
            if r["gene"] is None:
                via = "a site key, not a gene (ptm_site() lists its partners)" if re.match(r"^.+_[NSTY][0-9]+$", q.upper()) else "not a recognized symbol"
                out.append({"query": q, "gene": None, "resolved_via": via, "status_glyco": "not a recognized symbol", "status_phospho": "not a recognized symbol"})
                continue
            rec = next((g for g in self._index()["genes"] if g.get("gene") == r["gene"]), None)
            if rec is None:
                out.append({"query": q, "gene": r["gene"], "resolved_via": r["via"], "status_glyco": "not in the index", "status_phospho": "not in the index"})
                continue
            row = {"query": q, "resolved_via": r["via"]}
            row.update(rec)
            for k in ("status_glyco", "status_phospho"):   # the index's word, in the readers' terms
                if isinstance(row.get(k), str):
                    row[k] = row[k].replace("no served pair", "no reported pair")
            out.append(row)
        return out

    # ---- the default view ---------------------------------------------------------------------
    @staticmethod
    def n_total(packed):
        if _is_na(packed) or not packed:
            return None
        return sum(int(c.split(":")[-1]) for c in packed.split(";") if c)

    @staticmethod
    def _glyco_labelled(glyco_programs, rep_glycans):
        if _is_na(glyco_programs) or not glyco_programs:
            return False
        segs = [s for s in glyco_programs.split("|") if s]
        attributed = any("=" in s for s in segs)
        if not attributed or _is_na(rep_glycans) or not rep_glycans:
            return True
        used = [c.split(":", 1)[-1] for c in rep_glycans.split(";")]
        for s in segs:
            gly, _, body = s.partition("=")
            if gly in used and body:
                return True
        return False

    def _filter(self, rows, f=None):
        f = check_filters(f)
        withheld = {k: 0 for k in WITHHELD_KEYS}
        rows = [dict(r) for r in rows]
        for r in rows:
            if "cohorts_rho" in r and "n_total" not in r:
                r["n_total"] = self.n_total(r.get("cohorts_rho"))
        keep = [True] * len(rows)

        def charge(reason, pred):
            for i, r in enumerate(rows):
                if keep[i] and pred(r):
                    withheld[reason] += 1
                    keep[i] = False

        if not f["include_cis"]:
            charge("cis", lambda r: r.get("is_cis") is True)
        if f["shadow_clean_only"]:
            charge("shadow", lambda r: r.get("shadow_clean") is not True)
        if f["hide_self_coupled"]:
            charge("self", lambda r: not _is_na(r.get("self_coupling")) and r["self_coupling"] >= 0.6)
        charge("k", lambda r: _is_na(r.get("k_cohorts")) or r["k_cohorts"] < f["min_k"])
        if f["hide_composition"]:
            charge("composition", lambda r: self._glyco_labelled(r.get("glyco_programs"), r.get("rep_glycans")) and bool(r.get("phospho_programs") or ""))
        if f["anchored_only"]:
            charge("anchored", lambda r: r.get("anch_any") is not True)
        if f["significant_only"]:
            charge("significant", lambda r: _is_na(r.get("q_value")) or not (r["q_value"] < 0.05))
        if f["max_i2"] is not None:
            charge("heterogeneous", lambda r: not _is_na(r.get("I2")) and r["I2"] > f["max_i2"])
        if f["min_loco_pctile"] is not None:
            charge("loco", lambda r: _is_na(r.get("loco_pctile_k")) or not (r["loco_pctile_k"] >= f["min_loco_pctile"]))
        if f["adjusted_excludes_zero"]:
            charge("adjusted", lambda r: _is_na(r.get("rho_adjusted_lo")) or _is_na(r.get("rho_adjusted_hi")) or not (r["rho_adjusted_lo"] > 0 or r["rho_adjusted_hi"] < 0))
        out = [r for i, r in enumerate(rows) if keep[i]]
        out.sort(key=lambda r: -r["evidence"])   # stable, as R's order()
        return self._add_ids(out), withheld

    # ---- rows ---------------------------------------------------------------------------------
    def gene(self, gene, side="glyco", filters=None):
        if side not in ("glyco", "phospho"):
            raise ValueError("side is glyco or phospho, not %r" % (side,))
        genes = _names(gene, "gene", GENE_EXAMPLE)
        filters = check_filters(filters)
        if len(genes) > 1:   # several genes: stacked with a query_gene column; a name that is not recognized is said and skipped
            parts, seen = [], set()
            for q in genes:
                try:
                    r = self._resolve(q)
                except ValueError as e:
                    r = {"gene": None, "via": re.sub(r"^[^:]*: ", "", str(e))}
                if r["gene"] is not None and r["gene"] in seen:   # the same protein named twice: read once
                    continue
                seen.add(r["gene"])
                if r["gene"] is None:
                    why = ("a site key, not a gene (the site command lists its partners)" if re.match(r"^.+_[NSTY][0-9]+$", q.upper())
                           else "not a recognized symbol, alias or UniProt accession" if r["via"] == "unknown" else r["via"])
                    self._say("%s skipped: %s" % (q, why))
                    continue
                part = self.gene(q, side, filters)
                if part.get("note"):
                    self._say(part["note"])
                    continue
                parts.append(part)
            if not parts:
                raise ValueError("none of the %d names is a gene with a reported pair: %s" % (len(genes), "; ".join(self.messages[-len(genes):])))
            return {"genes": [p["gene"] for p in parts], "side": side, "rows": [dict(r, query_gene=p["gene"]) for p in parts for r in p["rows"]],
                    "withheld": {p["gene"]: p["withheld"] for p in parts}, "sites": {p["gene"]: p["sites"] for p in parts}, "filters": filters}
        gene = genes[0]
        r = self._resolve(gene)
        if r["gene"] is None:
            _not_symbol(gene)
        st = self.status(r["gene"])[0]
        if (st.get("n_served_as_glyco", 0) + st.get("n_served_as_phospho", 0)) == 0:
            return {"gene": r["gene"], "side": side, "resolved_via": r["via"], "rows": [], "withheld": {k: 0 for k in WITHHELD_KEYS}, "sites": [], "filters": filters, "note": no_pair_note(st)}
        sh = self._shard(r["gene"])
        rows = [x for x in sh["rows"] if (x["glyco_gene"] if side == "glyco" else x["phospho_gene"]) == r["gene"]]
        rows, withheld = self._filter(rows, filters)
        return {"inference": sh["inference"], "gene": r["gene"], "side": side, "resolved_via": r["via"], "rows": rows, "withheld": withheld,
                "sites": sh["glyco_sites" if side == "glyco" else "phospho_sites"], "filters": filters}

    def _canon_gene(self, gene):
        """A gene name as the index spells it (C10orf71, not C10ORF71); an alias or accession that names one symbol is
        read as that symbol (HER2_N530 is ERBB2_N530); an unknown name comes back in upper case."""
        up = str(gene).strip().upper()
        if getattr(self, "_gene_case", None) is None:
            try:
                idx = self._index()
                self._gene_case = {g["gene"].upper(): g["gene"] for g in reversed(idx["genes"]) if g.get("gene")}
                targets = {}
                for a in idx.get("aliases") or []:
                    if a.get("query") and a.get("resolves_to"):
                        targets.setdefault(a["query"].upper(), set()).add(a["resolves_to"])
                self._alias_one = {q: next(iter(t)) for q, t in targets.items() if len(t) == 1}
            except Exception:   # no index to read: the lookup that follows reports it
                return up
        return self._gene_case.get(up) or self._alias_one.get(up, up)

    def _canon_site(self, site):
        """A typed site key in the release's own letter case: the gene part as the index spells it, the residue in
        upper case. A site identifier (PTMG:..., PTMP:...) and anything that is not a site key come back in upper case."""
        up = str(site).strip().upper()
        m = re.match(r"^(.+)_([NSTY][0-9]+)$", up)
        if not m or re.match(r"^PTM[GP]:", up):
            return up
        return "%s_%s" % (self._canon_gene(m.group(1)), m.group(2))

    def _site_key(self, site):
        m = re.match(r"^(.+)_([NSTY])([0-9]+)$", self._canon_site(site))
        if not m:
            raise ValueError("%s is not a site key (GENE_N123 or GENE_S123)" % site)
        return {"gene": m.group(1), "key": "%s_%s%s" % (m.group(1), m.group(2), m.group(3)), "side": "glyco" if m.group(2) == "N" else "phospho"}

    def site(self, site, filters=None):
        sites = _names(site, "site", SITE_EXAMPLE)
        filters = check_filters(filters)
        if len(sites) > 1:   # several sites: stacked with a query_site column; a site in no reported pair is said and skipped
            parts, seen = [], set()
            for q in sites:
                key = self._site_key(q)["key"]
                if key in seen:   # the same site named twice: read once
                    continue
                seen.add(key)
                try:
                    parts.append(self.site(q, filters))
                except DataAbsent as e:
                    self._say(str(e))
            if not parts:
                raise DataAbsent("none of the %d sites is in a reported pair of this release" % len(seen))
            return {"sites": [p["site"] for p in parts], "rows": [dict(r, query_site=p["site"]) for p in parts for r in p["rows"]],
                    "withheld": {p["site"]: p["withheld"] for p in parts}, "filters": filters}
        k = self._site_key(sites[0])
        sh = self._shard_or_none(k["gene"])
        rows = [x for x in sh["rows"] if x["glyco_site"] == k["key"] or x["phospho_site"] == k["key"]] if sh else []
        if not rows:   # a site with no reported pair is said as such, never shown as "no pair passes the filters"
            raise DataAbsent("%s is in no reported pair of this release; %s" % (k["key"], ("the sites command lists the reported sites of %s" % k["gene"]) if sh else
                                                                                 ("status %s says what was measured" % k["gene"])))
        rows, withheld = self._filter(rows, filters)
        return {"inference": sh["inference"], "site": k["key"], "rows": rows, "withheld": withheld, "filters": filters}

    # ---- naming a pair ----------------------------------------------------------------------------
    @staticmethod
    def _split_pair_id(x, fn="pair"):
        """A pair named by its stable identifier, "PTMG:<accession>:N<position>|PTMP:<accession>:<residue><position>", as
        the one argument; split into its two site identifiers (read as sites by _pair_keys)."""
        xs = [x] if isinstance(x, str) else list(x or [])
        xs = [str(v).strip().upper() for v in xs]
        bad = [v for v in xs if not re.match(r"^PTMG:[^|]+\|PTMP:[^|]+$", v)]
        if not xs or bad:
            raise ValueError('%s takes a glycosite and a phosphosite (IL6ST_N131 TMX1_S270) or one pair_id '
                             '("PTMG:P40189:N131|PTMP:Q9H3N1:S270"); got %s' % (fn, bad[0] if bad else "nothing"))
        g = [v.split("|")[0] for v in xs]
        p = [v.split("|")[1] for v in xs]
        return (g[0], p[0]) if len(xs) == 1 else (g, p)

    @staticmethod
    def _key_side(k):
        return "glyco" if re.search(r"_N[0-9]+$", k) else ("phospho" if re.search(r"_[STY][0-9]+$", k) else None)

    def _pair_keys(self, gs, ps):
        """The site keys a pair argument stands for: a site key as given, or every key a site identifier (row_id) maps to
        in the shipped site list. A phosphosite given first is read in the right order, with a note."""
        def keys(x):
            if not re.match(r"^PTM[GP]:", x):
                return [x]
            self._add_ids([])
            k = sorted({site for site, rid in self._ids.items() if rid == x})
            if not k:
                raise ValueError("%s is not a site identifier of this release" % x)
            return k
        g, p = keys(gs), keys(ps)
        sg, sp = {self._key_side(k) for k in g}, {self._key_side(k) for k in p}
        if None in sg | sp:
            raise ValueError("%s is not a site key (GENE_N123 or GENE_S123)" % next(k for k in g + p if self._key_side(k) is None))
        note = None
        if sg == {"phospho"} and sp == {"glyco"}:
            g, p = p, g
            note = "the glycosite comes first: reading the pair as %s x %s" % (g[0], p[0])
        elif sg != {"glyco"} or sp != {"phospho"}:
            raise ValueError("a pair is one glycosite (GENE_N123) and one phosphosite (GENE_S123, GENE_T123 or GENE_Y123); got %s and %s" % (g[0], p[0]))
        return g, p, note

    def pair(self, glyco_site, phospho_site=None):
        if not isinstance(glyco_site, (str, list, tuple)) or (phospho_site is not None and not isinstance(phospho_site, (str, list, tuple))):
            raise ValueError('pair takes a glycosite and a phosphosite (IL6ST_N131 TMX1_S270) or one pair_id ("PTMG:P40189:N131|PTMP:Q9H3N1:S270")')
        glyco_site = (_names(glyco_site, "glyco_site", SITE_EXAMPLE) if isinstance(glyco_site, (list, tuple)) else _names(glyco_site, "glyco_site", SITE_EXAMPLE, one=True)[0])
        if phospho_site is None:
            glyco_site, phospho_site = self._split_pair_id(glyco_site, "pair")
        else:
            phospho_site = (_names(phospho_site, "phospho_site", SITE_EXAMPLE) if isinstance(phospho_site, (list, tuple)) else _names(phospho_site, "phospho_site", SITE_EXAMPLE, one=True)[0])
            if any("|" in v for v in ([glyco_site] if isinstance(glyco_site, str) else glyco_site)):
                raise ValueError("glyco_site holds a pair_id, which names both sites: leave phospho_site out")
        if isinstance(glyco_site, (list, tuple)) or isinstance(phospho_site, (list, tuple)):
            gl = list(glyco_site) if isinstance(glyco_site, (list, tuple)) else [glyco_site]
            pl = list(phospho_site) if isinstance(phospho_site, (list, tuple)) else [phospho_site]
            if not gl or not pl:
                raise ValueError("give at least one glycosite and one phosphosite")
            if len(gl) != len(pl) and 1 not in (len(gl), len(pl)):
                raise ValueError("glyco_site and phospho_site pair up element by element: give the same number of each, or one of either (got %d and %d)" % (len(gl), len(pl)))
            n = max(len(gl), len(pl))
            if n == 1:
                return self.pair(gl[0], pl[0])
            # several pairs: phospho_site is recycled to the length of glyco_site, as in the R package
            return {"pairs": [self.pair(gl[i % len(gl)], pl[i % len(pl)]) for i in range(n)]}
        gl, pl, input_note = self._pair_keys(self._canon_site(glyco_site), self._canon_site(phospho_site))
        gs, ps = gl[0], pl[0]
        gene = re.sub(r"_N[0-9]+$", "", gs)
        absent = {"glyco_site": gs, "phospho_site": ps, "row": None, "cohorts": None, "n_total": None, "rep_glycans": None, "glycoforms": None, "note": ABSENT_PAIR}
        if input_note:
            absent["input_note"] = input_note
        sh = self._shard_or_none(gene)
        if sh is None:
            return absent
        # the partner's data file teaches the partner site's identifier, so the pair id is complete
        pgene = self._site_key(ps)["gene"]
        if pgene != gene:
            self._shard_or_none(pgene)
        rows = [x for x in sh["rows"] if x["glyco_site"] in gl and x["phospho_site"] in pl]
        if not rows:
            return absent
        row = dict(max(rows, key=lambda x: x["evidence"]))
        gs, ps = row["glyco_site"], row["phospho_site"]
        self._add_ids([row])
        if "n_total" not in row:
            row["n_total"] = self.n_total(row.get("cohorts_rho"))
        cohorts = []
        kc = self._kcov()
        for cell in row["cohorts_rho"].split(";"):
            c, rho, n = cell.split(":")
            rho, n = float(rho), int(n)
            z = math.atanh(min(max(rho, -0.999999), 0.999999))
            k = kc.get(c, 1)   # the interval at the cohort's own weight n - 3 - k (cohorts.json k_cov; 1 when unstated)
            se = 1 / math.sqrt(n - 3 - k) if n - 3 - k > 0 else float("nan")
            lo, hi = math.tanh(z - 1.959964 * se), math.tanh(z + 1.959964 * se)   # to four decimals, the precision rho is released at
            cohorts.append({"cohort": c, "label": cohort_label(c), "rho": rho, "n": n, "lo": None if math.isnan(lo) else round(lo, 4),
                            "hi": None if math.isnan(hi) else round(hi, 4), "k_cov": k})
        reps = None
        if row.get("rep_glycans"):
            reps = [{"cohort": c.split(":")[0], "rep_glycan": c.split(":")[1]} for c in row["rep_glycans"].split(";")]
        glycoforms = None
        try:
            gf = self._json("glycoforms/%s.json.gz" % gs)
        except DataAbsent:
            gf = None
        if gf and gf.get("pairs", {}).get(ps):
            cols = gf["columns"]
            glycoforms = []
            for cells in gf["pairs"][ps]:
                m = dict(zip(cols, cells))
                glycoforms.append({**{k: m[k] for k in ("sum_weight", "q_value_anchored") if k in m}, "glycan": m["glycan"], "rho": m["rho"], "evidence": m["evidence"], "k": m["k"], "I2": m["I2"], "q": m["q"], "anchored": bool(m["anch"])})
        out = {"inference": sh["inference"], "row": row, "cohorts": cohorts, "n_total": sum(c["n"] for c in cohorts), "rep_glycans": reps, "glycoforms": glycoforms, "glycoform_inference": self._inference(gf) if gf else None, "note": NOTE_SYMMETRY}
        if input_note:
            out["input_note"] = input_note
        return out

    # ---- the network --------------------------------------------------------------------------
    def network(self, gene, filters=None, top_n=25, pin=()):
        if isinstance(gene, (list, tuple)) and len(gene) > 1:
            raise ValueError("network takes one gene, got %d; call it once per gene" % len(gene))
        gene = _names(gene, "gene", GENE_EXAMPLE, one=True)[0]
        filters = check_filters(filters)
        pin = [pin] if isinstance(pin, str) else list(pin or [])
        if top_n is not None and (isinstance(top_n, bool) or not isinstance(top_n, (int, float)) or top_n < 0):
            raise ValueError("top_n is a whole number, or null to mark every partner as drawn")
        r = self._resolve(gene)
        if r["gene"] is None:
            _not_symbol(gene)
        note = "an edge aggregates every reported site pair between two proteins on one side of the pair; the statistic is symmetric and licenses no direction, and the number on an edge is the strongest of its site pairs, not a consensus"
        st = self.status(r["gene"])[0]
        if (st.get("n_served_as_glyco", 0) + st.get("n_served_as_phospho", 0)) == 0:
            return {"gene": r["gene"], "resolved_via": r["via"], "edges": [], "withheld": None, "n_rows": 0, "n_cis_edges": 0, "top_n": top_n, "note": note,
                    "no_pair": no_pair_note(st), "release": self.release().get("release"), "filters": filters}
        sh = self._shard(r["gene"])
        res = self.network_edges(sh["rows"], r["gene"], filters, top_n, pin)
        res.update({"inference": sh["inference"], "gene": r["gene"], "resolved_via": r["via"], "n_edges": len(res["edges"]), "n_partner_genes": len(set(e["partner"] for e in res["edges"])),
                    "n_both_roles": len(res["edges"]) - len(set(e["partner"] for e in res["edges"])), "top_n": top_n, "note": note, "release": self.release().get("release"), "filters": filters})
        return res

    def network_edges(self, rows0, gene, filters=None, top_n=25, pin=()):
        rows, withheld = self._filter(rows0, filters)
        n_rows = len(rows)
        cis = [r for r in rows if r.get("is_cis") is True]
        rows = [r for r in rows if r.get("is_cis") is not True]
        groups = {}
        for r in rows:
            as_glyco = r["glyco_gene"] == gene
            role = "glyco" if as_glyco else "phospho"
            partner = r["phospho_gene"] if as_glyco else r["glyco_gene"]
            groups.setdefault((partner, role), []).append(r)
        edges = []
        for (partner, role), rs in groups.items():
            b = sorted(rs, key=lambda r: (-r["evidence"], r["glyco_site"], r["phospho_site"]))[0]
            own, other = ("glyco_site", "phospho_site") if role == "glyco" else ("phospho_site", "glyco_site")
            pos = sum(1 for r in rs if not _is_na(r.get("rho_combined")) and r["rho_combined"] >= 0)
            neg = sum(1 for r in rs if not _is_na(r.get("rho_combined")) and r["rho_combined"] < 0)
            n_tot = self.n_total(b.get("cohorts_rho")) if "cohorts_rho" in b else None
            edges.append({
                "gene": gene, "partner": partner, "role": role, "rank": 0, "drawn": False,
                "n_site_pairs": len(rs), "n_own_sites": len(set(r[own] for r in rs)), "n_partner_sites": len(set(r[other] for r in rs)),
                "n_site_pairs_positive": pos, "n_site_pairs_negative": neg, "sign_mixed": pos > 0 and neg > 0,
                "evidence_best": b["evidence"], "evidence_median": median(r["evidence"] for r in rs), "pctile_at_best": b.get("evidence_pctile_k"),
                "rho_at_best": b.get("rho_combined"), "k_at_best": b.get("k_cohorts"), "n_at_best": n_tot, "I2_at_best": b.get("I2"),
                "q_at_best": b.get("q_value"), "q_anchored_at_best": b.get("q_value_anchored"),
                "rho_prot": b.get("rho_prot"), "shadow_clean_at_best": b.get("shadow_clean"), "self_coupling_at_best": b.get("self_coupling"),
                "self_glyco_at_best": b.get("self_glyco"), "self_phospho_at_best": b.get("self_phospho"),
                "n_site_pairs_anchored": sum(1 for r in rs if r.get("anch_any") is True), "anch_any": any(r.get("anch_any") is True for r in rs),
                "anch_dose": any(r.get("anch_dose") is True for r in rs), "anch_gene_editing": any(r.get("anch_gene_editing") is True for r in rs), "anch_enzymatic": any(r.get("anch_enzymatic") is True for r in rs),
                "own_site_at_best": b[own], "best_glyco_site": b["glyco_site"], "best_phospho_site": b["phospho_site"],
            })
        edges.sort(key=lambda e: (-e["evidence_best"], e["partner"], 0 if e["role"] == "glyco" else 1))
        for i, e in enumerate(edges, 1):
            e["rank"] = i
            e["drawn"] = True if top_n is None else (i <= top_n or e["partner"] in pin)
        return {"edges": edges, "withheld": withheld, "n_rows": n_rows, "n_cis_edges": len(cis)}

    # ---- the reading rules --------------------------------------------------------------------
    def reference(self, **kw):
        """This reference as data; `as`="json" returns it as JSON text (the R package's ptm_reference(as =))."""
        if self._reference is None:
            with open(_reference_path(), encoding="utf-8") as fh:
                self._reference = json.load(fh)
        fmt = kw.pop("as", "list")
        if kw:
            raise TypeError("unexpected arguments: " + ", ".join(sorted(kw)))
        if fmt == "json":
            return json.dumps(self._reference)
        if fmt != "list":
            raise ValueError("as must be list or json")
        return self._reference

    def columns(self):
        return self.reference()["columns"]

    @staticmethod
    def _fmt(x, d=3):
        if _is_na(x):
            return "not measured"
        if isinstance(x, bool):
            return str(x).lower()
        if x != 0 and abs(x) < 1e-3:
            return "%.1e" % x
        return "%.*f" % (d, x)

    def _rule_value(self, rule, r, glyco_ambiguous=None):
        g = r.get
        rid = rule["id"]
        if rid == "direction":
            return None
        if rid == "significance":
            return g("q_value")
        if rid == "perturbation":
            return g("anch_any")
        if rid == "same_protein":
            return g("is_cis")
        if rid == "protein_shadow":
            return g("rho_prot")
        if rid == "self_coupling":
            return g("self_coupling")
        if rid == "tissue_composition":
            return self._glyco_labelled(g("glyco_programs"), g("rep_glycans")) and bool(g("phospho_programs") or "")
        if rid == "depth":
            return g("evidence_pctile_k")
        if rid == "heterogeneity":
            return g("I2")
        if rid == "one_cohort":
            return g("k_cohorts")
        if rid == "without_one_cohort":
            return g("loco_pctile_k")
        if rid == "adjusted":
            lo, hi = g("rho_adjusted_lo"), g("rho_adjusted_hi")
            return None if _is_na(lo) or _is_na(hi) else (lo > 0 or hi < 0)
        if rid == "glycoform_pooling":
            return g("n_rep_glycans")
        if rid == "sequon_ambiguity":
            v = g("glyco_site_ambiguous")
            return glyco_ambiguous if _is_na(v) else v
        return None

    @staticmethod
    def _rule_applies(rule, value):
        c = rule.get("comparator")
        if c is None or c == "always":
            return True
        if c == "informational":
            return None
        if _is_na(value):
            return None
        t = rule.get("threshold")
        return {"<": lambda: value < t, "<=": lambda: value <= t, ">": lambda: value > t, ">=": lambda: value >= t,
                "abs>=": lambda: abs(value) >= t, "==": lambda: value == t, "is_true": lambda: value is True}.get(c, lambda: None)()

    def _rule_statement(self, rule, value, applies, r):
        g, fmt = r.get, self._fmt
        rid, extra = rule["id"], ""
        # a rule that tests a value says so when the pair has none: its sentence would otherwise read as a finding
        if applies is None and rule.get("comparator") not in (None, "always", "informational"):
            return "Not assessed: this pair has no value for %s." % (", ".join(rule["field"]) if isinstance(rule.get("field"), list) else rule.get("field"))

        def cohorts(k):
            return "an unstated number of cohorts" if _is_na(k) else "1 cohort" if k == 1 else "%s cohorts" % fmt(k, 0)
        if rid == "perturbation" and applies is True:
            ch = []
            if g("anch_dose") is True:
                ch.append("dose-response (the pair was measured together)")
            if g("anch_gene_editing") is True:
                ch.append("gene editing (each site on its own)")
            if g("anch_enzymatic") is True:
                ch.append("sialidase (each site on its own)")
            extra = " Channels: " + "; ".join(ch) + "."
        elif rid == "self_coupling" and applies is True:
            # the side with the larger absolute correlation with its own protein, and that correlation with its sign
            sg, sp = g("self_glyco"), g("self_phospho")
            known = not _is_na(sg) and not _is_na(sp)
            side = ("glycosite" if abs(sg) >= abs(sp) else "phosphosite") if known else "stronger side"
            own = (sg if abs(sg) >= abs(sp) else sp) if known else value
            extra = (" The %s runs against its protein, at %s." if own < 0 else " The %s tracks its protein at %s.") % (side, fmt(own, 2))
        elif rid == "protein_shadow" and not _is_na(value):
            extra = " Protein rho %s." % fmt(value, 2)
        elif rid == "depth" and not _is_na(value):
            # evidence_pctile_k ranks among ALL scored site pairs with the same k, not only reported ones; the share from
            # the top is floored at the column's resolution (0.1 %) so a top pair never reads "top 0.0 %", as on the site
            extra = " This pair is in the top %s%% of all scored site pairs measured in %s." % (fmt(max(0.1, 100 * (1 - value)), 1), cohorts(g("k_cohorts")))
        elif rid == "heterogeneity" and not _is_na(value):
            extra = " I2 = %s over %s." % (fmt(value, 2), cohorts(g("k_cohorts")))
        elif rid == "without_one_cohort" and not _is_na(value):
            who = "" if _is_na(g("cohort_max_influence")) else " (%s)" % cohort_label(g("cohort_max_influence"))
            extra = " Without its most influential cohort%s, the pair keeps a larger share of its evidence than %s%% of reported pairs measured in %s." % (who, fmt(100 * value, 0), cohorts(g("k_cohorts")))
        elif rid == "adjusted" and not _is_na(value):
            extra = " Adjusted interval [%s, %s]." % (fmt(g("rho_adjusted_lo"), 2), fmt(g("rho_adjusted_hi"), 2))
        elif rid == "glycoform_pooling" and not _is_na(value):
            extra = " %s %s across the cohorts." % (fmt(value, 0), "representative glycoform" if value == 1 else "distinct representative glycoforms")
        elif rid == "significance" and not _is_na(value):
            extra = " q = %s; evidence %s over %s and %s tumors." % (fmt(value), fmt(g("evidence"), 1), cohorts(g("k_cohorts")), "an unstated number of" if _is_na(g("n_total")) else "{:,}".format(int(g("n_total"))))
        return rule["statement"] + extra

    def explain(self, glyco_site=None, phospho_site=None, rows=None):
        rules = self.reference()["rules"]
        if rows is None:
            if glyco_site is None:
                raise ValueError('explain takes a glycosite and a phosphosite (IL6ST_N131 TMX1_S270) or one pair_id ("PTMG:P40189:N131|PTMP:Q9H3N1:S270")')
            if phospho_site is None:
                glyco_site, phospho_site = self._split_pair_id(glyco_site, "explain")
            gl = [glyco_site] if isinstance(glyco_site, str) else list(glyco_site)
            ph = [phospho_site] if isinstance(phospho_site, str) else list(phospho_site)
            if not gl or not ph:
                raise ValueError("explain takes at least one glycosite and one phosphosite")
            if len(gl) != len(ph) and 1 not in (len(gl), len(ph)):
                raise ValueError("glyco_site and phospho_site pair up element by element: give the same number of each, or one of either (got %d and %d)" % (len(gl), len(ph)))
            n = max(len(gl), len(ph))
            rows = []
            for i in range(n):
                p = self.pair(gl[i % len(gl)], ph[i % len(ph)])
                if p.get("input_note"):
                    self._say(p["input_note"])
                if p["row"] is not None:
                    rows.append(p["row"])
                else:   # said, never dropped in silence: the answer then holds the rules of the reported pairs only
                    self._say("%s x %s is not a reported pair of this release, so no rule is read on it" % (gl[i % len(gl)], ph[i % len(ph)]))
        else:
            rows = [dict(r) for r in rows]
            for r in rows:
                if "n_total" not in r and "cohorts_rho" in r:
                    r["n_total"] = self.n_total(r.get("cohorts_rho"))
            self._add_ids(rows)
        out = []
        for r in rows:
            amb = None
            if "glyco_site_ambiguous" not in r:
                sh = self._shard_or_none(r["glyco_gene"])
                amb = next((s.get("glyco_site_ambiguous") for s in sh["glyco_sites"] if s["site"] == r["glyco_site"]), None) if sh else None
            for rule in rules:
                value = self._rule_value(rule, r, amb)
                applies = self._rule_applies(rule, value)
                out.append({"glyco_site": r["glyco_site"], "phospho_site": r["phospho_site"], "pair_id": r.get("pair_id"), "rule": rule["id"],
                            "field": ",".join(rule["field"]) if isinstance(rule.get("field"), list) else rule.get("field"),
                            "value": None if _is_na(value) else (str(value).upper() if isinstance(value, bool) else str(value)),
                            "threshold": None if rule.get("threshold") is None else "%s %s" % (rule["comparator"], str(rule["threshold"]).upper() if isinstance(rule["threshold"], bool) else rule["threshold"]),
                            "applies": applies, "statement": self._rule_statement(rule, value, applies, r), "filter": rule.get("filter")})
        return out

    # ---- sites and their modeled position ----------------------------------------------------
    def sites(self, x, structure=True):
        xs = _names(x, "x", GENE_EXAMPLE)
        out = []
        for q in xs:
            try:
                r = self._resolve(q)
            except AmbiguousName as e:
                if len(xs) == 1:
                    raise
                self._say("%s skipped: %s" % (q, str(e).split(": ", 1)[1]))
                continue
            if r["gene"] is None:
                self._say(("%s is a site key, not a gene; the site command lists its partners" if re.match(r"^.+_[NSTY][0-9]+$", q.upper())
                           else "%s is not a recognized symbol, alias or UniProt accession") % q)
                continue
            sh = self._shard_or_none(r["gene"])
            if sh is None:
                self._say("%s has no reported pair, so it has no reported site; status %s says what was measured" % (r["gene"], r["gene"]))
                continue
            recs = []
            for s in sh.get("glyco_sites") or []:
                recs.append({"query": q, "gene": r["gene"], "uniprot_acc": s.get("uniprot_acc"), "side": "glyco", "site": s["site"], "row_id": s.get("row_id"), "residue": "N",
                             "position_canonical": s.get("position"), "residue_ok": s.get("residue_ok"), "n_cohorts_served": s.get("n_cohorts_served"), "self_coupling_median": s.get("self_coupling_median"),
                             "window_15": s.get("window_15"), "frame_offset": None, "topology_class": s.get("topology_class"), "glyco_site_ambiguous": s.get("glyco_site_ambiguous")})
            for s in sh.get("phospho_sites") or []:
                pg, pc = s.get("position_gencode34"), s.get("position_canonical")
                recs.append({"query": q, "gene": r["gene"], "uniprot_acc": s.get("uniprot_acc"), "side": "phospho", "site": s["site"], "row_id": s.get("row_id"),
                             "residue": s.get("residue_canonical") if not _is_na(s.get("residue_canonical")) else s.get("residue_shipped"), "position_canonical": pc, "residue_ok": s.get("residue_ok"),
                             "n_cohorts_served": s.get("n_cohorts_served"), "self_coupling_median": s.get("self_coupling_median"), "window_15": s.get("window_15"),
                             "frame_offset": None if _is_na(pg) or _is_na(pc) else pg - pc, "topology_class": None, "glyco_site_ambiguous": None})
            # along the protein: glycosites, then phosphosites, each by position (the files list S1026 before S991)
            recs.sort(key=lambda s: (s["side"] != "glyco", s["position_canonical"] if s["position_canonical"] is not None else math.inf, s["site"].encode()))
            # the model: found as the website finds it; a failure to read it is said, never hidden, and never fails the call
            st, status = None, "not requested"
            if structure:
                try:
                    st = self._structure(r["gene"])
                    status = "model read" if st is not None else "no AlphaFold model in this data source"
                    if st is None:
                        self._say("%s: no AlphaFold model in this data source, so its structure columns are NA" % r["gene"])
                except Exception as e:   # noqa: BLE001 -- the site list stands without the model
                    st, status = None, "could not be read: %s" % e
                    self._say("the structure model of %s could not be read, so its structure columns are NA: %s" % (r["gene"], e))
            for rec in recs:
                if st is None:
                    rec.update({"model": None, "model_coverage": None, "modelled": None, "model_residue": None, "residue_matches_model": None, "plddt": None, "x": None, "y": None, "z": None,
                                "structure_status": status})
                    continue
                pos = rec["position_canonical"]
                inside = pos is not None and st["first"] <= pos <= st["last"]
                idx = pos - st["first"] if inside else None
                rec["model"] = "%s v%s" % (st["entity"], st["v"]) if st.get("entity") else None
                rec["model_coverage"] = st["n"] / st["seq_len"] if st.get("seq_len") else None
                rec["modelled"] = inside
                rec["model_residue"] = st["seq"][idx] if inside else None
                rec["residue_matches_model"] = (rec["model_residue"] == rec["residue"]) if inside else None
                rec["plddt"] = st["plddt"][idx] if inside else None
                rec["x"], rec["y"], rec["z"] = (st["xyz"][idx] if inside else (None, None, None))
                rec["structure_status"] = status
            out.extend(recs)
        return out

    # ---- the reported pairs ----------------------------------------------------------------------
    OPS = ("<", "<=", ">", ">=", "==", "!=", "in", "not_in", "is_na", "not_na")

    def _dataset(self):
        if not self.release_dir:
            raise PTMlinkError("query and list read the full table of reported pairs (site_served.parquet), which is released with "
                               "the published article; once you have it, give its release directory with --release-dir (release_dir= in Python).")
        try:
            import pyarrow.dataset as pads
        except ImportError as e:
            raise PTMlinkError("query and list need pyarrow: pip install pyarrow") from e
        pq = os.path.join(self.release_dir, "site_served.parquet")
        if not os.path.exists(pq):
            raise FileNotFoundError("no site_served.parquet in " + self.release_dir)
        return pads.dataset(pq)

    # the columns the default view, the identifiers and n_total need; always read whatever the projection asks
    NEED = ["glyco_site", "phospho_site", "glyco_gene", "phospho_gene", "evidence", "is_cis", "k_cohorts", "q_value", "I2", "anch_any", "shadow_clean", "self_coupling",
            "glyco_programs", "phospho_programs", "rep_glycans", "cohorts_rho", "rho_adjusted_lo", "rho_adjusted_hi", "loco_pctile_k", "glyco_row_id", "phospho_row_id", "pair_id"]

    def _provenance(self):
        versions = {}
        vp = os.path.join(self.release_dir, "versions.tsv")
        if os.path.exists(vp):
            with open(vp, encoding="utf-8") as fh:
                fh.readline()
                for line in fh:
                    cells = line.rstrip("\n").split("\t")
                    if len(cells) >= 2:
                        versions[cells[0]] = cells[1]
        return {"inference": self._inference(versions), "release": public_label(versions.get("release"), versions.get("public_release")), "built_utc": versions.get("built_utc"), "file": os.path.join(self.release_dir, "site_served.parquet"), "sha256": self.release_sha}

    # one scan of the reported pairs: the filter pushed down, the default view after it, the one sort key, provenance
    # on the result; query() and list() both come through here, as in the R package
    def _scan(self, ds, filt, columns, filters, limit, asked):
        schema = ds.schema
        if columns is None:
            proj = schema.names
        else:
            bad = [c for c in columns if c not in schema.names]
            if bad:
                raise ValueError("no column(s) %s; columns() lists them" % ", ".join(bad))
            proj = list(dict.fromkeys(list(columns) + [c for c in self.NEED if c in schema.names]))
        tb = ds.to_table(filter=filt, columns=proj)
        rows = tb.to_pylist()
        n_matched = len(rows)
        filters = check_filters(filters)
        out, withheld = self._filter(rows, filters)
        n_view = len(out)
        if limit is not None:
            out = out[:limit]
        if columns is not None:
            keep = list(dict.fromkeys(list(columns) + ["glyco_row_id", "phospho_row_id", "pair_id", "n_total"]))
            out = [{k: r.get(k) for k in keep if k in r} for r in out]
        import datetime
        prov = self._provenance()
        prov.update(asked)
        prov.update({"filters": filters, "limit": limit, "sorted_on": "evidence (descending); the total score is the only ranking",
                     "n_matched": n_matched, "n_in_view": n_view, "n_returned": len(out), "licence": "data CC BY 4.0", "queried_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                     "note": "the score is symmetric and licenses no direction; adjusted, parent and self-coupling columns are annotations, never a ranking"})
        return {"rows": out, "withheld": withheld, "n_matched": n_matched, "n_view": n_view, "filters": filters, "provenance": prov}

    # A query whose conditions match more rows than this is refused unless allow_full is set: every matching row is
    # read into memory before the view and the limit apply (3.4 million rows took minutes and tens of GB).
    QUERY_MAX_ROWS = 1000000

    @staticmethod
    def _condition_value(col, ty, v, text=False):
        """One value of a condition, as its column's type; `text` reads a number or true/false typed on the command line."""
        shown = json.dumps(v) if not isinstance(v, float) or v == v else "NaN"
        if v is None or (isinstance(v, float) and v != v):
            raise ValueError('%s: a condition compares with a value, got none; use the operator "is_na" or "not_na" to test for a missing value' % col)
        if ty.startswith(("int", "uint")) or ty in ("double", "float"):
            if text and isinstance(v, str):
                try:
                    v = float(v)
                except ValueError:
                    raise ValueError("%s holds numbers; got %s" % (col, shown)) from None
            if isinstance(v, bool) or not isinstance(v, (int, float)):
                raise ValueError("%s holds numbers; got %s" % (col, shown))
            if ty in ("double", "float"):
                return float(v)
            if v != int(v):
                raise ValueError("%s holds whole numbers; got %s" % (col, shown))
            return int(v)
        if ty == "bool":
            if text and isinstance(v, str) and v.lower() in ("true", "false"):
                v = v.lower() == "true"
            if not isinstance(v, bool):
                raise ValueError("%s holds true or false; got %s" % (col, shown))
            return v
        if not isinstance(v, str):
            raise ValueError("%s holds text; got %s" % (col, shown))
        return v

    def _query_genes(self, genes):
        """The genes of a query as the release spells them: an alias or accession resolved as the per-gene readers
        resolve it, a name the index does not know said and kept as typed (it then matches no pair)."""
        out = []
        for g in genes:
            try:
                r = self._resolve(g)
            except AmbiguousName:
                raise
            except Exception:   # no index to read: the name as typed, in the release's letter case
                out.append(self._canon_gene(g))
                continue
            if r["gene"] is None:
                self._say("%s is not a recognized symbol, alias or UniProt accession: it matches no pair" % g)
                out.append(self._canon_gene(g))
            else:
                out.append(r["gene"])
        return list(dict.fromkeys(out))

    def query(self, where=(), genes=None, sites=None, columns=None, filters=None, limit=None, allow_full=False, _text_values=False):
        filters = check_filters(filters)
        if not isinstance(where, (list, tuple)):
            raise ValueError('where is a list of conditions, each [column, operator, value], such as [["q_value", "<", 1e-10]]')
        genes = None if genes is None else _names(genes, "genes", GENE_EXAMPLE)
        sites = None if sites is None else _names(sites, "sites", SITE_EXAMPLE)
        columns = None if columns is None else _names(columns, "columns", "column names; the columns command lists them")
        if limit is not None and (isinstance(limit, bool) or not isinstance(limit, int) or limit < 0):
            raise ValueError("limit is a whole number, 0 or more")
        ds = self._dataset()
        import pyarrow.compute as pc
        schema = ds.schema
        where = list(where)
        if where and not isinstance(where[0], (list, tuple)):
            where = [where]
        if not where and genes is None and sites is None and not allow_full:
            raise ValueError("no condition: give where, genes or sites, or set allow_full to read every reported pair")
        if any(not isinstance(c, (list, tuple)) or not 2 <= len(c) <= 3 or not isinstance(c[0], str) or not isinstance(c[1], str) for c in where):
            raise ValueError('each condition is [column, operator, value], such as ["k_cohorts", ">=", 5]')

        def cond(c):
            col, op = str(c[0]), str(c[1])
            if col not in schema.names:
                raise ValueError("no column %s in site_served.parquet; columns() lists them" % col)
            if op not in self.OPS:
                raise ValueError("operator %s is not one of %s" % (op, ", ".join(self.OPS)))
            f = pc.field(col)
            if op == "is_na":
                return f.is_null()
            if op == "not_na":
                return ~f.is_null()
            if len(c) < 3:
                raise ValueError("operator %s needs a value" % op)
            ty = str(schema.field(col).type)
            v = c[2]
            if op in ("in", "not_in"):
                vals = [self._condition_value(col, ty, x, _text_values) for x in (v if isinstance(v, (list, tuple)) else [v])]
                if not vals:
                    raise ValueError("operator %s needs at least one value" % op)
                e = f.isin(vals)
                return e if op == "in" else ~e
            if isinstance(v, (list, tuple)):
                raise ValueError('operator %s takes one value; use "in" for a set' % op)
            v = self._condition_value(col, ty, v, _text_values)
            return {"<": f < v, "<=": f <= v, ">": f > v, ">=": f >= v, "==": f == v, "!=": f != v}[op]

        exprs = [cond(c) for c in where]
        if genes is not None:
            gs = self._query_genes(genes)
            exprs.append(pc.field("glyco_gene").isin(gs) | pc.field("phospho_gene").isin(gs))
        if sites is not None:
            ss = [self._canon_site(s) for s in sites]
            exprs.append(pc.field("glyco_site").isin(ss) | pc.field("phospho_site").isin(ss))
        filt = None
        for e in exprs:
            filt = e if filt is None else (filt & e)
        if not allow_full and filt is not None:
            n = ds.count_rows(filter=filt)
            if n > self.QUERY_MAX_ROWS:
                raise ValueError("the conditions match {:,} reported pairs, more than {:,}: reading them all takes minutes and tens of GB of memory "
                                 "(limit shortens the answer, not the read). Add a condition, or set allow_full to read them anyway".format(n, self.QUERY_MAX_ROWS))
        return self._scan(ds, filt, columns, filters, limit, {"where": [list(c) for c in where], "genes": genes, "sites": sites})

    # ---- a list of proteins or sites against the tier (ptm_list) ----------------------------------
    # The concepts, defined once in ptmlink/R/list.R: level, member, internal pair, neighborhood, partner site,
    # convergent partner, connectivity. The row tables obey the view; the tests count the significant reported pairs
    # between two different proteins.
    @classmethod
    def _entries(cls, x):
        """The entries of a list: one per element, or several separated by newlines, commas, semicolons or
        whitespace. Every token is an entry; a repeated entry is kept once, where it first appears."""
        xs = [x] if isinstance(x, str) else list(x)
        out, seen = [], set()
        for el in xs:
            for q in re.split(r"[\n\r,;\s]+", "" if el is None else str(el)):
                if q and q not in seen:
                    seen.add(q)
                    out.append({"query": q})
        return out

    TEST_UNIVERSE = "significant reported pairs (q_value < 0.05) between two different proteins"

    def _tier(self):
        """The pairs the list tests count, as the tests need them: every SIGNIFICANT reported row (q_value < 0.05) between two
        different proteins, its two sites and two proteins as integer codes, and, per unit on each side, its degree and
        reach, under the names ptmlink/R/list.R::.ptm_tier uses. A pair reported only for its perturbation support is not
        counted."""
        key = (self.release_dir, self.release_sha)
        if getattr(self, "_tier_cache", None) and self._tier_cache["key"] == key:
            return self._tier_cache
        import pyarrow as pa
        import pyarrow.compute as pc
        tb = self._dataset().to_table(columns=["glyco_site", "phospho_site", "glyco_gene", "phospho_gene"], filter=pc.field("q_value") < 0.05)
        tb = tb.filter(pc.not_equal(tb.column("glyco_gene"), tb.column("phospho_gene")))

        def encode(col):
            d = pc.dictionary_encode(tb.column(col)).combine_chunks()
            names = d.dictionary.to_pylist()
            order = sorted(range(len(names)), key=lambda i: names[i].encode())   # sorted levels, as R's factor gives them
            rank = [0] * len(names)
            for r, i in enumerate(order):
                rank[i] = r
            idx = pc.take(pa.array(rank, pa.int64()), d.indices)
            return [names[i] for i in order], idx
        levels, codes = {}, {}
        for col, short in (("glyco_site", "gs"), ("phospho_site", "ps"), ("glyco_gene", "gg"), ("phospho_gene", "pg")):
            levels[short], codes[short] = encode(col)

        def counts(idx, n):
            vc = pc.value_counts(idx)
            out = [0] * n
            for v, c in zip(vc.field("values").to_pylist(), vc.field("counts").to_pylist()):
                out[v] = c
            return out

        def distinct(a, b):
            nb = len(levels[b])
            u = pc.unique(pc.add(pc.multiply(codes[a], nb), codes[b]))
            ua = pc.divide(u, nb)
            return ua, pc.subtract(u, pc.multiply(ua, nb))
        d_gg_pg = distinct("gg", "pg"); d_ps_gg = distinct("ps", "gg"); d_gs_pg = distinct("gs", "pg")
        n = {k: len(v) for k, v in levels.items()}
        first_gene = {}
        for side, scol, gcol in (("glyco", "gs", "gg"), ("phospho", "ps", "pg")):
            sg = [None] * n[scol]
            s_list, g_list = codes[scol].to_pylist(), codes[gcol].to_pylist()
            for s, g in zip(s_list, g_list):
                if sg[s] is None:
                    sg[s] = levels[gcol][g]
            first_gene[side] = sg
        tier = {"key": key, "M": tb.num_rows, "rows": codes, "levels": levels,
                "sites_glyco": levels["gs"], "sites_phospho": levels["ps"], "genes_glyco": levels["gg"], "genes_phospho": levels["pg"],
                "site_gene_glyco": first_gene["glyco"], "site_gene_phospho": first_gene["phospho"],
                "gene_glyco": {"deg": counts(codes["gg"], n["gg"]), "reach_genes": counts(d_gg_pg[0], n["gg"]), "reach_sites": counts(d_ps_gg[1], n["gg"])},
                "gene_phospho": {"deg": counts(codes["pg"], n["pg"]), "reach_genes": counts(d_gg_pg[1], n["pg"]), "reach_sites": counts(d_gs_pg[1], n["pg"])},
                "site_glyco": {"deg": counts(codes["gs"], n["gs"]), "reach_genes": counts(d_gs_pg[0], n["gs"])},
                "site_phospho": {"deg": counts(codes["ps"], n["ps"]), "reach_genes": counts(d_ps_gg[0], n["ps"])},
                "index": {k: {name: i for i, name in enumerate(v)} for k, v in levels.items()}}
        self._tier_cache = tier
        return tier

    @staticmethod
    def _tier_key(level, side):
        return ("gs" if side == "glyco" else "ps") if level == "site" else ("gg" if side == "glyco" else "pg")

    def _codes(self, tier, level, side, units):
        ix = tier["index"][self._tier_key(level, side)]
        return [ix.get(u) for u in units]

    def _deg(self, tier, level, side, units):
        deg = tier[("site_" if level == "site" else "gene_") + side]["deg"]
        return [0 if c is None else deg[c] for c in self._codes(tier, level, side, units)]

    @staticmethod
    def _converge(pu, mr, reach, pnames, N, n, mnames):
        """One family: partner units reached by the members, with the members reached and the pair count."""
        if not pu or not n:
            return []
        seen, npairs, reached = {}, {}, {}
        for p, m in zip(pu, mr):
            npairs[p] = npairs.get(p, 0) + 1
            if (p, m) not in seen:
                seen[(p, m)] = True
                reached.setdefault(p, []).append(mnames[m])
        return [{"partner": pnames[p], "members_reached": len(ms), "members_testable": n, "tier_reach": reach[p], "tier_n": N, "n_pairs": npairs[p], "members": ";".join(sorted(ms))}
                for p, ms in reached.items()]

    def _list_tests(self, tier, level, units, neighbourhood):
        import pyarrow as pa
        import pyarrow.compute as pc
        r = tier["rows"]

        def one_side(member_side):
            partner_side = "phospho" if member_side == "glyco" else "glyco"
            ucol, ocol = self._tier_key(level, member_side), self._tier_key(level, partner_side)
            universe = tier["levels"][ucol]
            codes = self._codes(tier, level, member_side, units)
            present = [u for u, c in zip(units, codes) if c is not None]
            codes = [c for c in codes if c is not None]
            n, N = len(present), len(universe)
            if not n:
                return None, None
            other = [c for c in self._codes(tier, level, partner_side, units) if c is not None]
            keep = pc.is_in(r[ucol], value_set=pa.array(codes, pa.int64()))
            if other:
                keep = pc.and_(keep, pc.invert(pc.is_in(r[ocol], value_set=pa.array(other, pa.int64()))))
            rank = {c: i for i, c in enumerate(codes)}
            mr = [rank[c] for c in pc.filter(r[ucol], keep).to_pylist()]
            site_col, gene_col = ("gs", "gg") if partner_side == "glyco" else ("ps", "pg")
            st, gt = tier["site_" + partner_side], tier["gene_" + partner_side]
            reach_site = st["deg"] if level == "site" else st["reach_genes"]
            reach_gene = gt["reach_sites"] if level == "site" else gt["reach_genes"]
            s = self._converge(pc.filter(r[site_col], keep).to_pylist(), mr, reach_site, tier["levels"][site_col], N, n, present)
            g = self._converge(pc.filter(r[gene_col], keep).to_pylist(), mr, reach_gene, tier["levels"][gene_col], N, n, present)
            site_gene = tier["site_gene_" + partner_side]
            for o in s:
                o["gene"] = site_gene[tier["index"][site_col][o["partner"]]]
                o["side"] = partner_side
            for o in g:
                o["partner_side"] = partner_side
            return s, g
        a = one_side("glyco"); b = one_side("phospho")

        def finish(tab, key):
            side_col = "side" if key == "site" else "partner_side"
            view = self._partner_view(neighbourhood, level, units, key)
            for o in tab:
                if key == "site":
                    o["site"] = o.pop("partner")
                o["fold"] = (o["members_reached"] / o["members_testable"]) / (o["tier_reach"] / o["tier_n"])
                o["p"] = _phyper_upper(o["members_reached"], o["tier_reach"], o["tier_n"], o["members_testable"])
            for o, q in zip(tab, _bh([o["p"] for o in tab])):
                o["q"] = q
                v = view.get((o[key], o[side_col]))
                o["evidence_best"] = v["evidence_best"] if v else None
                o["best_pair_id"] = v["best_pair_id"] if v else None
            tab.sort(key=lambda o: (o["p"], -o["fold"], -(o["evidence_best"] if o["evidence_best"] is not None else -math.inf), o[key].encode()))
            cols = (["site", "gene", "side"] if key == "site" else ["partner", "partner_side"]) + ["members_reached", "members_testable", "tier_reach", "tier_n", "fold", "p", "q", "n_pairs", "members", "evidence_best", "best_pair_id"]
            return [{c: o[c] for c in cols} for o in tab]
        return {"sites": finish((a[0] or []) + (b[0] or []), "site"), "proteins": finish((a[1] or []) + (b[1] or []), "partner")}

    @staticmethod
    def _partner_view(nb, level, units, key):
        """What the view says about each partner: the strongest row (rows are sorted on evidence, so the first is the best)."""
        us = set(units)
        col = "site" if key == "site" else "gene"
        out = {}
        for r in nb:
            on_glyco = (r["glyco_site"] if level == "site" else r["glyco_gene"]) in us
            partner = r[("phospho_" if on_glyco else "glyco_") + col]
            side = "phospho" if on_glyco else "glyco"
            k = (partner, side)
            if k not in out:
                out[k] = {"evidence_best": r["evidence"], "best_pair_id": r.get("pair_id")}
        return out

    def _connectivity(self, tier, level, units, n_perm=200, seed=1):
        import random
        import pyarrow as pa
        import pyarrow.compute as pc
        r = tier["rows"]
        gcol, pcol = self._tier_key(level, "glyco"), self._tier_key(level, "phospho")

        def count(s):
            cg = [c for c in self._codes(tier, level, "glyco", s) if c is not None]
            cp = [c for c in self._codes(tier, level, "phospho", s) if c is not None]
            if not cg or not cp:
                return 0
            m = pc.and_(pc.is_in(r[gcol], value_set=pa.array(cg, pa.int64())), pc.is_in(r[pcol], value_set=pa.array(cp, pa.int64())))
            return pc.sum(m).as_py() or 0

        def expected(s):
            a = self._deg(tier, level, "glyco", s)
            b = self._deg(tier, level, "phospho", s)
            by_gene = {}
            for u, x, y in zip(s, a, b):
                gene = u.rsplit("_", 1)[0] if level == "site" else u
                d = by_gene.setdefault(gene, [0, 0])
                d[0] += x
                d[1] += y
            return (sum(a) * sum(b) - sum(x * y for x, y in by_gene.values())) / tier["M"] if tier["M"] else None
        kg, kp = self._deg(tier, level, "glyco", units), self._deg(tier, level, "phospho", units)
        present = [u for u, a, b in zip(units, kg, kp) if a + b > 0]
        obs, exp_ = count(units), expected(units)
        out = {"members_in_tier": len(present), "internal_pairs": obs, "expected_pairs": exp_, "ratio": (obs / exp_) if exp_ is not None and exp_ > 0 else None, "p_perm": None, "n_perm": 0, "perm_median": None}
        if n_perm > 0 and len(present) >= 2:
            rng = random.Random(seed)
            all_units = (tier["sites_glyco"] + tier["sites_phospho"]) if level == "site" else list(dict.fromkeys(tier["genes_glyco"] + tier["genes_phospho"]))
            tot = [a + b for a, b in zip(self._deg(tier, level, "glyco", all_units), self._deg(tier, level, "phospho", all_units))]
            # deciles as R's quantile (type 7) and cut(include.lowest = TRUE)
            srt = sorted(tot)
            def q7(p):
                h = (len(srt) - 1) * p
                lo = int(math.floor(h))
                return srt[lo] + (h - lo) * (srt[min(lo + 1, len(srt) - 1)] - srt[lo])
            br = sorted(set(q7(i / 10) for i in range(11)))
            def bin_of(v):
                if len(br) < 2:
                    return 0
                for i in range(1, len(br)):
                    if v <= br[i]:
                        return i - 1
                return len(br) - 2
            bins = {}
            for u, v in zip(all_units, tot):
                bins.setdefault(bin_of(v), []).append(u)
            need = {}
            tot_of = dict(zip(all_units, tot))
            for u in present:
                b = bin_of(tot_of[u])
                need[b] = need.get(b, 0) + 1
            perm = []
            for _ in range(n_perm):
                draw = []
                for b, k in need.items():
                    pool = bins[b]
                    draw.extend(pool if len(pool) <= k else rng.sample(pool, k))
                perm.append(count(draw))
            out["p_perm"] = (1 + sum(1 for v in perm if v >= obs)) / (1 + n_perm)
            out["n_perm"] = n_perm
            out["perm_median"] = median(perm)
        return out

    def list(self, x, level="auto", filters=None, n_perm=200, seed=1):
        if level not in ("auto", "protein", "site"):
            raise ValueError("level is auto, protein or site")
        filters = check_filters(filters)
        if isinstance(n_perm, bool) or not isinstance(n_perm, int) or n_perm < 0:
            raise ValueError("n_perm is a whole number of random lists, 0 or more")
        if not isinstance(x, (str, list, tuple)):
            raise ValueError("x is a list of gene symbols or site keys")
        ds = self._dataset()
        import pyarrow.compute as pc
        ent = self._entries(x)
        if not ent:
            raise ValueError("the list is empty: give gene symbols or site keys")
        is_site = [re.match(r"^.+_[NSTY][0-9]+$", e["query"].upper()) is not None for e in ent]
        mixed = level == "auto" and any(is_site) and not all(is_site)
        if level == "auto":
            level = "site" if any(is_site) else "protein"
        if level == "site" and not all(is_site) and not mixed:
            bad = [e["query"] for e, s in zip(ent, is_site) if not s]
            raise ValueError("at site level every entry must be a site key (GENE_N123 or GENE_S123); %s%s are not. Pass level = 'protein' to use their proteins" % (", ".join(bad[:5]), ", ..." if len(bad) > 5 else ""))
        members = []
        for e, s in zip(ent, is_site):
            q = e["query"]
            k = self._site_key(q) if s else None
            try:
                r = self._resolve(k["gene"] if s else q)
            except ValueError as err:
                r = {"gene": None, "via": str(err)}
            members.append({"query": q, "kind": "site" if s else "gene", "gene": r["gene"],
                            "resolved_via": "not a recognized symbol" if r["gene"] is None and r["via"] == "unknown" else r["via"],
                            "site": (r["gene"] + "_" + k["key"].rsplit("_", 1)[1]) if s and r["gene"] else None, "side_named": k["side"] if s else None})
        G = list(dict.fromkeys(m["gene"] for m in members if m["gene"]))
        S = list(dict.fromkeys(m["site"] for m in members if m["site"] and m["gene"]))
        if not G:
            raise ValueError("no entry of the list is a recognized symbol, alias, UniProt accession or site key; status says how names resolve")
        st = {s["gene"]: s for s in self.status(G)}
        for m in members:
            s = st.get(m["gene"]) if m["gene"] else None
            for nm in ("status_glyco", "status_phospho", "n_served_as_glyco", "n_served_as_phospho"):
                m[nm] = s.get(nm) if s else None
        # n_reported: the reported pairs an entry brings, before any filter. A gene entry: its protein's pairs on both
        # sides. A site entry: that site's own pairs; a site in no reported pair is said, as the site command says it.
        for m in members:
            m["n_reported"] = None if not m["gene"] else (m.get("n_served_as_glyco") or 0) + (m.get("n_served_as_phospho") or 0)
        named = [m for m in members if m["kind"] == "site" and m["site"]]
        if named:
            ss = list(dict.fromkeys(m["site"] for m in named))
            own = ds.to_table(columns=["glyco_site", "phospho_site"], filter=pc.field("glyco_site").isin(ss) | pc.field("phospho_site").isin(ss)).to_pydict()
            cnt = {}
            for v in own["glyco_site"] + own["phospho_site"]:
                cnt[v] = cnt.get(v, 0) + 1
            for m in named:
                m["n_reported"] = cnt.get(m["site"], 0)
            none = list(dict.fromkeys(m["site"] for m in named if m["n_reported"] == 0))
            if none:
                self._say("%s %s in no reported pair of this release" % (", ".join(none), "is" if len(none) == 1 else "are"))
        if mixed:
            expanded_genes = {m["gene"] for m in members if m["kind"] == "gene" and m["gene"]}
            expansion = ds.to_table(columns=["glyco_site", "phospho_site", "glyco_gene", "phospho_gene"],
                                    filter=pc.field("glyco_gene").isin(expanded_genes) | pc.field("phospho_gene").isin(expanded_genes)).to_pylist() if expanded_genes else []
            S = list(dict.fromkeys(S + [r[side + "_site"] for r in expansion for side in ("glyco", "phospho") if r[side + "_gene"] in expanded_genes]))
        unit_of = [m["site"] if level == "site" else m["gene"] for m in members]
        units = S if level == "site" else G
        if level == "site":
            filt = pc.field("glyco_site").isin(S) | pc.field("phospho_site").isin(S)
        else:
            filt = pc.field("glyco_gene").isin(G) | pc.field("phospho_gene").isin(G)
        res = self._scan(ds, filt, None, filters, None, {"list": [e["query"] for e in ent], "level": level, "genes": G, "sites": S or None})
        rows = res["rows"]
        Ss, Gs = set(S), set(G)
        for r in rows:
            r["site_named_glyco"] = r["glyco_site"] in Ss
            r["site_named_phospho"] = r["phospho_site"] in Ss
        mg = [(r["site_named_glyco"] if level == "site" else r["glyco_gene"] in Gs) for r in rows]
        mp = [(r["site_named_phospho"] if level == "site" else r["phospho_gene"] in Gs) for r in rows]
        same = [r["glyco_gene"] == r["phospho_gene"] for r in rows]
        internal = [r for r, a, b, s in zip(rows, mg, mp, same) if a and b and not s]
        same_protein = [r for r, a, b, s in zip(rows, mg, mp, same) if s and (a or b)]
        nb = [r for r, a, b, s in zip(rows, mg, mp, same) if (a != b) and not s]

        def in_rows(tab, u):
            return sum(1 for r in tab if u in ((r["glyco_site"], r["phospho_site"]) if level == "site" else (r["glyco_gene"], r["phospho_gene"])))
        for m, u in zip(members, unit_of):
            def count(tab):
                if mixed and m["kind"] == "gene" and m["gene"]:
                    return sum(m["gene"] in (r["glyco_gene"], r["phospho_gene"]) for r in tab)
                return in_rows(tab, u) if u else None
            m["n_internal"] = count(internal)
            m["n_neighbourhood"] = count(nb)
        tier = self._tier()
        tests = self._list_tests(tier, level, units, nb)
        connectivity = self._connectivity(tier, level, units, n_perm, seed)

        def share(tab):
            both = [self._glyco_labelled(r.get("glyco_programs"), r.get("rep_glycans")) and bool(r.get("phospho_programs")) for r in tab]
            return {"rows": len(tab), "both_labelled": sum(both), "share": (sum(both) / len(tab)) if tab else None}
        composition = [dict(table=t, **share(tab)) for t, tab in (("internal", internal), ("same_protein", same_protein), ("neighbourhood", nb))]
        return {"level": level, "mixed": mixed, "sites": S, "members": members, "internal": internal, "same_protein": same_protein, "neighbourhood": nb,
                "partner_sites": tests["sites"], "convergent": tests["proteins"], "connectivity": connectivity,
                "universe": {"rule": self.TEST_UNIVERSE, "pairs": tier["M"]},
                "composition": composition, "withheld": res["withheld"], "provenance": res["provenance"],
                "note": {"rows": "the row tables obey the view and are sorted on evidence, the only ranking; the score is symmetric and licenses no direction",
                         "tests": "the partner-site, convergent-partner and connectivity tests count the significant reported pairs (q < 0.05) between two different proteins; a pair reported only for its perturbation support is not counted. They say a site or protein is reached by more members, or the list is more connected, than the members' own numbers of significant pairs predict, and nothing about mechanism; a partner reached by one member is a pair-level finding, already in the neighbourhood table",
                         "composition": "a list from a tumor-versus-normal comparison is expected to be rich in tissue-composition labels; read the share before the pairs"}}

    # ---- tool descriptions for an agent framework (ptm_tools) -------------------------------------
    JSON_TYPES = {
        # one value or several: an agent may send "EGFR" or ["EGFR", "TFRC"], and both are read
        "character": {"anyOf": [{"type": "string"}, {"type": "array", "items": {"type": "string"}}]},
        "character(1)": {"type": "string"},
        "integer(1)": {"type": "integer"},
        "numeric(1)": {"type": "number"},
        "logical(1)": {"type": "boolean"},
        "list of list(column, operator, value)": {"type": "array", "items": {"type": "array", "minItems": 2, "maxItems": 3,
            "items": {"anyOf": [{"type": "string"}, {"type": "number"}, {"type": "boolean"}, {"type": "null"},
                      {"type": "array", "items": {"anyOf": [{"type": "string"}, {"type": "number"},
                                                           {"type": "boolean"}, {"type": "null"}]}}]}}},
    }

    def tools(self, format="mcp"):
        ref = self.reference()
        conv = ref["tool_conventions"]

        def frag(t):
            if t not in self.JSON_TYPES:
                raise ValueError("no JSON Schema mapping for the argument type " + t)
            return json.loads(json.dumps(self.JSON_TYPES[t]))

        words = conv.get("default_words") or {}

        def shown(v):   # a default as a tool caller writes it (true, null), not as R does (TRUE, NA)
            s = ("TRUE" if v else "FALSE") if isinstance(v, bool) else str(v)
            return words.get(s, s)

        def typed(a):   # an argument that also takes null (no limit, no cap) says so in its schema
            s = frag(a["type"])
            return {"anyOf": [s, {"type": "null"}]} if a.get("nullable") else s
        filter_schema = {"type": "object", "additionalProperties": False,
                         "properties": {f["name"]: dict(typed(f), description="%s (default %s)" % (f["meaning"], shown(f["default"]))) for f in ref["filters"]}}
        tools = []
        for f in ref["functions"]:
            if not f.get("tool"):
                continue
            props, req = {}, []
            for a in f["args"]:
                if a.get("tool") is False:
                    continue
                s = json.loads(json.dumps(filter_schema)) if a["name"] == "filters" else typed(a)
                s["description"] = "%s (default %s)" % (a["meaning"], shown(a["default"])) if a.get("default") is not None else a["meaning"]
                if a.get("enum"):
                    s["enum"] = list(a["enum"])
                props[a["name"]] = s
                if "default" not in a:
                    req.append(a["name"])
            if f.get("rows") and "limit" not in props:
                props["limit"] = {"type": "integer", "description": conv["row_limit"]["meaning"], "default": conv["row_limit"]["default"]}
            schema = {"type": "object", "properties": props, "additionalProperties": False}
            if req:
                schema["required"] = req
            desc = "%s Returns %s.%s" % (f["purpose"], f["returns"]["type"], (" " + conv["rows_note"]) if f.get("ranked") else "")
            if format == "mcp":
                tools.append({"name": f["name"], "description": desc, "inputSchema": schema})
            elif format == "anthropic":
                tools.append({"name": f["name"], "description": desc, "input_schema": schema})
            elif format == "openai":
                tools.append({"type": "function", "function": {"name": f["name"], "description": desc, "parameters": schema}})
            else:
                raise ValueError("format is mcp, anthropic or openai")
        return {"instructions": conv["instructions"], "tools": tools}

    def cite(self):
        p = os.path.join(self.release_dir, "CITATION.cff") if self.release_dir and os.path.exists(os.path.join(self.release_dir, "CITATION.cff")) else os.path.join(os.path.dirname(_reference_path()), "extdata", "ptmlink_citation.cff")
        lines = open(p, encoding="utf-8").read().splitlines()

        def scalar(key):
            for i, l in enumerate(lines):
                if l.startswith(key + ":"):
                    v = l[len(key) + 1:].strip().strip('"')
                    if v in (">-", ""):
                        buf = []
                        j = i + 1
                        while j < len(lines) and lines[j].startswith(" "):
                            buf.append(lines[j].strip())
                            j += 1
                        v = " ".join(buf)
                    return v
            return None
        lic = []
        li = next((i for i, l in enumerate(lines) if l.startswith("license:")), None)
        if li is not None:
            first = lines[li][len("license:"):].strip()
            lic = [first] if first else [re.sub(r"^\s+-\s*", "", l).strip() for l in lines[li + 1:li + 6] if re.match(r"^\s+-", l)]
        # the paper that describes the resource (its preprint, then its article), when the file names one
        paper = None
        p0 = next((i for i, l in enumerate(lines) if l.startswith("preferred-citation:")), None)
        if p0 is not None:
            block = []
            for l in lines[p0 + 1:]:
                if re.match(r"^[a-z-]+:", l):
                    break
                block.append(l)

            def field(key):
                v = next((l for l in block if l.startswith("  " + key + ":")), None)
                return v[len(key) + 3:].strip().strip('"') if v is not None else None
            paper = {k: field(k) for k in ("journal", "year", "doi", "url", "status")}
        c = {"title": scalar("title"), "version": scalar("version"), "date": scalar("date-released"), "license": ", ".join(lic) or None, "doi": scalar("doi"), "url": scalar("url"), "paper": paper}
        # the citation line, as the R package's .ptm_cite_text writes it
        resource = "ersion %s%s%s" % (c["version"], ", released " + c["date"] if c["date"] else "", ", doi:" + c["doi"] if c["doi"] else "")
        if paper and paper.get("doi"):
            c["text"] = "%s. %s%s. doi:%s. PTMlink release v%s%s%s%s. Data CC BY 4.0, code MIT." % (
                c["title"], paper.get("journal") or "Preprint", ", %s" % paper["year"] if paper.get("year") else "", paper["doi"], c["version"],
                ", released " + c["date"] if c["date"] else "", ", doi:" + c["doi"] if c["doi"] else "", ", " + c["url"] if c["url"] else "")
        else:
            c["text"] = "%s%s. V%s%s. Data CC BY 4.0, code MIT." % (c["title"], " (%s)" % c["date"][:4] if c["date"] else "", resource, ". " + c["url"] if c["url"] else "")
        try:
            c["release"] = self.release().get("release")
        except Exception:
            c["release"] = None
        return c


def _phyper_upper(x, K, N, n):
    """P(X >= x) for X hypergeometric: n draws without replacement from N items of which K are marked (R's phyper(x - 1, K, N - K, n, lower.tail = FALSE))."""
    def lchoose(a, b):
        return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)
    lo, hi = max(x, 0, n - (N - K)), min(K, n)
    if lo > hi:
        return 0.0
    denom = lchoose(N, n)
    return min(1.0, sum(math.exp(lchoose(K, i) + lchoose(N - K, n - i) - denom) for i in range(lo, hi + 1)))


def _bh(p):
    """Benjamini-Hochberg over one family, as R's p.adjust(method = 'BH')."""
    n = len(p)
    order = sorted(range(n), key=lambda i: p[i], reverse=True)
    q, cur = [0.0] * n, 1.0
    for rank, i in zip(range(n, 0, -1), order):
        cur = min(cur, p[i] * n / rank)
        q[i] = cur
    return q


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---- text output ------------------------------------------------------------------------------
# The command line prints these by default (--json prints the full result) and the MCP server sends them ahead of
# its JSON, so a person or an agent reads the answer first. They show the key columns only; nothing is recomputed
# except the pooled 95 % interval, +-1.96 (sum of weights)^-1/2 on the Fisher scale, as the website shows it. They
# follow the R package's print methods: a header that says what was asked, what the filters hid, a table of the key
# columns, and the next steps as aligned "label: how" lines. `style` picks the wording of those steps (command-line
# flags, or MCP tool names), and `width` the line length (the terminal's, at the command line).
WITHHELD_WORDS = {"cis": "same-protein", "k": "single-cohort", "composition": "tissue composition",
                  "anchored": "without perturbation support", "significant": "not significant",
                  "heterogeneous": "heterogeneous", "loco": "below the leave-one-cohort-out floor",
                  "adjusted": "adjusted interval includes zero", "shadow": "parent proteins co-vary",
                  "self": "site follows its own protein"}
SUPPORT_WORDS = (("anch_dose", "dose"), ("anch_gene_editing", "editing"), ("anch_enzymatic", "sialidase"))
# a rule's name in plain words, for the "Does not apply" line (as the R package prints it)
RULE_WORDS = {"protein_shadow": "parent-protein co-variation"}
RULE_LEADS = {"direction": "Co-variation only: neither site is shown to regulate the other.",
              "significance": "Significant after Benjamini-Hochberg correction:",
              "perturbation": "Perturbation support, which admits a pair but never raises its rank.",
              "same_protein": "Both sites are on one protein; read the pair as shared abundance first.",
              "protein_shadow": "The parent proteins co-vary, so the sites may simply follow them.",
              "self_coupling": "One site follows its own protein's abundance.",
              "tissue_composition": "Both sites track tissue composition, which alone can make them co-vary.",
              "heterogeneity": "The cohorts disagree beyond sampling noise.",
              "one_cohort": "Measured in one cohort only.",
              "adjusted": "The program-adjusted 95% interval excludes zero (an annotation, never a ranking).",
              "glycoform_pooling": "Different glycoforms are pooled across cohorts; see the pair's glycoforms.",
              "sequon_ambiguity": "The peptide carried more than one sequon; the site label is an attribution."}
# the next steps under a summary, by where it is read: the command line names flags, the MCP server names tools
HINTS = {"cli": {"more": "--rows 50", "columns": "--json", "pair": "pair {g} {p}", "filters": "--help",
                 "explain": "explain {g} {p}", "glycoforms": "--json", "wording": "--json", "fields": "--json",
                 "tables": "--json", "sites_of": "sites {gene}"},
         "mcp": {"more": "a higher limit; the JSON below holds {k}", "columns": "the JSON below; ptm_columns defines them",
                 "pair": 'ptm_pair with glyco_site "{g}" and phospho_site "{p}"', "filters": "the filters argument",
                 "explain": 'ptm_explain with glyco_site "{g}" and phospho_site "{p}"', "glycoforms": "the JSON below",
                 "wording": "the JSON below", "fields": "the JSON below", "tables": "the JSON below", "sites_of": 'ptm_sites with x "{gene}"'}}


def _num(x, d=2):
    return "-" if _is_na(x) else "%.*f" % (d, x)


def _q(x):
    return "-" if _is_na(x) else ("%.1e" % x if x < 0.001 else "%.3f" % x)


def _ci(lo, hi):
    return "-" if _is_na(lo) or _is_na(hi) else "[%s, %s]" % (_num(lo), _num(hi))


def _int(x):
    return "-" if _is_na(x) else "{:,}".format(int(x))


def _plural(n, one, many=None):
    return "%s %s" % (_int(n), one if n == 1 else (many or one + "s"))


def _cap(w):
    return w[:1].upper() + w[1:]


def pooled_ci(rho, sum_weight):
    if _is_na(rho) or _is_na(sum_weight) or sum_weight <= 0:
        return None, None
    z, se = math.atanh(min(max(rho, -0.999999), 0.999999)), 1 / math.sqrt(sum_weight)
    return math.tanh(z - 1.959964 * se), math.tanh(z + 1.959964 * se)


def _support(r):
    return "+".join(w for k, w in SUPPORT_WORDS if r.get(k) is True)


class _Text:
    """Lines of one summary, wrapped to `width`."""

    def __init__(self, width, style):
        self.width, self.style, self.lines = max(60, width) - 1, style, []

    def hint(self, key, **kw):
        return HINTS[self.style][key].format(**kw)

    def add(self, *lines):
        self.lines += list(lines)

    def say(self, text, indent=0):
        self.lines += textwrap.wrap(text, self.width, initial_indent=" " * indent, subsequent_indent=" " * (indent + 2)) or [""]

    def bullet(self, text):
        self.lines += textwrap.wrap(text, self.width, initial_indent="  - ", subsequent_indent="    ")

    def pack(self, units, sep="", indent=0, exdent=None):
        """Wrap between units, never inside one; `sep` ends every unit but the last."""
        units = [u for u in units if u]
        if not units:
            return
        exdent = indent + 2 if exdent is None else exdent
        units = [u + sep for u in units[:-1]] + units[-1:]
        cur = " " * indent + units[0]
        for u in units[1:]:
            if len(cur) + 1 + len(u) > self.width:
                self.lines.append(cur)
                cur = " " * exdent + u
            else:
                cur += " " + u
        self.lines.append(cur)

    def table(self, cols, right, drop=()):
        cols, right = list(cols), dict(zip([h for h, _ in cols], right))
        widths = lambda: [max([len(h)] + [len(v) for v in vals]) for h, vals in cols]
        for d in drop:
            if 2 + sum(widths()) + 2 * (len(cols) - 1) <= self.width:
                break
            cols = [c for c in cols if c[0] != d]
        w = widths()
        cell = lambda v, i: v.rjust(w[i]) if right[cols[i][0]] else v.ljust(w[i])
        self.lines.append(("  " + "  ".join(cell(h, i) for i, (h, _) in enumerate(cols))).rstrip())
        for r in range(len(cols[0][1])):
            self.lines.append(("  " + "  ".join(cell(vals[r], i) for i, (_, vals) in enumerate(cols))).rstrip())

    def definitions(self, pairs, cap=22):
        """key and value in two columns, the value wrapped under itself"""
        pairs = [(k, v) for k, v in pairs if v]
        if not pairs:
            return
        wk = min(max(len(k) for k, _ in pairs), cap)
        lead = 2 + wk + 2
        for k, v in pairs:
            body = textwrap.wrap(v, max(20, self.width - lead)) or [""]
            self.lines.append(("  " + k.ljust(wk) + "  " + body[0]).rstrip())
            self.lines += [" " * lead + b for b in body[1:]]

    def steps(self, hints, lead=()):
        """the next steps: one line each, a short label and what to do, in one column"""
        for l in lead:
            if l:
                self.lines += textwrap.wrap(l, self.width)
        hints = [(k, v) for k, v in hints if v]
        if not hints:
            return
        w = max(len(k) for k, _ in hints) + 1
        self.lines += ["  " + (k + ":").ljust(w) + "  " + v for k, v in hints]

    def text(self):
        return "\n".join(self.lines)


def _withheld_units(w, min_k=2):
    if not w:
        return []
    if all(isinstance(v, dict) for v in w.values()):   # several queries
        w = {k: sum(v.get(k, 0) for v in w.values()) for k in WITHHELD_WORDS}
    words = dict(WITHHELD_WORDS, k="single-cohort" if min_k in (None, 2) else "in fewer than %s cohorts" % min_k)
    u = ["{:,} {}".format(w[k], words[k]) for k in words if w.get(k)]
    if u:
        u[0] = "Hidden by filters: " + u[0]   # two statements: with one reason, u[0] and u[-1] are the same item
        u[-1] = u[-1] + "."
    return u


def _release_label(release):
    return "(PTMlink %s)" % release if release else "(PTMlink)"


def _desc(vals):
    v = [x for x in vals if not _is_na(x)]
    return all(a >= b for a, b in zip(v, v[1:]))


def text_pairs(rows, title, withheld=None, n=10, release=None, width=100, style="cli", n_view=None, n_matched=None, min_k=2, by=None, note=None, json_rows=None):
    """Ranked pairs. `rows` is every row of the result; `json_rows` the number the JSON beside the summary holds (the MCP
    server's limit); `n_view` the rows that passed before a query's own limit."""
    t = _Text(width, style)
    total = len(rows)
    qcol = "query_gene" if rows and "query_gene" in rows[0] else ("query_site" if rows and "query_site" in rows[0] else None)
    groups = list(dict.fromkeys(r[qcol] for r in rows)) if qcol else []
    stacked = len(groups) > 1
    ev = [r.get("evidence") for r in rows]
    order = "" if total < 2 else ", sorted by evidence" if _desc(ev) else (
        ", sorted by evidence within each %s" % (by or "query") if stacked and all(_desc([r.get("evidence") for r in rows if r[qcol] == g]) for g in groups) else "")
    count = "the top %s of %s" % (_int(total), _plural(n_view, "pair")) if n_view is not None and n_view > total else _plural(total, "pair")
    t.pack(["%s: %s%s" % (_cap(title), count, order), _release_label(release)])
    t.pack(_withheld_units(withheld, min_k), ",")
    if not rows and n_view:   # a query cut to no row by its own limit
        t.add("No rows returned: the limit was 0; a higher limit returns them.")
        return t.text()
    if not rows:
        hidden = sum(v for v in (withheld or {}).values() if isinstance(v, int)) if withheld and not all(isinstance(v, dict) for v in withheld.values()) else 0
        t.say(note if note else "No reported pair meets the conditions." if n_matched == 0 else
              "No pair passes the filters. To relax them: %s" % t.hint("filters") if hidden else "No reported pair.")
        return t.text()
    if json_rows == 0:
        t.add("No rows shown: the limit was 0; a higher limit returns them.")
        return t.text()
    if json_rows is not None:
        n = min(n, json_rows)
    # several genes or sites asked at once: the top rows of each, so no query is pushed off the screen by the first
    per = max(1, n // len(groups)) if stacked else n
    if stacked:
        s = [r for g in groups for r in [x for x in rows if x[qcol] == g][:per]][:max(n, len(groups))]
    else:
        s = rows[:n]
    # a field whose columns were not asked for (query --columns) is left out, never shown as "-"
    has = lambda *cols: any(c in r for r in s for c in cols)
    cis = [pooled_ci(r.get("rho_combined"), r.get("sum_weight")) for r in s]
    cols = [("glyco_site", [r["glyco_site"] for r in s], False), ("phospho_site", [r["phospho_site"] for r in s], False)]
    if has("rho_combined"):
        cols.append(("rho", [_num(r.get("rho_combined")) for r in s], True))
    if has("sum_weight") and has("rho_combined"):
        cols.append(("95% CI", [_ci(*c) for c in cis], True))
    if has("evidence"):
        cols.append(("evidence", [_num(r.get("evidence"), 1) for r in s], True))
    if has("q_value"):
        cols.append(("q", [_q(r.get("q_value")) for r in s], True))
    if has("k_cohorts"):
        cols.append(("cohorts", [_int(r.get("k_cohorts")) for r in s], True))
    if has("anch_dose", "anch_gene_editing", "anch_enzymatic"):
        cols.append(("perturbation", [_support(r) for r in s], False))
    if qcol:
        cols.insert(0, ("query", [r[qcol] for r in s], False))
    t.add("")
    t.table([(h, v) for h, v, _ in cols], [rt for _, _, rt in cols], drop=("95% CI", "q"))
    t.add("")
    held = total if json_rows is None else min(json_rows, total)
    more = total > len(s) if style == "cli" else (total > held or (n_view or 0) > total)
    t.steps([("more rows", t.hint("more", k=_int(held)) if more else ""), ("all columns", t.hint("columns")),
             ("one pair", t.hint("pair", g=s[0]["glyco_site"], p=s[0]["phospho_site"])), ("filters", t.hint("filters"))],
            lead=["Showing the top %d of each %s: %s of %s." % (per, by or "query", _int(len(s)), _int(total)) if stacked and total > len(s)
                  else "Showing %s of %s." % (_int(len(s)), _int(total)) if total > len(s) else ""])
    return t.text()


def text_status(recs):
    out = []
    for r in recs:
        if r.get("gene") is None:
            via = r.get("resolved_via") or ""
            out.append("%s: %s" % (r["query"], via if via.startswith("ambiguous") else "a site key, not a gene; the site command lists its partners" if via.startswith("a site key")
                                   else "not a recognized symbol, alias or UniProt accession"))
            continue
        how = re.sub(r" of .*$", "", r.get("resolved_via") or "")
        asked = r["query"].upper() != r["gene"].upper()
        out.append(r["gene"] + (" (asked as %s%s)" % (r["query"], ", " + how if how and how != "current symbol" else "") if asked else ""))
        for lab, side in (("glycosites:", "glyco"), ("phosphosites:", "phospho")):
            extra = ""
            if "n_cohorts_%s_detected" % side in r:
                nc, nr = r["n_cohorts_%s_detected" % side], r.get("n_served_as_" + side) or 0
                extra = "; detected in %s %s, {:,} %s".format(nr) % (nc, "cohort" if nc == 1 else "cohorts", "reported pair" if nr == 1 else "reported pairs")
            out.append("  %-14s%s%s" % (lab, r.get("status_" + side), extra))
    return "\n".join(out)


def text_pair(p, width=100, style="cli", release=None):
    if "pairs" in p:   # several pairs asked at once: each in full, in the order asked
        return "\n\n".join(text_pair(x, width, style, release) for x in p["pairs"])
    t = _Text(width, style)
    r = p.get("row")
    if p.get("input_note"):
        t.say(p["input_note"])
    if not r:
        t.say("%sot a reported pair: %s" % ("%s x %s: n" % (p["glyco_site"], p["phospho_site"]) if p.get("glyco_site") else "N", p.get("note", "")))
        return t.text()
    lo, hi = pooled_ci(r.get("rho_combined"), r.get("sum_weight"))
    words = {"dose": "tunicamycin dose response", "editing": "glycogene editing", "sialidase": "sialidase"}
    sup = _support(r)
    t.add("%s x %s" % (r["glyco_site"], r["phospho_site"]))
    t.definitions([("pair id", r.get("pair_id")),
                   ("pooled rho", "%s, 95%% CI %s" % (_num(r.get("rho_combined")), _ci(lo, hi))),
                   ("evidence", "%s, q = %s" % (_num(r.get("evidence"), 1), _q(r.get("q_value")))),
                   ("cohorts", "%s, %s tumors%s" % (r.get("k_cohorts"), _int(p.get("n_total")), "; I2 = " + _num(r["I2"]) if not _is_na(r.get("I2")) else "")),
                   ("perturbation", (", ".join(words[s] for s in sup.split("+")) if sup else "none") + (" (the pair itself was measured)" if r.get("anch_pair_observed") is True else "")),
                   ("parent proteins", "rho " + _num(r.get("rho_prot")) if "rho_prot" in r else ""),
                   ("own-protein coupling", _num(r.get("self_coupling")) if "self_coupling" in r else ""),
                   ("adjusted 95% CI", _ci(r.get("rho_adjusted_lo"), r.get("rho_adjusted_hi")) + " (an annotation; never used to rank)" if "rho_adjusted_lo" in r else "")])
    co = p.get("cohorts") or []
    if co:
        t.add("")
        t.table([("cohort", [c["label"] for c in co]), ("rho", [_num(c["rho"]) for c in co]),
                 ("95% CI", [_ci(c["lo"], c["hi"]) for c in co]), ("tumors", [_int(c["n"]) for c in co])], [False, True, True, True])
    t.add("")
    t.steps([("cautions", t.hint("explain", g=r["glyco_site"], p=r["phospho_site"])), ("glycoforms", t.hint("glycoforms")), ("all columns", t.hint("columns"))])
    return t.text()


_RULE_BASES = {}


def _rule_extra(rule, statement):
    """The part of a rule's statement that carries this pair's values (the statement after the rule's fixed sentence)."""
    if not _RULE_BASES:
        with open(_reference_path(), encoding="utf-8") as fh:
            _RULE_BASES.update({x["id"]: x["statement"] for x in json.load(fh)["rules"]})
    b = _RULE_BASES.get(rule)
    return statement[len(b):].strip() if b and statement.startswith(b) else ""


def _rule_short(rule, statement):
    extra = _rule_extra(rule, statement)
    if rule not in RULE_LEADS:
        return extra or statement
    return (RULE_LEADS[rule] + " " + extra).strip()


def text_explain(e, width=100, style="cli"):
    rows = e if isinstance(e, list) else e.get("rules", e.get("rows", []))
    t = _Text(width, style)
    if not rows:   # the pair is not reported: say so, as the pair view does, instead of printing nothing
        t.say("Not a reported pair: " + ABSENT_PAIR)
        return t.text()
    for key in dict.fromkeys((r["glyco_site"], r["phospho_site"]) for r in rows):
        rs = [r for r in rows if (r["glyco_site"], r["phospho_site"]) == key]
        t.add("%s x %s" % key)
        blank = [r for r in rs if r.get("applies") is None and str(r.get("statement") or "").startswith("Not assessed")]   # a tested rule whose value this pair lacks
        for title, sel in (("Applies:", [r for r in rs if r.get("applies") is True]), ("Context:", [r for r in rs if r.get("applies") is None and r not in blank])):
            if sel:
                t.add(title)
                for r in sel:
                    t.bullet(_rule_short(r["rule"], r["statement"]))
        no = [RULE_WORDS.get(r["rule"], r["rule"].replace("_", " ")) for r in rs if r.get("applies") is False]
        if no:
            no[0] = "Does not apply: " + no[0]
            no[-1] = no[-1] + "."
            t.pack(no, ",")
        na = [RULE_WORDS.get(r["rule"], r["rule"].replace("_", " ")) for r in blank]
        if na:
            na[0] = "Not assessed, no value: " + na[0]
            na[-1] = na[-1] + "."
            t.pack(na, ",")
    t.add("")
    t.steps([("full wording", t.hint("wording"))])
    return t.text()


def text_network(net, n=10, width=100, style="cli", json_rows=None):
    e = net.get("edges") or []
    total, gene = len(e), net.get("gene")
    t = _Text(width, style)
    t.pack(["Network of %s: %s from %s" % (gene, _plural(total, "partner link"), _plural(sum(x.get("n_site_pairs") or 0 for x in e), "site pair")), _release_label(net.get("release"))])
    t.add("One row per partner, with its strongest site pair; links have no direction.")
    t.pack(_withheld_units(net.get("withheld"), (net.get("filters") or {}).get("min_k", 2)), ",")
    if not e:
        if net.get("no_pair"):
            t.say(net["no_pair"])
        return t.text()
    if json_rows is not None:
        n = min(n, json_rows)
    if n == 0:
        t.add("No rows shown: the limit was 0; a higher limit returns them.")
        return t.text()
    s = e[:n]
    short = lambda site: site.rsplit("_", 1)[-1]   # the gene heads the column, so the site alone is shown
    own = [short(x["best_glyco_site"] if x["role"] == "glyco" else x["best_phospho_site"]) for x in s]
    other = [short(x["best_phospho_site"] if x["role"] == "glyco" else x["best_glyco_site"]) for x in s]
    t.add("")
    t.table([("partner", [x["partner"] for x in s]), ("site pairs", [_int(x["n_site_pairs"]) for x in s]), ("%s site" % gene, own), ("partner site", other),
             ("rho", [_num(x["rho_at_best"]) for x in s]), ("evidence", [_num(x["evidence_best"], 1) for x in s]),
             ("perturbation", ["yes" if x.get("anch_any") is True else "" for x in s])],
            [False, True, False, False, True, True, False], drop=("perturbation", "site pairs"))
    t.add("")
    held = total if json_rows is None else min(json_rows, total)
    more = total > len(s) if style == "cli" else total > held
    t.steps([("more rows", t.hint("more", k=_int(held)) if more else ""), ("all columns", t.hint("columns"))],
            lead=["Showing %s of %s." % (_int(len(s)), _int(total)) if total > len(s) else ""])
    return t.text()


def text_release(b, width=100, style="cli"):
    t = _Text(width, style)
    t.add("PTMlink release %s%s%s" % (b.get("release"), ", released " + b["released"] if b.get("released") else "",
                                      " (data built %s)" % b["built_utc"][:10] if b.get("built_utc") else ""))
    has = lambda k: isinstance(b.get(k), (int, float)) and not isinstance(b.get(k), bool)
    inf = b.get("inference") or {}
    pairs = ""
    if has("n_served_rows"):
        pairs = _int(b["n_served_rows"]) + (" on " + _plural(b["n_genes"], "gene") if has("n_genes") else "") + (
            "; %s with perturbation support" % _int(b["n_served_anchored"]) if has("n_served_anchored") else "")
    cohorts = ""
    if has("n_cohorts"):
        cohorts = _plural(b["n_cohorts"], "tumor cohort") + (", " + _plural(b["n_samples"], "tumor") if has("n_samples") else "") + (
            "; %s, never pooled" % _plural(b["n_normal_cohorts"], "normal adjacent tissue cohort") if has("n_normal_cohorts") and b["n_normal_cohorts"] > 0 else "")
    t.definitions([("reported pairs", pairs), ("cohorts", cohorts),
                   ("statistics", ("Student tests pooled by weighted Stouffer" if inf.get("inference_revision") == "student_normal_bh" else "Fisher-normal pooling") if inf else ""),
                   ("q values", "Benjamini-Hochberg over all scored site pairs" if inf else "")])
    t.add("")
    t.steps([("all fields", t.hint("fields"))])
    return t.text()


def text_cohorts(cohorts, release=None, width=100, style="cli"):
    t = _Text(width, style)
    tumor = [c for c in cohorts if (c.get("role") or "tumour") != "normal"]
    t.add("%s %s" % (_cap(_plural(len(cohorts), "cohort")), _release_label(release)))
    t.say("%d tumor cohorts are pooled%s." % (len(tumor), "; %d normal adjacent tissue cohorts are shown beside them and never pooled" % (len(cohorts) - len(tumor))
                                              if len(cohorts) > len(tumor) else ""))
    t.add("")
    kind = lambda c: re.sub(r",\s*normal adjacent tissue.*$", "", c.get("tumour_type") or "-")   # the role column already says normal tissue
    t.table([("cohort", [str(c.get("label") or c.get("cohort")) for c in cohorts]), ("tumor type", [kind(c) for c in cohorts]),
             ("samples", [_int(c.get("n_samples_glyco_phospho")) for c in cohorts]),
             ("role", ["normal tissue" if c.get("role") == "normal" else "tumor" for c in cohorts]),
             ("covariates", [c.get("covariates") or "-" for c in cohorts])], [False, False, True, False, False], drop=("covariates", "tumor type"))
    t.add("")
    t.steps([("sources, enrichment", t.hint("fields"))], lead=["Samples: those with both a glycoproteome and a phosphoproteome."])
    return t.text()


def text_sites(sites, release=None, n=30, width=100, style="cli", json_rows=None):
    t = _Text(width, style)
    if not sites:
        t.add("No reported site.")
        return t.text()
    genes = list(dict.fromkeys(s["gene"] for s in sites))
    acc = {s["gene"]: s.get("uniprot_acc") for s in sites}
    who = ", ".join("%s (%s)" % (g, acc[g]) if acc.get(g) else g for g in genes)
    t.say("%s: %s, %s reported %s" % (who, _plural(sum(s["side"] == "glyco" for s in sites), "glycosite"),
                                      _plural(sum(s["side"] == "phospho" for s in sites), "phosphosite"), _release_label(release)))
    if json_rows is not None:
        n = min(n, json_rows)
    if n == 0:
        t.add("No rows shown: the limit was 0; a higher limit returns them.")
        return t.text()
    s = sites[:n]
    match = lambda v: "matches" if v is True else "differs" if v is False else "-"
    t.add("")
    t.table([("site", [x["site"] for x in s]), ("side", [x["side"] for x in s]), ("cohorts", [_int(x.get("n_cohorts_served")) for x in s]),
             ("own-protein coupling", [_num(x.get("self_coupling_median")) for x in s]), ("pLDDT", [_num(x.get("plddt"), 0) for x in s]),
             ("model", [match(x.get("residue_matches_model")) for x in s])], [False, False, True, True, True, False], drop=("model", "own-protein coupling"))
    t.add("")
    total = len(sites)
    held = total if json_rows is None else min(json_rows, total)
    t.steps([("more rows", t.hint("more", k=_int(held)) if (total > len(s) if style == "cli" else total > held) else ""), ("coordinates, all columns", t.hint("columns"))],
            lead=["pLDDT: the AlphaFold model's confidence at the site.", "Showing %s of %s." % (_int(len(s)), _int(total)) if total > len(s) else ""])
    return t.text()


def text_glycans(glycans, release=None, n=10, width=100, style="cli", json_rows=None):
    t = _Text(width, style)
    total = len(glycans)
    t.add("%s glycan compositions, the most often measured first %s" % (_int(total), _release_label(release)))
    if json_rows is not None:
        n = min(n, json_rows)
    if n == 0:
        t.add("No rows shown: the limit was 0; a higher limit returns them.")
        return t.text()
    s = glycans[:n]
    t.add("")
    t.table([("composition", [str(g.get("composition")) for g in s]), ("class", [str(g.get("glycan_class") or "-") for g in s]),
             ("notation", [str(g.get("community_notation") or "-") for g in s]), ("mass (Da)", [_num(g.get("monoisotopic_mass"), 1) for g in s]),
             ("cohorts", [_int(g.get("n_cohorts_observed")) for g in s])], [False, False, False, True, True], drop=("mass (Da)", "notation"))
    t.add("")
    held = total if json_rows is None else min(json_rows, total)
    t.steps([("more rows", t.hint("more", k=_int(held)) if (total > len(s) if style == "cli" else total > held) else ""), ("all columns", t.hint("columns"))],
            lead=["Composition: N HexNAc, H hexose, F fucose, S NeuAc, G NeuGc.", "Showing %s of %s." % (_int(len(s)), _int(total)) if total > len(s) else ""])
    return t.text()


def text_cite(c, width=100, style="cli"):
    t = _Text(width, style)
    t.add("To cite PTMlink:", "")
    t.lines += textwrap.wrap(c.get("text") or "", t.width, initial_indent="  ", subsequent_indent="  ")
    t.add("")
    if c.get("release"):
        t.add("Name the data release you used: %s." % c["release"])
    return t.text()


def text_list(l, n=5, width=100, style="cli"):
    """A list of proteins or sites: what it holds, how its members pair with each other, and what they reach."""
    t = _Text(width, style)
    m, cn, site = l["members"], l["connectivity"], l["level"] == "site"
    release = (l.get("provenance") or {}).get("release")
    n_units = len(l["sites"]) if site else len(dict.fromkeys(x["gene"] for x in m if x["gene"]))
    t.add("List of %s against the reported pairs %s" % (_plural(n_units, "site" if site else "protein"), _release_label(release)))
    t.say("Pairs in view: %s within the list, %s with a partner outside it%s." % (
        _int(len(l["internal"])), _int(len(l["neighbourhood"])), ", %s on one protein" % _int(len(l["same_protein"])) if l["same_protein"] else ""))
    t.pack(_withheld_units(l.get("withheld"), ((l.get("provenance") or {}).get("filters") or {}).get("min_k", 2)), ",")
    if l.get("mixed"):
        k = sum(1 for x in m if x["kind"] == "gene" and x["gene"])
        t.say("Read at site level: %s expanded to reported sites." % _plural(k, "gene entry", "gene entries"))
    lost = [x["query"] for x in m if not x["gene"]]
    if lost:
        t.pack(["Not recognized: " + lost[0]] + lost[1:] if len(lost) > 1 else ["Not recognized: " + lost[0] + "."], ",")
        if len(lost) > 1:
            t.lines[-1] += "."
    t.add("", "Members")
    s = m[:10]
    # counted by site, an entry shows its own reported pairs; counted by protein, the protein's pairs on each side
    if site:
        t.table([("entry", [x["query"] for x in s]), ("protein", [x["gene"] or "-" for x in s]), ("reported pairs", [_int(x.get("n_reported")) for x in s]),
                 ("within list", [_int(x.get("n_internal")) for x in s]), ("outside list", [_int(x.get("n_neighbourhood")) for x in s])], [False, False, True, True, True])
    else:
        t.table([("entry", [x["query"] for x in s]), ("protein", [x["gene"] or "-" for x in s]), ("glyco pairs", [_int(x.get("n_served_as_glyco")) for x in s]),
                 ("phospho pairs", [_int(x.get("n_served_as_phospho")) for x in s]), ("within list", [_int(x.get("n_internal")) for x in s]),
                 ("outside list", [_int(x.get("n_neighbourhood")) for x in s])], [False, False, True, True, True, True], drop=("phospho pairs", "glyco pairs"))
    if len(m) > len(s):
        t.add("  %s more entries" % _int(len(m) - len(s)))
    t.add("", "Connectivity among the members")
    exp_, ratio, pp = cn.get("expected_pairs"), cn.get("ratio"), cn.get("p_perm")
    if exp_ == 0 and cn.get("internal_pairs") == 0:
        one = len(dict.fromkeys(x["gene"] for x in m if x["gene"])) < 2
        t.say("No pair within the list is possible: " + ("it holds one protein, and pairs within a protein are not counted." if one else
              "its members with significant pairs are all glycosites, or all phosphosites."), indent=2)
    else:
        t.say("%s within the list; %s expected from the members' own pair counts%s.%s" % (
            _plural(cn.get("internal_pairs") or 0, "significant pair"), "-" if _is_na(exp_) else "{:,.1f}".format(exp_),
            "" if _is_na(ratio) else " (%.2f times as many)" % ratio,
            "" if _is_na(pp) else " Permutation P = %s over %s random lists." % (_q(pp), _int(cn.get("n_perm")))), indent=2)

    def top(tab, title, key, head):
        t.add("", title)
        if not tab:
            t.add("  none")
            return
        d = [o for o in tab if o["members_reached"] >= 2]
        # a partner reached by one member is a single pair, already in the row tables: the summary shows the shared ones,
        # and says so when there is none (a fold computed on one member is not a finding)
        if not d:
            t.say("None is reached by two or more members; %s reached by one member each (single pairs, in %s)." % (_int(len(tab)), t.hint("tables")), indent=2)
            return
        k = d[:n]
        t.table([(head, [o[key] for o in k]), ("side", [o["side" if key == "site" else "partner_side"] for o in k]),
                 ("members", ["%s of %s" % (_int(o["members_reached"]), _int(o["members_testable"])) for o in k]), ("fold", [_num(o["fold"], 1) for o in k]),
                 ("q", [_q(o["q"]) for o in k]), ("best evidence", [_num(o.get("evidence_best"), 1) for o in k])], [False, False, True, True, True, True], drop=("best evidence",))
        sig = sum(1 for o in tab if not _is_na(o.get("q")) and o["q"] < 0.05)
        t.say("%s%s%s; %s at q < 0.05." % ("Top %d by P of " % len(k) if len(d) > len(k) else "", _int(len(d)),
                                           " reached by two or more members (%s reached in all)" % _int(len(tab)) if len(d) < len(tab) else " reached", _int(sig)), indent=2)
    top(l["partner_sites"], "Partner sites reached by more members than expected", "site", "site")
    top(l["convergent"], "Partner proteins reached by more members than expected", "partner", "protein")
    t.add("")
    uni = (l.get("universe") or {}).get("pairs")
    t.say("The three tests count %ssignificant pairs (q < 0.05) between two different proteins." % (_int(uni) + " " if uni is not None else ""))
    share = next((c.get("share") for c in l.get("composition") or [] if c.get("table") == "neighbourhood"), None)
    if not _is_na(share) and share >= 0.05:
        t.say("%.0f%% of the pairs with a partner outside the list carry a tissue-composition label on both sites." % (100 * share))
    t.add("")
    t.steps([("every table", t.hint("tables"))])
    return t.text()


def to_text(cmd, out, release=None, n=10, width=100, style="cli", json_rows=None):
    """The readable form of a result; None when the command has no text form (its JSON is printed instead). `out` is
    the full result; `json_rows` is the number of rows the JSON beside the summary holds, where it is capped."""
    if cmd == "status":
        return text_status(out)
    if cmd in ("gene", "site", "query"):
        rows = out.get("rows", [])
        min_k = (out.get("filters") or {}).get("min_k", 2)
        if cmd == "gene":
            genes = out.get("genes") or [out.get("gene", "")]
            title, by = "%s %s" % (", ".join(genes), "glycosites" if out.get("side", "glyco") == "glyco" else "phosphosites"), "gene"
        elif cmd == "site":
            title, by = "partners of " + ", ".join(out.get("sites") or [out.get("site", "")]), "site"
        else:
            title, by = "whole-release query", None
        return text_pairs(rows, title, out.get("withheld"), n=n, release=release, width=width, style=style, n_view=out.get("n_view"),
                          n_matched=out.get("n_matched"), min_k=min_k, by=by, note=out.get("note"), json_rows=json_rows)
    if cmd == "pair":
        return text_pair(out, width, style, release)
    if cmd == "explain":
        return text_explain(out, width, style)
    if cmd == "network":
        return text_network(out, n, width, style, json_rows)
    if cmd == "release":
        return text_release(out, width, style)
    if cmd == "cohorts":
        return text_cohorts(out, release, width, style)
    if cmd == "sites":
        return text_sites(out, release, max(n, 30), width, style, json_rows)
    if cmd == "glycans":
        return text_glycans(out, release, n, width, style, json_rows)
    if cmd == "cite":
        return text_cite(out, width, style)
    if cmd == "list":
        return text_list(out, 5, width, style)
    return None


def _filters_from_args(a):
    f = dict(DEFAULT_FILTERS)
    if getattr(a, "include_cis", False):
        f["include_cis"] = True
    if getattr(a, "min_k", None) is not None:
        f["min_k"] = a.min_k
    if getattr(a, "hide_composition", False):
        f["hide_composition"] = True
    if getattr(a, "anchored_only", False):
        f["anchored_only"] = True
    if getattr(a, "significant_only", False):
        f["significant_only"] = True
    if getattr(a, "no_significant_only", False):   # the q filter is on by default; this unticks it, as the website's box does
        f["significant_only"] = False
    if getattr(a, "max_i2", None) is not None:
        f["max_i2"] = a.max_i2
    if getattr(a, "adjusted_excludes_zero", False):
        f["adjusted_excludes_zero"] = True
    if getattr(a, "shadow_clean_only", False):
        f["shadow_clean_only"] = True
    if getattr(a, "hide_self_coupled", False):
        f["hide_self_coupled"] = True
    if getattr(a, "min_loco_pctile", None) is not None:
        f["min_loco_pctile"] = a.min_loco_pctile
    return f


def _rows_arg(v):
    try:
        n = int(v)
    except ValueError:
        n = 0
    if n < 1:
        raise argparse.ArgumentTypeError("a whole number of rows, 1 or more; got %s" % v)
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(description="PTMlink reader: prints a readable summary, or the full result with --json")
    ap.add_argument("--json", action="store_true", help="print the full result as JSON instead of the readable summary")
    ap.add_argument("--data-dir", help="a local copy of the website's data directory (the one holding build.json)")
    ap.add_argument("--data-url", help="the data address to download from (default %s)" % DEFAULT_DATA_URL)
    ap.add_argument("--release-dir", help="the Parquet release, needed by query and list")
    ap.add_argument("--no-verify", action="store_true", help="accept the release directory without checking its checksum")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def add_filters(p):
        p.add_argument("--include-cis", action="store_true"); p.add_argument("--min-k", type=int); p.add_argument("--hide-composition", action="store_true")
        p.add_argument("--anchored-only", action="store_true", help="keep only pairs with perturbation evidence")
        p.add_argument("--significant-only", action="store_true", help="keep only pairs at q < 0.05 on the total score (the default view already does)")
        p.add_argument("--no-significant-only", action="store_true", help="turn off the q < 0.05 filter: return every reported pair the other filters keep")
        p.add_argument("--max-i2", type=float)
        p.add_argument("--adjusted-excludes-zero", action="store_true"); p.add_argument("--shadow-clean-only", action="store_true"); p.add_argument("--hide-self-coupled", action="store_true")
        p.add_argument("--min-loco-pctile", type=float)

    def add_rows(p):
        p.add_argument("--rows", type=_rows_arg, default=10, help="rows to show in the summary, 1 or more (the --json result has them all)")
    p = sub.add_parser("status"); p.add_argument("gene", nargs="+")
    p = sub.add_parser("gene"); p.add_argument("gene", nargs="+"); p.add_argument("--side", default="glyco"); add_filters(p); add_rows(p)
    p = sub.add_parser("site"); p.add_argument("site", nargs="+"); add_filters(p); add_rows(p)
    p = sub.add_parser("pair", help="one pair: two site keys, or one pair_id"); p.add_argument("glyco_site", help="glycosite key (IL6ST_N131), or a pair_id (\"PTMG:P40189:N131|PTMP:Q9H3N1:S270\")"); p.add_argument("phospho_site", nargs="?", help="phosphosite key (TMX1_S270)")
    p = sub.add_parser("network"); p.add_argument("gene"); p.add_argument("--top-n", type=int, default=25); p.add_argument("--pin", nargs="*", default=[]); add_filters(p); add_rows(p)
    p = sub.add_parser("explain", help="the reading rules on one pair: two site keys, or one pair_id"); p.add_argument("glyco_site"); p.add_argument("phospho_site", nargs="?")
    p = sub.add_parser("sites"); p.add_argument("x", nargs="+"); p.add_argument("--no-structure", action="store_true")
    p = sub.add_parser("query"); p.add_argument("--where", nargs=3, action="append", metavar=("COLUMN", "OP", "VALUE"), default=[]); p.add_argument("--genes", nargs="*"); p.add_argument("--sites", nargs="*")
    p.add_argument("--columns", nargs="*"); p.add_argument("--limit", type=int); p.add_argument("--allow-full", action="store_true"); add_filters(p); add_rows(p)
    p = sub.add_parser("list"); p.add_argument("x", nargs="+"); p.add_argument("--level", default="auto", choices=["auto", "protein", "site"]); p.add_argument("--n-perm", type=int, default=200); p.add_argument("--seed", type=int, default=1); add_filters(p)
    p = sub.add_parser("tools"); p.add_argument("--format", default="mcp", choices=["mcp", "anthropic", "openai"])
    sub.add_parser("columns"); sub.add_parser("reference"); sub.add_parser("cite"); sub.add_parser("release"); sub.add_parser("cohorts")
    p = sub.add_parser("glycans"); add_rows(p)
    for sp in sub.choices.values():   # --json is read before or after the command
        sp.add_argument("--json", action="store_true", default=argparse.SUPPRESS, help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    try:
        return _run(a)
    except KeyboardInterrupt:
        return 130
    except Exception as e:   # noqa: BLE001 -- one line for the user; PTMLINK_DEBUG=1 shows the traceback
        if os.environ.get("PTMLINK_DEBUG"):
            raise
        msg = str(e) if isinstance(e, (PTMlinkError, ValueError, LookupError, FileNotFoundError)) else "%s: %s" % (type(e).__name__, e)
        sys.stderr.write("ptmlink: %s\n" % msg)
        return 1


def _run(a):
    db = PTMlink(data_dir=a.data_dir, data_url=a.data_url, release_dir=a.release_dir, verify_release=not a.no_verify)
    if a.cmd == "status":
        out = db.status(a.gene)
    elif a.cmd == "gene":
        out = db.gene(a.gene if len(a.gene) > 1 else a.gene[0], a.side, _filters_from_args(a))
    elif a.cmd == "site":
        out = db.site(a.site if len(a.site) > 1 else a.site[0], _filters_from_args(a))
    elif a.cmd == "pair":
        out = db.pair(a.glyco_site, a.phospho_site)
    elif a.cmd == "network":
        out = db.network(a.gene, _filters_from_args(a), a.top_n, tuple(a.pin))
    elif a.cmd == "explain":
        out = db.explain(a.glyco_site, a.phospho_site)
    elif a.cmd == "sites":
        out = db.sites(a.x, structure=not a.no_structure)
    elif a.cmd == "query":
        where = [(c, op, (v.split(",") if op in ("in", "not_in") else v)) for c, op, v in a.where]
        out = db.query(where, a.genes, a.sites, a.columns, _filters_from_args(a), a.limit, a.allow_full, _text_values=True)
    elif a.cmd == "list":
        out = db.list(a.x, a.level, _filters_from_args(a), a.n_perm, a.seed)
    elif a.cmd == "tools":
        out = db.tools(a.format)
    elif a.cmd == "columns":
        out = db.columns()
    elif a.cmd == "reference":
        out = db.reference()
    elif a.cmd == "cite":
        out = db.cite()
    elif a.cmd == "cohorts":
        out = db.cohorts()
    elif a.cmd == "glycans":
        out = db.glycans()
    else:
        out = db.release()
    text = None if a.json else to_text(a.cmd, out, release=None if a.cmd in ("status", "release", "tools", "columns", "reference") else db.release().get("release"),
                                       n=getattr(a, "rows", 10), width=shutil.get_terminal_size((100, 24)).columns)
    if text is not None:
        sys.stdout.write(text + "\n")
        return 0
    json.dump(out, sys.stdout, indent=1, default=str)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
