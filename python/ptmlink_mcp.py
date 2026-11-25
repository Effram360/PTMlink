#!/usr/bin/env python3
"""ptmlink_mcp -- PTMlink as a Model Context Protocol server over stdio, with no dependency beyond the reader.

The tools are the ones ptmlink/inst/reference.json marks as tools, described by PTMlink.tools() from that
file (the R package's ptm_tools() generates the same descriptions). Every call runs the reader's method of
the same name without the ptm_ prefix and returns a readable summary followed by its JSON. Functions that
return rows are capped at `limit` rows (default from the reference's tool conventions) with the full count
reported beside them, and each row keeps its key columns, so a tool result stays a readable size.

    python3 python/ptmlink_mcp.py                       # reads the public website's data
    python3 python/ptmlink_mcp.py --data-dir path/to/data --release-dir path/to/release   # optional local copies

An MCP client configuration names that command; the server's `instructions` carry the release contract
(one ranking, no direction, perturbation never in the score, the default view) so the agent reads it
before its first call. The transport is newline-delimited JSON-RPC 2.0 on stdin and stdout, as the MCP
stdio transport specifies; nothing else is written to stdout.
"""
import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ptmlink import PTMlink, PTMlinkError, pooled_ci, to_text  # noqa: E402

PROTOCOL = "2025-06-18"
# the versions this stdio server speaks; a client asking another one is answered with PROTOCOL, as the spec says
SUPPORTED_PROTOCOLS = ("2024-11-05", "2025-03-26", PROTOCOL)   # protocolVersion values

# The tables a result may hold, each capped at `limit` with its full count reported beside it.
CAPPED = ("rows", "edges", "internal", "same_protein", "neighbourhood", "partner_sites", "convergent")
LIST_PAIR_TABLES = ("internal", "same_protein", "neighbourhood")   # ptm_list's tables of site pairs

# The columns a row keeps in a tool result: identity, the score, its interval, the annotations and perturbation support.
# A row's packed per-cohort strings and its working columns are left out to keep a result readable; ptm_pair returns
# one pair with every column and its cohorts, and ptm_query takes `columns` for any other set.
ROW_COLUMNS = ("query_gene", "query_site", "pair_id", "glyco_site", "phospho_site", "glyco_gene", "phospho_gene", "rho_combined",
               "rho_lo", "rho_hi", "evidence", "q_value", "k_cohorts", "n_total", "I2", "anch_any", "anch_dose",
               "anch_gene_editing", "anch_enzymatic", "anch_pair_observed", "rho_prot", "self_coupling", "glyco_programs",
               "phospho_programs", "rho_adjusted_lo", "rho_adjusted_hi", "loco_pctile_k", "evidence_pctile_k", "is_cis")
# The columns a network edge keeps: the partner, how many site pairs stand behind the link, and its strongest pair.
EDGE_COLUMNS = ("gene", "partner", "role", "rank", "drawn", "n_site_pairs", "sign_mixed", "evidence_best", "rho_at_best", "k_at_best", "n_at_best",
                "I2_at_best", "q_at_best", "rho_prot", "self_coupling_at_best", "anch_any", "best_glyco_site", "best_phospho_site")
# the tools that have a readable summary, and the reader's name for each
TEXT_FORM = {"ptm_status": "status", "ptm_gene": "gene", "ptm_site": "site", "ptm_pair": "pair", "ptm_explain": "explain",
             "ptm_network": "network", "ptm_query": "query", "ptm_release": "release", "ptm_cohorts": "cohorts",
             "ptm_sites": "sites", "ptm_glycans": "glycans", "ptm_cite": "cite", "ptm_list": "list"}
KEY_COLUMNS_NOTE = "key columns only; ptm_pair gives one pair in full, ptm_columns defines each column"

# The reader words its messages for the command line; an assistant is pointed at the tool of the same job instead.
MCP_WORDS = (
    (re.compile(r"the site command \(site\(\) in Python, ptm_site in MCP\)"), "ptm_site"),
    (re.compile(r"the site command lists its partners"), "ptm_site lists its partners"),
    (re.compile(r"the sites command lists"), "ptm_sites lists"),
    (re.compile(r"\bstatus \S+ says"), "ptm_status says"),
    (re.compile(r"\bstatus says"), "ptm_status says"),
    (re.compile(r"; see status\(\)"), "; ptm_status says how names resolve"),
    (re.compile(r"\bcolumns\(\) lists them|the columns command lists them"), "ptm_columns lists them"),
    (re.compile(r"^query and list read the full table of reported pairs \(site_served\.parquet\), which is released with the published article; once you have it, .*$"),
     "ptm_query and ptm_list read the full table of reported pairs, which is released with the published article. This server was started "
     "without it (its --release-dir option), so these two tools are not available here; the other tools are."),
    (re.compile(r"^query and list need pyarrow: pip install pyarrow$"),
     "ptm_query and ptm_list need the pyarrow package on the computer that runs this server (pip install pyarrow); the other tools do not."),
    (re.compile(r"^pair takes "), "ptm_pair takes "), (re.compile(r"^explain takes "), "ptm_explain takes "), (re.compile(r"^network takes "), "ptm_network takes "),
    (re.compile(r"\(IL6ST_N131 TMX1_S270\)"), '(glyco_site "IL6ST_N131", phospho_site "TMX1_S270")'),
    (re.compile(r"one pair_id \("), "one pair_id as glyco_site ("),
    (re.compile(r"--data-dir \(data_dir= in Python\)"), "the server's --data-dir option"),
    (re.compile(r"--data-url \(data_url= in Python\)"), "the server's --data-url option"),
    (re.compile(r", or set allow_full to read every reported pair"), ", or set allow_full to true to read every reported pair"),
)


def for_mcp(text):
    for pat, new in MCP_WORDS:
        text = pat.sub(new, text)
    return text


class InvalidParams(ValueError):
    """A request whose params are not what the method takes: JSON-RPC error -32602."""


class InvalidRequest(ValueError):
    """A message that is not a JSON-RPC 2.0 request: error -32600."""


# An argument that keeps its meaning when it is null; any other null argument is read as "not given".
NULL_MEANS_SOMETHING = ("top_n",)
MAX_PAIRS = 25   # pairs in one ptm_pair or ptm_explain call: each comes back in full


def _type_ok(value, schema):
    """Does a value fit a JSON Schema fragment of the tool descriptions (type, anyOf, items)?"""
    if "anyOf" in schema:
        return any(_type_ok(value, s) for s in schema["anyOf"])
    t = schema.get("type")
    if t == "string":
        return isinstance(value, str)
    if t == "integer":
        return isinstance(value, int) and not isinstance(value, bool)
    if t == "number":
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if t == "boolean":
        return isinstance(value, bool)
    if t == "null":
        return value is None
    if t == "array":
        return (isinstance(value, list) and len(value) >= schema.get("minItems", 0) and len(value) <= schema.get("maxItems", len(value))
                and all(_type_ok(v, schema.get("items", {})) for v in value))
    if t == "object":
        return isinstance(value, dict)
    return True


def _type_words(schema):
    if "anyOf" in schema:
        kinds = {s.get("type") for s in schema["anyOf"]}
        if kinds == {"string", "array"}:
            return "text, or a list of text"
        if "null" in kinds and len(kinds) == 2:
            return _type_words(next(s for s in schema["anyOf"] if s.get("type") != "null")) + ", or null"
        return "one of several kinds of value"
    return {"string": "text", "integer": "a whole number", "number": "a number", "boolean": "true or false",
            "array": "a list", "object": "a set of named settings"}.get(schema.get("type"), "another kind of value")


def release_line(db):
    """The release this server reads, said from its data: the label and the cohorts, counted from cohorts.json (a
    normal adjacent tissue (NAT) cohort, role "normal", is view-only and never counted as pooled). Empty when the data source cannot
    be read at start-up; the reference's instructions stand alone then and ptm_release / ptm_cohorts still say it."""
    try:
        rel = db.release()
        cohorts = db.cohorts() or []
    except Exception:
        return ""
    # display names (GC, ccRCC, Rare kidney cancers); ptm_cohorts gives each one's analysis key, which arguments use
    tumor = [c.get("label") or c.get("cohort") for c in cohorts if (c.get("role") or "tumour") != "normal"]
    normal = [c.get("label") or c.get("cohort") for c in cohorts if c.get("role") == "normal"]
    line = " This data source: release %s; %d tumor cohorts pooled (%s)" % (rel.get("release"), len(tumor), ", ".join(map(str, tumor)))
    if normal:
        line += "; %d normal adjacent tissue (NAT) cohorts shown beside them, view-only (%s)" % (len(normal), ", ".join(map(str, normal)))
    return line + "."


class Server:
    def __init__(self, db):
        self.db = db
        spec = db.tools("mcp")
        self.instructions = spec["instructions"] + release_line(db)
        if db.messages:   # said while reading the release at start-up (an offline cache, for example)
            self.instructions += " Note from the reader: " + " ".join(db.messages)
            db.messages.clear()
        self.tools = spec["tools"]
        self.rows = {f["name"]: bool(f.get("rows")) for f in db.reference()["functions"]}
        conv = db.reference()["tool_conventions"]["row_limit"]
        self.default_limit, self.max_limit = conv["default"], conv.get("max", 500)
        for t in self.tools:   # ptm_query's own limit follows this server's cap, and its schema says so
            lim = t["inputSchema"]["properties"].get("limit")
            if lim is not None and "default" not in lim:
                lim["default"] = self.default_limit
                lim["description"] = "%s (the MCP server's cap: %d unless given, at most %d)" % (lim.get("description", "").replace(" (default None)", ""), self.default_limit, self.max_limit)
        try:
            self.release = db.release().get("release")
        except Exception:
            self.release = None

    def _check(self, name, args):
        """The arguments against the tool's schema: an unknown argument or filter, a missing required one, or a value of
        the wrong kind is a clear error rather than a silent default or a Python exception."""
        schema = next(t["inputSchema"] for t in self.tools if t["name"] == name)
        props = schema.get("properties", {})
        unknown = [k for k in args if k not in props]
        if unknown:
            raise ValueError("%s has no argument %s; it takes %s" % (name, ", ".join('"%s"' % k for k in unknown), ", ".join(props) or "no arguments"))
        missing = [k for k in schema.get("required", []) if args.get(k) in (None, "", [])]
        if missing:
            raise ValueError("%s needs %s" % (name, ", ".join(missing)))
        if "limit" in args and (not isinstance(args["limit"], int) or isinstance(args["limit"], bool) or not 0 <= args["limit"] <= self.max_limit):
            raise ValueError("limit is a whole number from 0 to %d; a narrower question (filters, ptm_site, ptm_query with conditions) reaches the rows a longer answer would hold" % self.max_limit)
        if name in ("ptm_pair", "ptm_explain"):
            n = max(len(args[k]) for k in ("glyco_site", "phospho_site") if isinstance(args.get(k), list)) if any(isinstance(args.get(k), list) for k in ("glyco_site", "phospho_site")) else 1
            if n > MAX_PAIRS:
                raise ValueError("%s reads at most %d pairs in one call; got %d" % (name, MAX_PAIRS, n))
        for k, v in args.items():
            if v is None or k == "filters":
                continue
            if not _type_ok(v, props[k]):
                if k == "where":
                    raise ValueError('where is a list of conditions, each [column, operator, value], such as [["q_value", "<", 1e-10], ["k_cohorts", ">=", 5]]; '
                                     'the operators is_na and not_na take no value')
                raise ValueError("%s: %s is %s; got %s" % (name, k, _type_words(props[k]), json.dumps(v)[:60]))
            if "enum" in props[k] and v not in props[k]["enum"]:
                raise ValueError("%s: %s is one of %s; got %s" % (name, k, ", ".join(map(str, props[k]["enum"])), json.dumps(v)[:60]))
        f = args.get("filters")
        if f is not None:
            if not isinstance(f, dict):
                raise ValueError('filters is a set of named settings, such as {"min_k": 1}')
            known = props["filters"].get("properties", {})
            bad = [k for k in f if k not in known]
            if bad:
                raise ValueError("unknown filter %s; the filters are %s" % (", ".join('"%s"' % k for k in bad), ", ".join(known)))
            for k, v in f.items():
                if v is not None and not _type_ok(v, known[k]):
                    raise ValueError("filter %s is %s; got %s" % (k, _type_words(known[k]), json.dumps(v)[:60]))

    def call(self, name, args):
        """Run one tool: the full result, and the row cap that applies to it (None for a tool without rows)."""
        if not isinstance(name, str) or name not in {t["name"] for t in self.tools}:
            raise InvalidParams("no tool named %s; tools/list names them" % (name if isinstance(name, str) else json.dumps(name)[:60]))
        if args is not None and not isinstance(args, dict):
            raise InvalidParams("arguments is a set of named values, such as {\"gene\": \"EGFR\"}")
        # null is "not given", except where null itself is a setting
        args = {k: v for k, v in (args or {}).items() if v is not None or k in NULL_MEANS_SOMETHING}
        if name == "ptm_query" and isinstance(args.get("where"), list) and args["where"] and not isinstance(args["where"][0], list):
            args["where"] = [args["where"]]   # one bare condition, as the reader accepts it
        self._check(name, args)
        limit = None
        if self.rows.get(name):
            limit = args.pop("limit", None)
            if limit is None:
                limit = self.default_limit
            if name == "ptm_query":   # the reader cuts the ranking itself and reports how many rows passed
                args["limit"] = limit
        method = getattr(self.db, name[len("ptm_"):])
        return method(**args), limit, args

    @staticmethod
    def _pairs(rows, also=()):
        """site-pair rows cut to ROW_COLUMNS, with the pooled 95 % interval (rho_lo, rho_hi) beside rho"""
        out = []
        for r in rows:
            r = dict(r)
            r["rho_lo"], r["rho_hi"] = pooled_ci(r.get("rho_combined"), r.get("sum_weight"))
            out.append({k: r[k] for k in ROW_COLUMNS + tuple(also) if k in r})
        return out

    @staticmethod
    def _top(rows, limit):
        """The first `limit` rows; when several genes or sites were asked at once, the top rows of each in turn, so the
        answer holds every query and not the first one alone (the summary shows the same rows)."""
        qcol = next((c for c in ("query_gene", "query_site") if rows and c in rows[0]), None)
        groups = list(dict.fromkeys(r[qcol] for r in rows)) if qcol else []
        if len(groups) < 2:
            return rows[:limit]
        per = max(1, limit // len(groups))
        return [r for g in groups for r in [x for x in rows if x[qcol] == g][:per]][:limit]

    def shape(self, name, result, limit, args):
        """The JSON part of a tool result: every table capped at `limit` with its full count beside it, site pairs
        and network links cut to their key columns, and the bulky side tables left to the tools that own them."""
        if isinstance(result, list):   # a plain table: sites, glycans (capped: always rows and their full count), cohorts, columns, an explanation
            if limit is not None:
                return {"rows": result[:limit], "n_rows_total": len(result)}
            return result
        if not isinstance(result, dict):
            return result
        out = dict(result)
        if limit is not None:
            for key in CAPPED:
                if isinstance(out.get(key), list) and len(out[key]) > limit:
                    out["n_%s_total" % key] = len(out[key])
                    out[key] = self._top(out[key], limit) if key == "rows" else out[key][:limit]
        if name in ("ptm_gene", "ptm_site") or (name == "ptm_query" and not args.get("columns")):
            out["rows"] = self._pairs(out.get("rows") or [])
            out["row_columns"] = KEY_COLUMNS_NOTE
        if name == "ptm_gene" and "sites" in out:   # the gene's site table belongs to ptm_sites, not to every ranked result
            del out["sites"]
        if name == "ptm_network":
            out["edges"] = [{k: e[k] for k in EDGE_COLUMNS if k in e} for e in out.get("edges") or []]
            out["edge_columns"] = "key columns only; ptm_pair gives the strongest pair of a link in full"
        if name == "ptm_list":
            for t in LIST_PAIR_TABLES:
                if isinstance(out.get(t), list):
                    out[t] = self._pairs(out[t], ("site_named_glyco", "site_named_phospho"))
            out["row_columns"] = KEY_COLUMNS_NOTE
        return out

    def summary(self, name, result, limit):
        """The readable summary sent ahead of the JSON; None for a tool that has none."""
        form = TEXT_FORM.get(name)
        if form:
            text = to_text(form, result, release=self.release, style="mcp", json_rows=limit)
            return for_mcp(text) if text else text
        if name == "ptm_columns":
            return "%d columns of a reported pair; the JSON below gives each one's type and meaning." % len(result)
        if name == "ptm_reference" and isinstance(result, dict):
            return "The PTMlink reference: %d reading rules, %d functions, %d columns and %d filters, in the JSON below." % (
                len(result.get("rules", [])), len(result.get("functions", [])), len(result.get("columns", [])), len(result.get("filters", [])))
        return None

    def handle(self, msg):
        method, params = msg.get("method"), msg.get("params")
        if params is None:
            params = {}
        if not isinstance(params, dict):
            raise InvalidParams("params is a JSON object of named values")
        if method == "initialize":
            asked = params.get("protocolVersion")
            return {"protocolVersion": asked if asked in SUPPORTED_PROTOCOLS else PROTOCOL, "capabilities": {"tools": {}},
                    "serverInfo": {"name": "ptmlink", "version": self.db.reference().get("version", "")}, "instructions": self.instructions}
        if method == "tools/list":
            return {"tools": self.tools}
        if method == "tools/call":
            name = params.get("name")
            self.db.messages.clear()
            try:
                result, limit, args = self.call(name, params.get("arguments"))
                # a readable summary first, then the result as JSON
                content = []
                try:
                    text = self.summary(name, result, limit)
                except Exception:   # a summary that cannot be drawn never hides the result
                    text = None
                if text:
                    content.append({"type": "text", "text": text})
                # what the reader said along the way (a name it skipped, a structure model it could not read, a swapped pair)
                said = list(self.db.messages)
                if said:
                    content.insert(0, {"type": "text", "text": for_mcp("\n".join(said))})
                content.append({"type": "text", "text": json.dumps(self.shape(name, result, limit, args), default=str)})
                return {"content": content, "isError": False}
            except InvalidParams:   # an unknown tool or arguments that are not an object: a protocol error
                raise
            except (PTMlinkError, ValueError, LookupError, FileNotFoundError) as e:   # said to the agent in the reader's own words
                return {"content": [{"type": "text", "text": for_mcp(str(e))}], "isError": True}
            except MemoryError:
                return {"content": [{"type": "text", "text": "the answer did not fit in this computer's memory; ask a narrower question"}], "isError": True}
            except Exception as e:  # anything else goes back as text, never as a crash, and without the program's internals
                sys.stderr.write("ptmlink_mcp: %s in %s: %s\n" % (type(e).__name__, name, e))
                return {"content": [{"type": "text", "text": "%s could not answer this call (an internal error, written to the server's log); "
                                                             "check the arguments against the tool's description" % name}], "isError": True}
        if method == "ping":
            return {}
        raise LookupError(method)

    def answer(self, msg):
        """The reply to one message, or None when it takes none (a notification, or a response sent to this server)."""
        def err(code, text, rid=None):
            return {"jsonrpc": "2.0", "id": rid, "error": {"code": code, "message": text}}
        if not isinstance(msg, dict):
            return err(-32600, "a request is one JSON object")
        rid = msg.get("id")
        id_ok = "id" not in msg or (isinstance(rid, (str, int)) and not isinstance(rid, bool))
        if "method" not in msg:
            if "result" in msg or "error" in msg:   # a response: this server sends no request, so it awaits none
                return None
            return err(-32600, "a request names a method", rid if id_ok else None)
        if msg.get("jsonrpc") != "2.0" or not isinstance(msg["method"], str) or not id_ok:
            return err(-32600, 'a request carries "jsonrpc": "2.0", a method name and, unless it is a notification, a text or whole-number id', rid if id_ok else None)
        if "id" not in msg:   # a notification: nothing is answered
            return None
        try:
            return {"jsonrpc": "2.0", "id": rid, "result": self.handle(msg)}
        except InvalidParams as e:
            return err(-32602, str(e), rid)
        except LookupError as e:
            return err(-32601, "method not found: %s" % e, rid)
        except Exception as e:
            sys.stderr.write("ptmlink_mcp: %s in %s: %s\n" % (type(e).__name__, msg["method"], e))
            return err(-32603, "internal error in %s" % msg["method"], rid)

    def serve(self, inp=None, out=sys.stdout):
        def send(reply):
            out.write(json.dumps(reply, default=str) + "\n")
            out.flush()
        for line in (sys.stdin.buffer if inp is None else inp):
            if isinstance(line, bytes):   # stdin is read as bytes, so text that is not UTF-8 is an answer, not a crash
                try:
                    line = line.decode("utf-8")
                except UnicodeDecodeError:
                    send({"jsonrpc": "2.0", "id": None, "error": {"code": -32700, "message": "parse error: the message is not UTF-8 text"}})
                    continue
            line = line.strip()
            if not line:
                continue
            try:
                msg = json.loads(line)
            except ValueError:
                send({"jsonrpc": "2.0", "id": None, "error": {"code": -32700, "message": "parse error"}})
                continue
            if isinstance(msg, list):   # a batch: one reply per request, in one array; none for notifications alone
                if not msg:
                    send({"jsonrpc": "2.0", "id": None, "error": {"code": -32600, "message": "an empty batch holds no request"}})
                    continue
                replies = [r for r in (self.answer(m) for m in msg) if r is not None]
                if replies:
                    send(replies)
                continue
            reply = self.answer(msg)
            if reply is not None:
                send(reply)


def main(argv=None):
    ap = argparse.ArgumentParser(description="PTMlink as an MCP server on stdio")
    ap.add_argument("--data-dir", default=os.environ.get("PTMLINK_DATA_DIR"), help="a local copy of the website's data directory (the one holding build.json)")
    ap.add_argument("--data-url", default=os.environ.get("PTMLINK_DATA_URL"), help="the data address to download from (default: the PTMlink website)")
    ap.add_argument("--release-dir", default=os.environ.get("PTMLINK_RELEASE_DIR"), help="the Parquet release, needed by ptm_query and ptm_list")
    ap.add_argument("--no-verify", action="store_true", help="accept the release directory without checking its checksum")
    a = ap.parse_args(argv)
    try:
        db = PTMlink(data_dir=a.data_dir, data_url=a.data_url, release_dir=a.release_dir, verify_release=not a.no_verify)
    except (PTMlinkError, ValueError, OSError) as e:   # a wrong folder is one line, before any client connects
        sys.stderr.write("ptmlink_mcp: cannot start: %s\n" % e)
        return 2
    db.echo_messages = False   # messages go back to the agent inside the tool result; stdout carries JSON-RPC only
    Server(db).serve()
    return 0


if __name__ == "__main__":
    sys.exit(main())
