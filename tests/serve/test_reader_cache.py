# Checks the Python reader's remote cache using local file URLs.
"""Remote-reader cache regressions using local file URLs; no external requests."""
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr

# parents[2]: tests/serve/ -> repository root.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
from ptmlink import PTMlink, DataAbsent, DataUnreachable

ROOT = Path(__file__).resolve().parents[2]


class CacheTest(unittest.TestCase):
    def test_source_switch_correction_and_corruption(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            a, b = root / "a", root / "b"
            cache = root / "cache"

            def publish(path, value):
                path.mkdir(exist_ok=True)
                payload = json.dumps({"value": value}).encode()
                (path / "build.json").write_bytes(payload)
                (path / "SHA256SUMS.json").write_text(json.dumps({"build.json": hashlib.sha256(payload).hexdigest()}))

            def reader(path):
                return PTMlink(data_url=path.as_uri(), cache_dir=str(cache))

            publish(a, 1)
            publish(b, 2)
            first, second = reader(a), reader(b)
            self.assertEqual(first._json("build.json")["value"], 1)
            self.assertEqual(second._json("build.json")["value"], 2)
            self.assertNotEqual(first.cache_dir, second.cache_dir)
            publish(a, 3)
            self.assertEqual(reader(a)._json("build.json")["value"], 3)
            # Current instances remain consistent with their original session snapshot.
            self.assertEqual(first._json("build.json")["value"], 1)
            third = reader(a)
            Path(third._path("build.json")).write_text("corrupt cache")
            self.assertEqual(reader(a)._json("build.json")["value"], 3)
            (a / "build.json").write_text("corrupt upstream")
            Path(third._path("build.json")).unlink()
            with self.assertRaisesRegex(RuntimeError, "checksum mismatch"):
                reader(a)._json("build.json")
            with self.assertRaisesRegex(RuntimeError, "not in this data source"):
                reader(a)._path("unlisted.json")
            with self.assertRaises(DataAbsent):
                reader(a)._path("unlisted.json")


def _publish(path, files):
    """Write files {relative path: bytes} and the SHA256SUMS.json that lists them, as the website does."""
    path.mkdir(parents=True, exist_ok=True)
    sums = {}
    for rel, payload in files.items():
        (path / rel).parent.mkdir(parents=True, exist_ok=True)
        (path / rel).write_bytes(payload)
        sums[rel] = hashlib.sha256(payload).hexdigest()
    (path / "SHA256SUMS.json").write_text(json.dumps(sums))


class SourceTest(unittest.TestCase):
    """What a failure is called: absent (the source does not carry the file) versus unreachable (the host)."""

    def test_unreachable_host_is_not_absence(self):
        with tempfile.TemporaryDirectory() as td:
            db = PTMlink(data_url=(Path(td) / "nowhere").as_uri(), cache_dir=str(Path(td) / "cache"))
            with self.assertRaisesRegex(DataUnreachable, "could not reach"):
                db.pair("IL6ST_N131", "TMX1_S270")

    def test_cached_manifest_is_used_offline_with_a_warning(self):
        with tempfile.TemporaryDirectory() as td:
            src, cache = Path(td) / "src", Path(td) / "cache"
            _publish(src, {"build.json": b'{"release": "x"}', "cohorts.json": b"[]"})
            PTMlink(data_url=src.as_uri(), cache_dir=str(cache))._json("build.json")
            (src / "SHA256SUMS.json").unlink()
            (src / "cohorts.json").unlink()
            err = io.StringIO()
            with redirect_stderr(err):
                db = PTMlink(data_url=src.as_uri(), cache_dir=str(cache))
                self.assertEqual(db._json("build.json")["release"], "x")
            self.assertIn("checksum list cached", err.getvalue())
            with self.assertRaises(DataUnreachable):   # listed but never cached: the host is needed, and it is gone
                db._json("cohorts.json")

    def test_structure_through_bundles_then_per_gene(self):
        rec = {"gene": "EGFR", "first": 1, "last": 3, "n": 3, "seq_len": 3, "seq": "ANT", "plddt": [90, 80, 70], "xyz": [[0, 0, 0], [1, 0, 0], [2, 0, 0]]}
        with tempfile.TemporaryDirectory() as td:
            src = Path(td) / "src"
            _publish(src, {"structure/bundles.json": json.dumps({"format": 1, "genes": {"EGFR": "structure/bundles/0000.json.gz"}}).encode(),
                           "structure/bundles/0000.json.gz": gzip.compress(json.dumps({"EGFR": rec}).encode())})
            db = PTMlink(data_url=src.as_uri(), cache_dir=str(Path(td) / "cache"))
            self.assertEqual(db._structure("EGFR")["plddt"], [90, 80, 70])
            self.assertIsNone(db._structure("TFRC"))
            local = Path(td) / "local"
            (local / "structure").mkdir(parents=True)
            (local / "structure" / "EGFR.json.gz").write_bytes(gzip.compress(json.dumps(rec).encode()))
            self.assertEqual(PTMlink(data_dir=str(local))._structure("EGFR")["plddt"], [90, 80, 70])
            self.assertIsNone(PTMlink(data_dir=str(local))._structure("TFRC"))

    def test_command_line_errors_are_one_line_with_a_nonzero_exit(self):
        with tempfile.TemporaryDirectory() as td:
            env = dict(os.environ, HOME=td)
            r = subprocess.run([sys.executable, str(ROOT / "python" / "ptmlink.py"), "--data-url", (Path(td) / "nowhere").as_uri(), "pair", "IL6ST_N131", "TMX1_S270"],
                               capture_output=True, text=True, env=env)
            self.assertEqual(r.returncode, 1)
            self.assertEqual(r.stdout, "")
            self.assertEqual(len(r.stderr.strip().splitlines()), 1, r.stderr)
            self.assertIn("could not reach", r.stderr)
            self.assertNotIn("Traceback", r.stderr)


class SiteKeyTest(unittest.TestCase):
    def test_site_keys_keep_the_release_letter_case(self):
        # C10orf71 and its like have lowercase letters; a typed key must reach them whatever its case
        with tempfile.TemporaryDirectory() as td:
            r = PTMlink(data_dir=td)
            r._cache["index.json.gz"] = {"genes": [{"gene": "EGFR"}, {"gene": "C10orf71"}], "aliases": []}
            self.assertEqual(r._canon_site(" c10ORF71_s677 "), "C10orf71_S677")
            self.assertEqual(r._canon_site("C10orf71_S677"), "C10orf71_S677")
            self.assertEqual(r._canon_site("egfr_n603"), "EGFR_N603")
            self.assertEqual(r._canon_site("nope_s1"), "NOPE_S1")                      # unknown gene: upper case, as before
            self.assertEqual(r._canon_site("ptmg:p00533:n603"), "PTMG:P00533:N603")    # an identifier is not a site key
            self.assertEqual(r._canon_site("EGFR"), "EGFR")
            self.assertEqual(r._canon_gene("c10orf71"), "C10orf71")
            self.assertEqual(r._site_key("c10orf71_s677"), {"gene": "C10orf71", "key": "C10orf71_S677", "side": "phospho"})

    def test_one_filtered_reason_keeps_its_prefix_and_an_unreported_pair_is_named(self):
        import ptmlink
        self.assertEqual(ptmlink._withheld_units({"cis": 1}), ["Hidden by filters: 1 same-protein."])
        two = ptmlink._withheld_units({"cis": 1, "k": 2})
        self.assertTrue(two[0].startswith("Hidden by filters: ") and two[-1].endswith(".") and len(two) == 2)
        self.assertEqual(two[-1], "2 single-cohort.")                                   # the label follows the cohort minimum in force
        self.assertEqual(ptmlink._withheld_units({"k": 2}, min_k=6), ["Hidden by filters: 2 in fewer than 6 cohorts."])
        self.assertIn("Not a reported pair", ptmlink.text_explain([]))


class McpAnswerTest(unittest.TestCase):
    def test_the_list_answer_is_capped_and_cut_to_key_columns(self):
        # an assistant cannot read a table of 25,000 partner sites: every table of the list answer is capped, and the
        # site-pair tables keep the key columns only
        import ptmlink_mcp
        row = {"pair_id": "a|b", "glyco_site": "A_N1", "phospho_site": "B_S2", "rho_combined": 0.5, "sum_weight": 100.0,
               "evidence": 9.0, "cohorts_rho": "x" * 500, "site_named_glyco": True}
        res = {"internal": [dict(row) for _ in range(3)], "same_protein": [], "neighbourhood": [dict(row) for _ in range(7)],
               "partner_sites": [{"site": "S%d" % i} for i in range(50)], "convergent": [{"partner": "P"}], "note": {"a": 1}}
        server = ptmlink_mcp.Server.__new__(ptmlink_mcp.Server)   # shaping a result needs no data source
        out = server.shape("ptm_list", res, 5, {})
        self.assertEqual((len(out["partner_sites"]), out["n_partner_sites_total"]), (5, 50))
        self.assertEqual((len(out["neighbourhood"]), out["n_neighbourhood_total"]), (5, 7))
        self.assertNotIn("cohorts_rho", out["internal"][0])
        self.assertIn("rho_lo", out["internal"][0])
        self.assertTrue(out["internal"][0]["site_named_glyco"])
        self.assertEqual(out["partner_sites"], res["partner_sites"][:5])
        self.assertEqual(out["note"], {"a": 1})
        self.assertEqual(len(res["partner_sites"]), 50)                                   # the reader's own result is never cut
        # a ranked gene result keeps key columns and leaves the gene's site table to ptm_sites; a query with chosen
        # columns keeps exactly those; a network link keeps its key columns; a plain table is capped with its count
        gene = server.shape("ptm_gene", {"rows": [dict(row) for _ in range(4)], "sites": [{"site": "A_N1"}], "withheld": {"cis": 1}}, 2, {})
        self.assertEqual((len(gene["rows"]), gene["n_rows_total"]), (2, 4))
        self.assertNotIn("sites", gene); self.assertNotIn("cohorts_rho", gene["rows"][0]); self.assertIn("row_columns", gene)
        chosen = server.shape("ptm_query", {"rows": [{"glyco_site": "A_N1", "cohorts_rho": "x"}]}, 5, {"columns": ["cohorts_rho"]})
        self.assertEqual(chosen["rows"], [{"glyco_site": "A_N1", "cohorts_rho": "x"}])
        net = server.shape("ptm_network", {"edges": [{"partner": "B", "evidence_best": 3.0, "n_own_sites": 2}]}, 5, {})
        self.assertEqual(net["edges"], [{"partner": "B", "evidence_best": 3.0}])
        self.assertEqual(server.shape("ptm_sites", [{"site": "S%d" % i} for i in range(3)], 2, {}), {"rows": [{"site": "S0"}, {"site": "S1"}], "n_rows_total": 3})

    def test_arguments_are_checked_against_the_tool_description(self):
        # a value of the wrong kind is refused in words, never passed on to fail as a Python exception
        import ptmlink_mcp
        text = {"anyOf": [{"type": "string"}, {"type": "array", "items": {"type": "string"}}]}
        self.assertTrue(ptmlink_mcp._type_ok("EGFR", text) and ptmlink_mcp._type_ok(["EGFR", "TFRC"], text))
        self.assertFalse(ptmlink_mcp._type_ok(123, text) or ptmlink_mcp._type_ok(["EGFR", 1], text))
        self.assertFalse(ptmlink_mcp._type_ok(True, {"type": "integer"}) or ptmlink_mcp._type_ok("5", {"type": "integer"}))
        self.assertTrue(ptmlink_mcp._type_ok(0.5, {"type": "number"}) and ptmlink_mcp._type_ok(False, {"type": "boolean"}))
        self.assertEqual(ptmlink_mcp._type_words(text), "text, or a list of text")
        # the reader's command-line wording is turned to tool names for an assistant
        self.assertEqual(ptmlink_mcp.for_mcp("X is not a recognized symbol, alias or UniProt accession; status X says how symbols resolve"),
                         "X is not a recognized symbol, alias or UniProt accession; ptm_status says how symbols resolve")
        self.assertIn("not available here", ptmlink_mcp.for_mcp("query and list read the full table of reported pairs (site_served.parquet), which is released with "
                                                                 "the published article; once you have it, give its release directory with --release-dir (release_dir= in Python)."))


class FilterCheckTest(unittest.TestCase):
    def test_filters_are_completed_and_checked(self):
        import ptmlink
        f = ptmlink.check_filters({"min_k": 1})
        self.assertEqual(f["min_k"], 1); self.assertTrue(f["significant_only"]); self.assertEqual(set(f), set(ptmlink.DEFAULT_FILTERS))
        self.assertEqual(ptmlink.check_filters(None), ptmlink.DEFAULT_FILTERS)
        for bad, word in (({"min_kk": 1}, 'unknown filter "min_kk"'), ({"min_k": "two"}, "min_k is a whole number"), ({"significant_only": "no"}, "true or false"),
                          ({"max_i2": 5}, "between 0 and 1"), ("glyco", "named settings")):
            with self.assertRaises(ValueError) as cm:
                ptmlink.check_filters(bad)
            self.assertIn(word, str(cm.exception))
        self.assertIsNone(ptmlink.check_filters({"max_i2": None})["max_i2"])

    def test_names_are_checked_before_anything_is_read(self):
        import ptmlink
        self.assertEqual(ptmlink._names(" EGFR ", "gene", "x"), ["EGFR"])
        for bad, word in ((123, "gene is int, not text"), ([], "gene is empty"), (["EGFR", None], "a missing entry"), ([""], "missing or empty entry"), (None, "gene is missing")):
            with self.assertRaises(ValueError) as cm:
                ptmlink._names(bad, "gene", "a gene symbol")
            self.assertIn(word, str(cm.exception))


if __name__ == "__main__":
    unittest.main()


def _tiny(td):
    """A reader over a two-gene data source held in memory: no file is read, nothing is downloaded."""
    def row(gs, ps, ev, q=1e-6, k=5, **kw):
        gg, pg = gs.rsplit("_", 1)[0], ps.rsplit("_", 1)[0]
        return dict({"glyco_site": gs, "phospho_site": ps, "glyco_gene": gg, "phospho_gene": pg, "is_cis": gg == pg, "rho_combined": 0.3, "sum_weight": 100.0,
                     "evidence": ev, "q_value": q, "k_cohorts": k, "I2": 0.1, "cohorts_rho": "BRCA:0.3:60;OV:0.3:50", "anch_any": False, "rho_prot": 0.1,
                     "self_glyco": 0.2, "self_phospho": 0.3, "self_coupling": 0.3}, **kw)
    rows = [row("EGFR_N603", "TFRC_S34", 9.0), row("EGFR_N603", "COIL_T303", 7.0), row("EGFR_N196", "TFRC_S34", 5.0), row("TFRC_N727", "COIL_T303", 12.0)]
    r = PTMlink(data_dir=td)
    genes = [{"gene": g, "status_glyco": "usable", "status_phospho": "usable", "n_served_as_glyco": sum(x["glyco_gene"] == g for x in rows),
              "n_served_as_phospho": sum(x["phospho_gene"] == g for x in rows)} for g in ("EGFR", "TFRC", "COIL")]
    r._cache["index.json.gz"] = {"genes": genes, "aliases": [
        {"query": "ERBB1", "resolves_to": "EGFR", "kind": "alias"}, {"query": "P00533", "resolves_to": "EGFR", "kind": "accession"},
        {"query": "CD71", "resolves_to": "TFRC", "kind": "alias"}, {"query": "CD71", "resolves_to": "TFRC", "kind": "previous symbol"},
        {"query": "G2", "resolves_to": "EGFR", "kind": "alias"}, {"query": "G2", "resolves_to": "TFRC", "kind": "alias"}]}
    r._cache["build.json"] = {"release": "fixture"}
    r._cache["cohorts.json"] = []
    for g in ("EGFR", "TFRC", "COIL"):
        r._cache["genes/%s.json.gz" % g] = {"gene": g, "rows": [dict(x) for x in rows if g in (x["glyco_gene"], x["phospho_gene"])], "glyco_sites": [], "phospho_sites": []}
    r.echo_messages = False
    return r


class NameTest(unittest.TestCase):
    def test_an_alias_listed_twice_resolves_and_an_ambiguous_one_says_what_to_pass(self):
        with tempfile.TemporaryDirectory() as td:
            r = _tiny(td)
            self.assertEqual(r._resolve("cd71"), {"gene": "TFRC", "via": "alias of TFRC"})
            with self.assertRaises(ValueError) as cm:
                r.gene("G2")
            self.assertEqual(str(cm.exception), "G2 is ambiguous: it is an alias of EGFR and TFRC; give one of these symbols")
            st = r.status(["EGFR", "G2", "NOPE", "EGFR_N603"])   # one row per query: an ambiguous alias does not lose the others
            self.assertEqual([s["gene"] for s in st], ["EGFR", None, None, None])
            self.assertEqual(st[1]["status_glyco"], "ambiguous alias"); self.assertTrue(st[1]["resolved_via"].startswith("ambiguous: it is an alias of EGFR and TFRC"))
            self.assertIn("a site key", st[3]["resolved_via"])
            both = r.gene(["EGFR", "G2"])
            self.assertEqual(both["genes"], ["EGFR"]); self.assertIn("G2 skipped: it is an alias of EGFR and TFRC", r.messages[-1])

    def test_a_name_given_twice_is_read_once_and_a_site_key_may_use_an_alias(self):
        with tempfile.TemporaryDirectory() as td:
            r = _tiny(td)
            one = r.gene("EGFR")
            many = r.gene(["EGFR", "egfr", "ERBB1", "P00533"])
            self.assertEqual(many["genes"], ["EGFR"]); self.assertEqual(len(many["rows"]), len(one["rows"]))
            self.assertEqual(len(r.site(["EGFR_N603", "egfr_n603", "ERBB1_N603"])["rows"]), len(r.site("EGFR_N603")["rows"]))
            self.assertEqual(r.site("ERBB1_N603")["site"], "EGFR_N603"); self.assertEqual(r._canon_site("p00533_n603"), "EGFR_N603")
            with self.assertRaises(ValueError) as cm:
                r.gene(["NOPE1", "NOPE2"])
            self.assertIn("none of the 2 names is a gene", str(cm.exception))

    def test_a_typed_name_never_leaves_the_data_source(self):
        with tempfile.TemporaryDirectory() as td:
            data = Path(td) / "data"; (data / "genes").mkdir(parents=True)
            with gzip.open(Path(td) / "EVIL.json.gz", "wt") as fh:
                json.dump({"gene": "EVIL", "rows": []}, fh)
            r = PTMlink(data_dir=str(data))
            for rel in ("genes/../../EVIL.json.gz", "/etc/hostname", "genes\\..\\x"):
                with self.assertRaises(DataAbsent):
                    r._path(rel)

    def test_pairs_are_asked_element_by_element_and_an_unreported_one_is_said(self):
        with tempfile.TemporaryDirectory() as td:
            r = _tiny(td)
            for call in (r.pair, r.explain):
                with self.assertRaises(ValueError) as cm:
                    call(["EGFR_N603", "EGFR_N196", "TFRC_N727"], ["TFRC_S34", "COIL_T303"])
                self.assertIn("pair up element by element", str(cm.exception))
            with self.assertRaises(ValueError) as cm:
                r.pair("PTMG:P00533:N603|PTMP:P02786:S34", "TFRC_S34")
            self.assertIn("leave phospho_site out", str(cm.exception))
            with self.assertRaises(ValueError):
                r.explain([], [])
            with self.assertRaises(ValueError) as cm:
                r.pair("", "TFRC_S34")
            self.assertIn("glyco_site holds a missing or empty entry", str(cm.exception))
            r.messages.clear()
            e = r.explain(["EGFR_N603", "TFRC_N727"], ["TFRC_S34", "TFRC_S34"])   # the second is a same-protein pair that is not reported
            self.assertEqual({x["glyco_site"] for x in e}, {"EGFR_N603"})
            self.assertTrue(any("TFRC_N727 x TFRC_S34 is not a reported pair" in m for m in r.messages))
            absent = r.pair("TFRC_N727", "TFRC_S34")
            self.assertEqual((absent["glyco_site"], absent["phospho_site"], absent["row"]), ("TFRC_N727", "TFRC_S34", None))

    def test_the_own_protein_rule_names_the_larger_absolute_side_and_a_missing_value_is_said(self):
        import ptmlink
        with tempfile.TemporaryDirectory() as td:
            r = _tiny(td)
            row = dict(r.gene("EGFR")["rows"][0], self_glyco=0.46, self_phospho=-0.68, self_coupling=0.68, rho_prot=None)
            e = {x["rule"]: x for x in r.explain(rows=[row])}
            self.assertIn("The phosphosite runs against its protein, at -0.68.", e["self_coupling"]["statement"])
            self.assertIsNone(e["protein_shadow"]["applies"]); self.assertEqual(e["protein_shadow"]["statement"], "Not assessed: this pair has no value for rho_prot.")
            text = ptmlink.text_explain(list(e.values()))
            self.assertNotIn("parent proteins co-vary", text); self.assertIn("Not assessed, no value: parent-protein co-variation", text)
            row.update(self_glyco=0.75, self_coupling=0.75)
            self.assertIn("The glycosite tracks its protein at 0.75.", {x["rule"]: x for x in r.explain(rows=[row])}["self_coupling"]["statement"])

    def test_a_condition_value_has_its_column_type(self):
        cv = PTMlink._condition_value
        self.assertEqual((cv("k", "int32", 5), cv("k", "int32", 5.0), cv("q", "double", 1), cv("a", "bool", True), cv("g", "string", "EGFR")), (5, 5, 1.0, True, "EGFR"))
        self.assertEqual((cv("k", "int32", "5", text=True), cv("q", "double", "1e-10", text=True), cv("a", "bool", "TRUE", text=True)), (5, 1e-10, True))   # typed on the command line
        for args, word in ((("k_cohorts", "int32", 2.5), "k_cohorts holds whole numbers; got 2.5"), (("q_value", "double", "small"), "q_value holds numbers"),
                           (("anch_any", "bool", "yes"), "anch_any holds true or false"), (("q_value", "double", None), "compares with a value, got none"),
                           (("glyco_gene", "string", 5), "glyco_gene holds text"), (("k_cohorts", "int32", True), "k_cohorts holds numbers")):
            with self.assertRaises(ValueError) as cm:
                cv(*args)
            self.assertIn(word, str(cm.exception))

    def test_a_reader_reads_the_checksum_list_again_when_the_site_was_updated(self):
        with tempfile.TemporaryDirectory() as td:
            site, cache = Path(td) / "site", Path(td) / "cache"
            site.mkdir()

            def publish(files):
                for name, value in files.items():
                    (site / name).write_bytes(json.dumps(value).encode())
                (site / "SHA256SUMS.json").write_text(json.dumps({n: hashlib.sha256((site / n).read_bytes()).hexdigest() for n in files}))
            publish({"build.json": {"v": 1}, "cohorts.json": [{"cohort": "A"}]})
            r = PTMlink(data_url=site.as_uri(), cache_dir=str(cache))
            self.assertEqual(r._json("build.json"), {"v": 1})
            publish({"build.json": {"v": 2}, "cohorts.json": [{"cohort": "B"}]})   # the site is updated while this reader lives
            self.assertEqual(r._json("cohorts.json"), [{"cohort": "B"}])           # the list is read again, once, and the new file accepted
            (site / "cohorts.json").write_bytes(b'[{"cohort": "tampered"}]')       # a file that no list names is still refused
            r2 = PTMlink(data_url=site.as_uri(), cache_dir=str(Path(td) / "cache2"))
            with self.assertRaises(Exception) as cm:
                r2._json("cohorts.json")
            self.assertIn("checksum mismatch", str(cm.exception))


class McpProtocolTest(unittest.TestCase):
    def _server(self, td):
        import ptmlink_mcp
        return ptmlink_mcp.Server(_tiny(td))

    def _talk(self, server, lines):
        out = io.StringIO()
        server.serve(inp=lines, out=out)
        return [json.loads(l) for l in out.getvalue().splitlines()]

    def test_messages_that_are_not_requests_get_the_right_answer_or_none(self):
        with tempfile.TemporaryDirectory() as td:
            s = self._server(td)
            ping = lambda i: json.dumps({"jsonrpc": "2.0", "id": i, "method": "ping"})
            r = self._talk(s, [b'\xff\xfe{"x":1}\n', ping(1) + "\n"])   # bytes that are not UTF-8 are answered, and the server lives on
            self.assertEqual((r[0]["error"]["code"], r[1]["result"]), (-32700, {}))
            r = self._talk(s, ["[%s, %s, %s]" % (ping(1), json.dumps({"jsonrpc": "2.0", "method": "notifications/initialized"}), ping("b")), "[]"])
            self.assertEqual([x["id"] for x in r[0]], [1, "b"]); self.assertEqual(r[1]["error"]["code"], -32600)   # a batch: one array of replies
            r = self._talk(s, [json.dumps({"jsonrpc": "2.0", "id": 22, "result": {}}), json.dumps({"jsonrpc": "2.0", "id": True, "method": "ping"}),
                               json.dumps({"id": 3, "method": "ping"}), json.dumps({"jsonrpc": "2.0", "id": 4, "method": ["ping"]}), ping(5)])
            self.assertEqual([x.get("error", {}).get("code") for x in r], [-32600, -32600, -32600, None])   # the response is not answered
            self.assertEqual(r[0]["id"], None)

    def test_tool_calls_are_checked_and_nothing_of_python_leaks(self):
        with tempfile.TemporaryDirectory() as td:
            s = self._server(td)
            call = lambda i, name, args: json.dumps({"jsonrpc": "2.0", "id": i, "method": "tools/call", "params": {"name": name, "arguments": args}})
            r = self._talk(s, [call(1, "ptm_nope", {}), call(2, ["ptm_gene"], {}), call(3, "ptm_gene", ["EGFR"]),
                               json.dumps({"jsonrpc": "2.0", "id": 4, "method": "tools/call", "params": ["ptm_gene"]}),
                               call(5, "ptm_gene", {"gene": "EGFR", "side": None}), call(6, "ptm_gene", {"gene": "EGFR", "limit": 10 ** 30}),
                               call(7, "ptm_pair", {"glyco_site": ["EGFR_N603"] * 26, "phospho_site": "TFRC_S34"}),
                               call(8, "ptm_gene", {"gene": ["EGFR", "TFRC"], "limit": 2}), call(9, "ptm_sites", {"x": "EGFR", "limit": 100}),
                               call(10, "ptm_status", {"gene": ["EGFR", "G2"]})])
            self.assertEqual([x.get("error", {}).get("code") for x in r[:4]], [-32602] * 4)   # an unknown tool, a name or arguments of the wrong shape
            self.assertFalse(r[4]["result"]["isError"])                                       # null is "not given"
            self.assertIn("limit is a whole number from 0 to 500", r[5]["result"]["content"][0]["text"])
            self.assertIn("at most 25 pairs", r[6]["result"]["content"][0]["text"])
            rows = json.loads(r[7]["result"]["content"][-1]["text"])["rows"]
            self.assertEqual([x["query_gene"] for x in rows], ["EGFR", "TFRC"])               # the top row of each gene, not the first gene alone
            self.assertEqual(set(json.loads(r[8]["result"]["content"][-1]["text"])), {"rows", "n_rows_total"})   # the same shape whether or not the cap cut it
            self.assertFalse(r[9]["result"]["isError"])
            for x in r:
                text = json.dumps(x)
                self.assertFalse(any(w in text for w in ("Traceback", "Error:", "NoneType", "object has no attribute", "unhashable")), text[:300])

    def test_descriptions_speak_json_and_a_link_keeps_its_drawn_mark(self):
        import ptmlink_mcp
        with tempfile.TemporaryDirectory() as td:
            s = self._server(td)
            text = json.dumps(s.tools) + s.instructions
            for word in ("default TRUE", "default FALSE", "default NA", "character(0)", "list()", "ptm_default_filters()", "to FALSE", "withheld$"):
                self.assertNotIn(word, text)
            props = {t["name"]: t["inputSchema"]["properties"] for t in s.tools}
            self.assertEqual(props["ptm_query"]["limit"]["default"], 25)
            self.assertIn({"type": "null"}, props["ptm_network"]["top_n"]["anyOf"]); self.assertIn({"type": "null"}, props["ptm_gene"]["filters"]["properties"]["max_i2"]["anyOf"])
            self.assertFalse(ptmlink_mcp._type_ok(["q_value", "<", 0.01, 5], props["ptm_query"]["where"]["items"]))   # a condition has at most three parts
            self.assertIn("drawn", ptmlink_mcp.EDGE_COLUMNS)
            net = s.shape("ptm_network", {"edges": [{"partner": "B", "drawn": True, "n_own_sites": 2}]}, 5, {})
            self.assertEqual(net["edges"], [{"partner": "B", "drawn": True}])
