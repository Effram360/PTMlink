# Checks that the Python reader's connectivity excludes same-protein pairs.
"""Connectivity excludes same-protein pairs even when its sites differ."""
from pathlib import Path
import sys
import unittest
try:
    import pyarrow as pa
except ImportError:   # the whole-release tools need pyarrow; without it this test is skipped, not broken
    pa = None
# parents[2]: tests/serve/ -> repository root.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
from ptmlink import PTMlink


@unittest.skipIf(pa is None, "pyarrow is not installed")
class ConnectivityTest(unittest.TestCase):
    def test_same_protein_cross_site_exclusion(self):
        class Tiny(PTMlink):
            def _deg(self, tier, level, side, units):
                return [int(u == ("EGFR_N1" if side == "glyco" else "EGFR_S1")) for u in units]
            def _codes(self, tier, level, side, units):
                return [0 if u == ("EGFR_N1" if side == "glyco" else "EGFR_S1") else None for u in units]
        # Two trans rows: EGFR_N1--TFRC_S1 and TFRC_N1--EGFR_S1.
        tier = {"M": 2, "rows": {"gs": pa.array([0, 1]), "ps": pa.array([1, 0])}}
        result = Tiny.__new__(Tiny)._connectivity(tier, "site", ["EGFR_N1", "EGFR_S1"], 0)
        self.assertEqual(result["internal_pairs"], 0)
        self.assertEqual(result["expected_pairs"], 0)
        self.assertIsNone(result["ratio"])


if __name__ == "__main__":
    unittest.main()
