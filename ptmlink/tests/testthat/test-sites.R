# ptm_sites: every served site with its identifier and canonical position, and the modeled position
# where the structure sidecar exists
test_that("ptm_sites lists both sides with identifiers and reads the structure sidecar", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d)
  s <- ptm_sites("EGFR")
  expect_equal(nrow(s), 4L); expect_equal(sort(s$site), c("EGFR_N196", "EGFR_N352", "EGFR_N603", "EGFR_T693"))
  expect_equal(s$row_id[s$site == "EGFR_T693"], "PTMP:P00533:T693"); expect_equal(s$side[s$site == "EGFR_T693"], "phospho")
  expect_equal(s$position_canonical[s$site == "EGFR_N603"], 603L); expect_equal(s$residue[s$site == "EGFR_T693"], "T")
  expect_equal(s$frame_offset[s$site == "EGFR_T693"], 0L); expect_true(is.na(s$frame_offset[s$site == "EGFR_N603"]))
  expect_true(s$glyco_site_ambiguous[s$site == "EGFR_N352"])
  # the sidecar: a straight-line trace at 3.8 A per residue, confident to 650, disordered after
  expect_true(all(s$modelled)); expect_true(all(s$residue_matches_model)); expect_equal(s$model[1], "AF-P00533-F1 v6"); expect_equal(s$model_coverage[1], 1)
  expect_equal(s$plddt[s$site == "EGFR_N603"], 90); expect_equal(s$plddt[s$site == "EGFR_T693"], 40)
  expect_equal(s$x[s$site == "EGFR_N196"], 196 * 3.8); expect_equal(s$y[s$site == "EGFR_N196"], 0)
  # no sidecar for TFRC: structure columns NA, the rest complete
  expect_message(t <- ptm_sites("TFRC"), "no AlphaFold model"); expect_equal(nrow(t), 2L); expect_true(all(is.na(t$plddt))); expect_true(all(is.na(t$modelled))); expect_equal(t$row_id[t$side == "glyco"], "PTMG:P02786:N727")
  both <- suppressMessages(ptm_sites(c("EGFR", "P02786"))); expect_equal(nrow(both), 6L); expect_equal(unique(both$query), c("EGFR", "P02786"))
  expect_message(expect_null(ptm_sites("NOTAGENE")), "not a recognized symbol")
  off <- ptm_sites("EGFR", structure = FALSE); expect_true(all(is.na(off$plddt)))
})
