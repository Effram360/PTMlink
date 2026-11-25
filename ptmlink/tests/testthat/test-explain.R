# ptm_explain: one row per rule per pair, every threshold read from inst/reference.json
test_that("ptm_explain reads every rule from the reference and judges the fixture's pairs as the site would", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d); ptm_set_release_dir(file.path(d, "release"), verify = FALSE)   # the feature tables complete the partner ids
  rules <- vapply(ptm_reference()$rules, function(r) r$id, character(1))
  e <- ptm_explain("TFRC_N727", "COIL_T303")
  expect_equal(e$rule, rules); expect_equal(unique(e$pair_id), "PTMG:P02786:N727|PTMP:P38432:T303")
  ap <- setNames(e$applies, e$rule)
  expect_true(ap[["direction"]]); expect_true(ap[["significance"]]); expect_true(ap[["perturbation"]]); expect_false(ap[["same_protein"]])
  expect_false(ap[["protein_shadow"]]); expect_true(ap[["self_coupling"]]); expect_false(ap[["tissue_composition"]])   # the glycosite is labeled on the used glycoform, COIL_T303 is not
  expect_true(is.na(ap[["depth"]])); expect_false(ap[["heterogeneity"]]); expect_false(ap[["one_cohort"]]); expect_true(is.na(ap[["without_one_cohort"]]))
  expect_true(ap[["adjusted"]]); expect_false(ap[["glycoform_pooling"]]); expect_false(ap[["sequon_ambiguity"]])
  st <- setNames(e$statement, e$rule)
  expect_match(st[["perturbation"]], "gene editing"); expect_match(st[["self_coupling"]], "phosphosite tracks its protein at 0.77"); expect_match(st[["significance"]], "130 tumors")
  expect_equal(e$threshold[e$rule == "protein_shadow"], "abs>= 0.15"); expect_equal(e$filter[e$rule == "self_coupling"], "hide_self_coupled")
  expect_equal(e$value[e$rule == "heterogeneity"], "0.44")
  # a cis pair with a missing self-coupling on another row, and the sequon flag from the site record
  e2 <- ptm_explain("EGFR_N352", "SLC20A2_S256"); ap2 <- setNames(e2$applies, e2$rule)
  expect_true(is.na(ap2[["self_coupling"]])); expect_true(ap2[["sequon_ambiguity"]]); expect_true(is.na(ap2[["adjusted"]]))
  e4 <- ptm_explain("TFRC_N727", "TFRC_S34"); expect_true(setNames(e4$applies, e4$rule)[["tissue_composition"]])   # labeled on both sides
  e3 <- ptm_explain("EGFR_N603", "EGFR_T693"); expect_true(setNames(e3$applies, e3$rule)[["same_protein"]]); expect_true(setNames(e3$applies, e3$rule)[["glycoform_pooling"]])
  # several pairs, and rows already in hand
  both <- ptm_explain(c("TFRC_N727", "EGFR_N603"), c("COIL_T303", "EGFR_T693")); expect_equal(nrow(both), 2L * length(rules))
  g <- ptm_gene("EGFR", "glyco"); r <- ptm_explain(rows = g); expect_equal(nrow(r), nrow(g) * length(rules)); expect_equal(unique(r$glyco_site), unique(g$glyco_site))
  expect_message(expect_null(ptm_explain("EGFR_N603", "NOPE_S1")), "not among the reported pairs")
  expect_error(ptm_explain(), "takes a glycosite and a phosphosite")
})
