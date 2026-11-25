# ptm_network() reproduces the website's edge table on the shared fixture, row for row and column for column
test_that("the network roll-up matches the site's on the shared fixture", {
  sh <- jsonlite::fromJSON(testthat::test_path("fixtures", "net_shard.json"), simplifyVector = TRUE)
  rows <- sh$rows
  cols <- c("gene", "partner", "role", "rank", "drawn", "n_site_pairs", "n_own_sites", "n_partner_sites",
            "n_site_pairs_positive", "n_site_pairs_negative", "sign_mixed", "evidence_best", "evidence_median", "pctile_at_best",
            "rho_at_best", "k_at_best", "n_at_best", "I2_at_best", "q_at_best", "q_anchored_at_best",
            "rho_prot", "shadow_clean_at_best", "self_coupling_at_best", "self_glyco_at_best", "self_phospho_at_best",
            "n_site_pairs_anchored", "anch_any", "anch_dose", "anch_gene_editing", "anch_enzymatic",
            "own_site_at_best", "best_glyco_site", "best_phospho_site")
  res <- ptmlink:::.ptm_network_edges(rows, "G", ptm_default_filters(), 4L, "Q", cols)
  e <- res$edges
  expect_equal(names(e), cols); expect_equal(nrow(e), 5L)
  want <- utils::read.delim(testthat::test_path("fixtures", "net_edges_expected.tsv"), stringsAsFactors = FALSE)
  for (cn in names(want)) expect_equal(e[[cn]], want[[cn]], info = cn, tolerance = 1e-12)
  expect_equal(e$drawn, c(TRUE, TRUE, TRUE, TRUE, TRUE))   # the cap of four plus the pinned Q
  expect_equal(res$withheld[["cis"]], 1L); expect_equal(res$withheld[["k"]], 1L); expect_equal(res$n_cis, 0L)
  res2 <- ptmlink:::.ptm_network_edges(rows, "G", modifyList(ptm_default_filters(), list(include_cis = TRUE)), 25L, character(0), cols)
  expect_equal(res2$n_cis, 1L); expect_false(any(res2$edges$partner == "G"))
})

test_that("the composition filter reads the labels of the glycoforms the score used, as the site does", {
  sh <- jsonlite::fromJSON(testthat::test_path("fixtures", "net_shard.json"), simplifyVector = TRUE)
  f <- modifyList(ptm_default_filters(), list(hide_composition = TRUE))
  res <- ptmlink:::.ptm_filter(sh$rows, f)
  expect_equal(res$withheld[["composition"]], 1L)   # G_N1 x P1_S1: labeled on the used glycoform and on the phosphosite
  expect_false(any(res$rows$glyco_site == "G_N1" & res$rows$phospho_site == "P1_S1"))
  expect_true(any(res$rows$glyco_site == "G_N2" & res$rows$phospho_site == "P1_S2"))   # the label sits on an unused glycoform
})
