# Names, values and shapes a caller may get wrong: each is read as meant or refused in plain words, never bent to fit.
fx <- function(parquet = FALSE) { d <- ptmlink:::.ptm_write_fixture(parquet = parquet); ptm_set_data_dir(d); assign("release_dir", NULL, envir = ptmlink:::.ptm_env); d }

test_that("an alias listed twice for one symbol resolves; a truly ambiguous one says what to pass", {
  fx()
  expect_equal(ptm_status("CD71")$gene, "TFRC"); expect_equal(ptm_status("cd71")$resolved_via, "alias of TFRC")
  expect_error(ptm_gene("G2"), "G2 is ambiguous: it is an alias of EGFR and TFRC; give one of these symbols", fixed = TRUE)
  st <- ptm_status(c("EGFR", "G2", "NOPE", "EGFR_N603"))
  expect_equal(nrow(st), 4L); expect_equal(st$gene, c("EGFR", NA, NA, NA))
  expect_equal(st$status_glyco[2], "ambiguous alias"); expect_match(st$resolved_via[2], "^ambiguous: it is an alias of EGFR and TFRC")
  expect_match(st$resolved_via[4], "a site key, not a gene"); expect_equal(st$resolved_via[3], "not a recognized symbol")
  # the same columns in the same order, whatever the first query is
  expect_equal(names(ptm_status(c("NOPE", "EGFR"))), names(ptm_status(c("EGFR", "NOPE"))))
  expect_message(g <- ptm_gene(c("EGFR", "G2"), filters = list(min_k = 1)), "G2 skipped: it is an alias of EGFR and TFRC")
  expect_equal(unique(g$query_gene), "EGFR")
  expect_message(expect_message(s <- ptm_sites(c("EGFR", "G2", "EGFR_N603"), structure = FALSE), "G2 skipped"), "EGFR_N603 is a site key, not a gene")
  expect_equal(unique(s$gene), "EGFR"); expect_error(ptm_sites("EGFR", structure = NA), "structure is TRUE or FALSE")
})

test_that("a gene or a site named twice, in any spelling, is read once", {
  fx()
  one <- ptm_gene("EGFR"); many <- ptm_gene(c("EGFR", "egfr", "ERBB1", "P00533"))
  expect_equal(nrow(many), nrow(one)); expect_equal(sum(attr(many, "withheld")), sum(attr(one, "withheld")))
  expect_equal(nrow(ptm_site(c("EGFR_N603", "egfr_n603", "ERBB1_N603"))), nrow(ptm_site("EGFR_N603")))
})

test_that("a site key may name its protein by an alias or accession", {
  fx()
  expect_equal(attr(ptm_site("ERBB1_N603"), "query"), "partners of EGFR_N603")
  expect_equal(ptm_pair("P02786_N727", "COIL_T303")$row$glyco_site, "TFRC_N727")
  expect_error(ptm_site("ERBB1_N9999"), "EGFR_N9999 is in no reported pair", fixed = TRUE)
})

test_that("a name typed by the caller never reads a file outside the data source", {
  d <- fx()
  writeLines("{}", file.path(dirname(d), "EVIL.json"))
  expect_error(ptm_site("../../EVIL_N1"), "no reported pair|not in this data source")
  expect_error(ptmlink:::.ptm_path("genes/../../EVIL.json"), class = "ptm_absent"); expect_error(ptmlink:::.ptm_path("/etc/hostname"), class = "ptm_absent")
})

test_that("pairs are asked element by element: equal lengths, or one of either", {
  fx()
  expect_error(ptm_pair(c("EGFR_N603", "EGFR_N196", "EGFR_N352"), c("EGFR_T693", "ALG3_S13")), "pair up element by element", fixed = TRUE)
  expect_error(ptm_explain(c("EGFR_N603", "EGFR_N196", "EGFR_N352"), c("EGFR_T693", "ALG3_S13")), "pair up element by element", fixed = TRUE)
  p <- ptm_pair("EGFR_N603", c("EGFR_T693", "ALG3_S13")); expect_equal(nrow(p$row), 2L)
  none <- ptm_pair(c("EGFR_N603", "TFRC_N727"), c("COIL_T303", "ALG3_S13")); expect_null(none$row)
  out <- capture.output(print(none)); expect_match(out[1], "Not reported: EGFR_N603 x COIL_T303, TFRC_N727 x ALG3_S13.", fixed = TRUE)
  g <- ptm_gene("EGFR")
  expect_error(ptm_explain("EGFR_N352", "XXX_S1", rows = g), "rows holds none of the pairs asked for", fixed = TRUE)
})

test_that("the own-protein rule names the side with the larger absolute correlation, with its sign", {
  fx()
  g <- as.data.frame(ptm_gene("EGFR", filters = list(min_k = 1, include_cis = TRUE)))
  r <- g[1, ]; r$self_glyco <- 0.46; r$self_phospho <- -0.68; r$self_coupling <- 0.68
  s <- ptm_explain(rows = r); s <- s$statement[s$rule == "self_coupling"]
  expect_match(s, "The phosphosite runs against its protein, at -0.68.", fixed = TRUE)
  r$self_glyco <- 0.75; r$self_coupling <- 0.75
  s <- ptm_explain(rows = r); expect_match(s$statement[s$rule == "self_coupling"], "The glycosite tracks its protein at 0.75.", fixed = TRUE)
})

test_that("a rule whose value the pair lacks says so, and is never printed as a finding", {
  fx()
  r <- as.data.frame(ptm_gene("EGFR"))[1, ]; r$rho_prot <- NA; r$I2 <- NA
  e <- ptm_explain(rows = r)
  expect_true(all(is.na(e$applies[e$rule %in% c("protein_shadow", "heterogeneity")])))
  expect_match(e$statement[e$rule == "protein_shadow"], "^Not assessed: this pair has no value for rho_prot")
  out <- paste(capture.output(print(e)), collapse = "\n")
  expect_false(grepl("The parent proteins co-vary", out, fixed = TRUE)); expect_false(grepl("The cohorts disagree", out, fixed = TRUE))
  expect_match(out, "Not assessed, no value: parent-protein co-variation, heterogeneity", fixed = TRUE)
  expect_match(out, "Co-variation only", fixed = TRUE)   # the rules that test nothing stay as context
})

test_that("a query's values are checked against their columns, and its genes resolve as everywhere else", {
  testthat::skip_if_not_installed("arrow")
  d <- fx(parquet = TRUE); ptm_set_release_dir(file.path(d, "release"))
  all_rows <- list(min_k = 1, include_cis = TRUE)
  expect_error(ptm_query(list(list("k_cohorts", ">=", 5.5))), "k_cohorts holds whole numbers; got 5.5", fixed = TRUE)
  expect_error(ptm_query(list(list("q_value", "<", "abc"))), "q_value holds numbers; got \"abc\"", fixed = TRUE)
  expect_error(ptm_query(list(list("q_value", "<", NA))), "compares with a value, got none", fixed = TRUE)
  expect_error(ptm_query(list(list("anch_any", "==", "yes"))), "anch_any holds TRUE or FALSE", fixed = TRUE)
  expect_error(ptm_query(list(list("glyco_gene", "==", 5))), "glyco_gene holds text", fixed = TRUE)
  expect_error(ptm_query(list(list("q_value", "<", 0.01, 5))), "a condition is list(column, operator, value)", fixed = TRUE)
  expect_equal(nrow(ptm_query(list(list("k_cohorts", ">=", 6)), filters = all_rows)), 4L); expect_equal(nrow(ptm_query(list(list("k_cohorts", "in", c(5, 9))), filters = all_rows)), 4L)
  # an alias and an accession are the gene; a name that is not one is said
  egfr <- nrow(ptm_query(genes = "EGFR", filters = all_rows))
  expect_equal(nrow(ptm_query(genes = "ERBB1", filters = all_rows)), egfr); expect_equal(nrow(ptm_query(genes = c("P00533", "egfr"), filters = all_rows)), egfr)
  expect_message(q <- ptm_query(genes = "NOTAGENE"), "NOTAGENE is not a recognized symbol"); expect_equal(nrow(q), 0L)
  expect_error(ptm_query(genes = "G2"), "G2 is ambiguous", fixed = TRUE)
  expect_equal(nrow(ptm_query(sites = "ERBB1_N603", filters = all_rows)), 2L)
  out <- capture.output(print(ptm_query(genes = "EGFR", limit = 0)))
  expect_true(any(grepl("No rows returned: the limit was 0", out, fixed = TRUE)))
  # conditions that match more rows than the reader will hold are refused, with the way through
  testthat::local_mocked_bindings(.PTM_QUERY_MAX_ROWS = 3L, .package = "ptmlink")
  expect_error(ptm_query(list(list("k_cohorts", ">=", 1))), "the conditions match 6 reported pairs, more than 3", fixed = TRUE)
  expect_equal(nrow(ptm_query(list(list("k_cohorts", ">=", 1)), filters = all_rows, allow_full = TRUE)), 6L)
  expect_equal(nrow(ptm_query(list(list("k_cohorts", ">=", 9)), filters = all_rows)), 2L)
})

test_that("a list's site entry shows its own pairs, and a site in no reported pair is said", {
  testthat::skip_if_not_installed("arrow")
  d <- fx(parquet = TRUE); ptm_set_release_dir(file.path(d, "release"))
  expect_message(l <- ptm_list(c("EGFR_N603", "EGFR_N9999", "TFRC_N727"), n_perm = 0), "EGFR_N9999 is in no reported pair of this release", fixed = TRUE)
  expect_equal(l$members$n_reported, c(2L, 0L, 2L))
  out <- capture.output(print(l)); expect_true(any(grepl("reported pairs", out, fixed = TRUE))); expect_false(any(grepl("glyco pairs", out, fixed = TRUE)))
  one <- paste(capture.output(print(ptm_list(c("EGFR", "egfr", "P00533"), n_perm = 0))), collapse = " ")
  expect_match(gsub(" +", " ", one), "it holds one protein, and pairs within a protein are not counted", fixed = TRUE)
  expect_error(ptm_list("EGFR", seed = "a"), "seed is one whole number", fixed = TRUE)
})

test_that("print methods refuse a row count that is not one, and stacked networks name every gene", {
  fx()
  g <- ptm_gene("EGFR")
  expect_error(print(g, n = 0), "n is the number of rows to show", fixed = TRUE); expect_error(print(g, n = -1), "n is the number of rows to show", fixed = TRUE)
  both <- rbind(ptm_network("EGFR"), ptm_network("TFRC"))
  out <- capture.output(print(both))
  expect_match(out[1], "^Network of EGFR, TFRC: "); expect_match(out[1], paste0(sum(both$n_site_pairs), " site pair"))
  net <- ptm_network("EGFR", top_n = 0)
  expect_error(ptm_plot_network(net), "marks no partner for drawing", fixed = TRUE)
  expect_error(ptm_plot_network(ptm_network("EGFR"), max_edges = 2.5), "max_edges is a whole number", fixed = TRUE)
  expect_error(ptm_write_network(net, NA), "path is the file name to write", fixed = TRUE)
  expect_error(ptm_write_network(net, file.path(tempdir(), "no_such_folder_here", "x")), "no such directory", fixed = TRUE)
})
