# Printed summaries must fit the console: at 80 columns no line is wider, and the objects keep every column.
test_that("every print method fits an 80-column console and leaves the object whole", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d)
  old <- options(width = 80); on.exit(options(old))
  g <- ptm_gene("EGFR", "glyco"); p <- ptm_pair("TFRC_N727", "COIL_T303")
  shown <- list(ptm_status(c("EGFR", "NOTAGENE")), g, ptm_site("EGFR_T693"), p, ptm_explain("TFRC_N727", "COIL_T303"),
                ptm_network("EGFR"), ptm_release())
  for (x in shown) {
    out <- utils::capture.output(print(x))
    expect_gt(length(out), 0L)
    expect_true(all(nchar(out) <= 80L), info = paste(class(x)[1], out[nchar(out) > 80L], collapse = "\n"))
    expect_false(any(grepl(" $", out)), info = class(x)[1])
  }
  expect_true(all(c("glyco_site", "phospho_site", "evidence", "q_value") %in% names(g)))
  expect_identical(utils::capture.output(print(g, n = 1))[1], utils::capture.output(print(g))[1])
})

test_that("a result cut down to other columns prints as a plain data.frame, never empty and never an error", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d)
  e <- ptm_explain("TFRC_N727", "COIL_T303")
  out <- utils::capture.output(print(e[e$applies %in% TRUE, c("rule", "value", "threshold", "statement")]))
  expect_true(any(grepl("significance", out))); expect_false(any(grepl("^Full wording", out)))
  net <- ptm_network("EGFR")
  out <- utils::capture.output(expect_no_error(print(head(net[net$drawn, c("partner", "role", "n_site_pairs", "evidence_best", "rho_at_best")]))))
  expect_true(any(grepl("VDAC2", out)))
  st <- ptm_status("EGFR"); out <- utils::capture.output(print(st[, c("query", "gene")])); expect_true(any(grepl("EGFR", out)))
  g <- ptm_gene("EGFR", "glyco"); out <- utils::capture.output(print(g[, c("glyco_site", "phospho_site")])); expect_true(any(grepl("VDAC2_S115", out)))
})

test_that("a long table cut down to other columns shows its first rows and says how many are left", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d)
  g <- ptm_gene("EGFR", "glyco"); long <- g[rep(seq_len(nrow(g)), length.out = 40L), c("glyco_site", "phospho_site", "evidence")]
  out <- utils::capture.output(print(long))
  expect_lte(length(out), 12L); expect_true(any(grepl("^30 more rows: print\\(x, n = 50\\)", out)))
  expect_length(utils::capture.output(print(long, n = 50)), 41L)          # a header line and all 40 rows, no footer
  expect_equal(nrow(as.data.frame(long)), 40L)                            # the object itself is never cut
})

test_that("a summary leaves out a field whose columns were not asked for, and counts one cohort as one", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d)
  g <- ptm_gene("EGFR", "glyco")
  narrow <- g[, c("glyco_site", "phospho_site", "rho_combined", "evidence", "q_value", "k_cohorts")]
  out <- utils::capture.output(print(narrow))
  expect_false(any(grepl("perturbation|95% CI", out))); expect_true(any(grepl("VDAC2_S115", out)))
  full <- utils::capture.output(print(g)); expect_true(any(grepl("perturbation", full)))
  st <- ptm_status("EGFR"); st$n_cohorts_glyco_detected <- 1L; st$n_served_as_glyco <- 1L
  out <- utils::capture.output(print(st))
  expect_true(any(grepl("detected in 1 cohort, 1 reported pair$", out))); expect_false(any(grepl("1 cohorts", out)))
})

test_that("the explanation names the parent-protein rule in plain words", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d)
  out <- paste(utils::capture.output(print(ptm_explain("TFRC_N727", "COIL_T303"))), collapse = " ")
  expect_false(grepl("protein shadow", out)); expect_true(grepl("parent-protein co-variation", out))
})

# The reference tables print as short summaries too: no "[1]" or "$name" lines of a raw R object, every line within
# 80 columns, and the object underneath unchanged (a data.frame or list with every field).
test_that("the cohorts, sites, columns, glycans, filters, citation, reference and tools print as summaries", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d)
  old <- options(width = 80); on.exit(options(old))
  shown <- list(ptm_cohorts(), ptm_sites("EGFR"), ptm_columns(), ptm_glycans(), ptm_default_filters(), ptm_cite(), ptm_reference(), ptm_tools())
  for (x in shown) {
    out <- utils::capture.output(print(x))
    expect_gt(length(out), 1L)
    expect_false(any(grepl("^\\s*(\\[[0-9]+\\]|\\$[a-z_]+)", out)), info = paste(class(x)[1], collapse = " "))
    expect_true(all(nchar(out) <= 80L), info = paste(class(x)[1], out[nchar(out) > 80L], collapse = "\n"))
  }
  expect_s3_class(ptm_cohorts(), "data.frame"); expect_true(is.list(ptm_default_filters())); expect_equal(ptm_default_filters()$min_k, 2L)
  expect_true(is.list(ptm_cite())); expect_match(ptm_cite()$text, "released 20[0-9]{2}-[0-9]{2}-[0-9]{2}")
  expect_match(paste(utils::capture.output(print(ptm_release())), collapse = "\n"), "released 20[0-9]{2}-[0-9]{2}-[0-9]{2}")
  # a cut-down table prints without row numbers
  expect_false(any(grepl("^\\s*1\\s", utils::capture.output(print(ptm_cohorts()[, c("cohort", "k_cov")])))))
})

test_that("summaries read as a header, a table and aligned next steps", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d)
  old <- options(width = 80); on.exit(options(old))
  out <- utils::capture.output(print(ptm_network("EGFR"), n = 1))
  expect_match(out[1], "^Network of EGFR: [0-9,]+ partner links? from [0-9,]+ site pairs? \\(PTMlink ")
  expect_false(any(grepl("Filtered out|edges from|Draw it|\\|", out)))
  expect_match(out[grepl("^  partner", out)], "EGFR site +partner site")
  expect_false(any(grepl("_[NSTY][0-9]+ x ", out)))                       # sites in their own columns, not "A x B"
  expect_true(any(grepl("^Showing 1 of [0-9,]+\\.$", out)))
  for (h in c("more rows: +print\\(x, n = 50\\)", "all columns: +View\\(x\\)", "plot: +ptm_plot_network\\(x\\)"))
    expect_true(any(grepl(paste0("^  ", h, "$"), out)), info = h)
  g <- ptm_gene("EGFR", "glyco"); out <- utils::capture.output(print(g, n = 1))
  expect_match(out[1], "^EGFR glycosites: [0-9,]+ pairs?, sorted by evidence \\(PTMlink ")
  expect_true(any(grepl("^  one pair: +ptm_pair\\(\"EGFR_N[0-9]+\", \"[A-Z0-9]+_[STY][0-9]+\"\\)$", out)))
  expect_true(any(grepl("^  filters: +\\?ptm_default_filters$", out)))
  g$k_cohorts[1] <- 1L; attr(g, "withheld") <- c(cis = 2, k = 3)
  expect_true(any(grepl("^Hidden by filters: 2 same-protein, 3 single-cohort\\.$", utils::capture.output(print(g)))))
  out <- utils::capture.output(print(ptm_pair("TFRC_N727", "COIL_T303")))
  expect_true(any(grepl("^  cautions: +ptm_explain\\(\"TFRC_N727\", \"COIL_T303\"\\)$", out)))
  expect_true(any(grepl("^  glycoforms: +x\\$glycoforms$", out)))
  for (k in c("pair id", "pooled rho", "evidence", "cohorts", "perturbation")) expect_true(any(grepl(paste0("^  ", k, "  "), out)), info = k)
  out <- utils::capture.output(print(ptm_sites("EGFR"), n = 1))
  expect_match(out[1], "^EGFR \\(P00533\\): [0-9]+ glycosites?, [0-9]+ phosphosites? reported \\(PTMlink ")
  expect_true(any(grepl("^Showing 1 of [0-9]+\\.$", out)))
})
