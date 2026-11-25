# The package must install and check with NO release present: these tests build a tiny
# fixture data directory (two genes, six rows) in tempdir() and point the package at it.
make_fixture <- function() ptmlink:::.ptm_write_fixture(parquet = FALSE)
test_that("the package works from a local data directory and needs no release to check", {
  d <- make_fixture(); ptm_set_data_dir(d)
  expect_equal(ptm_release()$release, "fixture")
  st <- ptm_status("egfr"); expect_equal(st$gene, "EGFR"); expect_equal(st$status_glyco, "usable")
  expect_equal(ptm_status("ERBB1")$gene, "EGFR"); expect_match(ptm_status("ERBB1")$resolved_via, "alias")
  expect_equal(ptm_status("NOTAGENE")$status_glyco, "not a recognized symbol")
  expect_equal(ptm_status("CD274")$status_glyco, "protein measured, no site")
})
test_that("ptm_gene obeys the default view: cis hidden with a count, sorted on evidence, k >= 2, q < 0.05", {
  d <- make_fixture(); ptm_set_data_dir(d)
  g <- ptm_gene("EGFR", "glyco")
  expect_false(any(g$is_cis)); expect_equal(unname(attr(g, "withheld")["cis"]), 1L)
  expect_true(all(diff(g$evidence) <= 0))
  # the q filter is ON by default, as on the website; every fixture row is
  # significant, so the default view withholds none of them for that reason
  expect_true(ptm_default_filters()$significant_only)
  expect_true(all(g$q_value < 0.05)); expect_equal(unname(attr(g, "withheld")["significant"]), 0L)
  expect_equal(nrow(ptm_gene("EGFR", "glyco", modifyList(ptm_default_filters(), list(include_cis = TRUE)))), 4L)
  expect_equal(nrow(ptm_gene("EGFR", "glyco", modifyList(ptm_default_filters(), list(anchored_only = TRUE)))), 1L)
  expect_true(is.data.frame(attr(g, "sites")))
  # the label filter reads the PAIR (both sides labeled); TFRC_N727 is labeled, COIL_T303 is not, so the row stays
  tf <- ptm_gene("TFRC", "glyco", modifyList(ptm_default_filters(), list(hide_composition = TRUE)))
  expect_equal(nrow(tf), 1L); expect_equal(unname(attr(tf, "withheld")["composition"]), 0L)
  expect_null(suppressMessages(ptm_gene("CD274")))
})
test_that("an I2 ceiling withholds heterogeneous rows and NA applies none", {
  d <- make_fixture(); ptm_set_data_dir(d)
  g <- ptm_gene("EGFR", "glyco")
  h <- ptm_gene("EGFR", "glyco", modifyList(ptm_default_filters(), list(max_i2 = 0.5)))
  expect_true(all(is.na(h$I2) | h$I2 <= 0.5))
  expect_equal(unname(attr(h, "withheld")["heterogeneous"]), sum(g$I2 > 0.5, na.rm = TRUE))
  expect_equal(nrow(h) + unname(attr(h, "withheld")["heterogeneous"]), nrow(g))
  expect_true(is.na(ptm_default_filters()$max_i2))
})

test_that("the adjusted-interval filter is honoured when the columns exist and withholds everything when they do not", {
  d <- make_fixture(); ptm_set_data_dir(d)
  g <- ptm_gene("EGFR", "glyco", modifyList(ptm_default_filters(), list(adjusted_excludes_zero = TRUE)))
  expect_equal(nrow(g), 0L)   # the fixture has no adjusted block, so no row can qualify
  expect_true(unname(attr(g, "withheld")["adjusted"]) > 0)
  expect_false(ptm_default_filters()$adjusted_excludes_zero)
})

test_that("the q filter is on by default and unticking it returns the served rows whole", {
  d <- make_fixture(); ptm_set_data_dir(d)
  # the fixture carries no q >= 0.05 row, so raise one and charge the default view against it
  raw <- ptm_gene("EGFR", "glyco", modifyList(ptm_default_filters(), list(significant_only = FALSE)))
  raw$q_value[1] <- 0.2
  on <- ptmlink:::.ptm_filter(raw, ptm_default_filters())
  off <- ptmlink:::.ptm_filter(raw, modifyList(ptm_default_filters(), list(significant_only = FALSE)))
  expect_equal(nrow(on$rows), nrow(raw) - 1L); expect_equal(unname(on$withheld["significant"]), 1L)
  expect_equal(nrow(off$rows), nrow(raw)); expect_equal(unname(off$withheld["significant"]), 0L)
})

test_that("ptm_site and ptm_pair read the packed per-cohort breakdown", {
  d <- make_fixture(); ptm_set_data_dir(d)
  s <- ptm_site("EGFR_N603"); expect_true(all(s$glyco_site == "EGFR_N603")); expect_equal(nrow(s), 1L)   # cis hidden by default
  p <- ptm_pair("EGFR_N603", "EGFR_T693")
  expect_equal(nrow(p$cohorts), 2L); expect_equal(p$cohorts$cohort, c("BRCA", "OV")); expect_true(p$cohorts$lo[1] < 0.40 && p$cohorts$hi[1] > 0.40)
  expect_equal(nrow(p$rep_glycans), 2L)
  expect_null(ptm_pair("EGFR_N603", "NOPE_S1")$row)
  # the glycoform layer of the pair: one row per glycoform, evidence descending, the layer's own q
  g <- p$glycoforms
  expect_equal(nrow(g), 2L); expect_equal(g$glycan, c("N2H8F0S0G0", "N2H9F0S0G0")); expect_equal(g$rho[1], 0.4012); expect_equal(g$k, c(5L, 4L)); expect_equal(g$anchored, c(TRUE, FALSE))
  expect_equal(nrow(ptm_pair("EGFR_N603", "ALG3_S13")$glycoforms), 1L)
  expect_null(ptm_pair("TFRC_N727", "COIL_T303")$glycoforms)
})

test_that("the shadow-clean filter keeps only pairs whose parents do not co-vary, and n_total rides every row", {
  d <- make_fixture(); ptm_set_data_dir(d)
  g <- ptm_gene("EGFR", "glyco")
  expect_true("n_total" %in% names(g)); expect_equal(g$n_total[g$phospho_site == "ALG3_S13"], 40L)
  sc <- ptm_gene("EGFR", "glyco", modifyList(ptm_default_filters(), list(shadow_clean_only = TRUE)))
  expect_true(all(sc$shadow_clean)); expect_equal(nrow(sc), 3L)
  both <- ptm_gene("EGFR", "glyco", modifyList(ptm_default_filters(), list(include_cis = TRUE, shadow_clean_only = TRUE)))
  expect_equal(nrow(both), 3L); expect_equal(unname(attr(both, "withheld")["shadow"]), 1L)   # the cis pair is never shadow-clean
  expect_false(ptm_default_filters()$shadow_clean_only)
  expect_equal(ptm_pair("EGFR_N603", "EGFR_T693")$n_total, 75L)
})
test_that("hide_self_coupled withholds pairs at self_coupling >= 0.6 and keeps pairs with no value", {
  d <- make_fixture(); ptm_set_data_dir(d)
  expect_false(ptm_default_filters()$hide_self_coupled)
  g <- ptm_gene("TFRC", "glyco", modifyList(ptm_default_filters(), list(include_cis = TRUE)))
  expect_equal(nrow(g), 2L)
  h <- ptm_gene("TFRC", "glyco", modifyList(ptm_default_filters(), list(include_cis = TRUE, hide_self_coupled = TRUE)))
  expect_equal(nrow(h), 0L); expect_equal(unname(attr(h, "withheld")["self"]), 2L)
  e <- ptm_gene("EGFR", "glyco", modifyList(ptm_default_filters(), list(hide_self_coupled = TRUE)))
  expect_true("EGFR_N352" %in% e$glyco_site, "a pair with no self_coupling value is kept"); expect_equal(unname(attr(e, "withheld")["self"]), 0L)
  # the cutoff is the article's candidate screen's, 0.6: a pair at 0.65 is withheld, one at 0.55 kept
  rows <- g[c(1L, 1L), ]; rows$self_coupling <- c(0.65, 0.55)
  kept <- ptmlink:::.ptm_filter(rows, modifyList(ptm_default_filters(), list(include_cis = TRUE, hide_self_coupled = TRUE)))$rows
  expect_equal(kept$self_coupling, 0.55)
})

test_that("ptm_status, ptm_gene, ptm_site and ptm_pair take vectors and stack their results", {
  d <- make_fixture(); ptm_set_data_dir(d)
  st <- ptm_status(c("EGFR", "P00533", "NOTAGENE")); expect_equal(nrow(st), 3L); expect_equal(st$gene, c("EGFR", "EGFR", NA)); expect_match(st$resolved_via[2], "accession")
  g <- ptm_gene(c("EGFR", "TFRC"), "glyco"); expect_equal(names(g)[1], "query_gene"); expect_equal(sort(unique(g$query_gene)), c("EGFR", "TFRC"))
  expect_true(is.matrix(attr(g, "withheld"))); expect_equal(unname(attr(g, "withheld")["EGFR", "cis"]), 1L)
  s <- ptm_site(c("EGFR_N603", "TFRC_N727")); expect_equal(names(s)[1], "query_site"); expect_equal(nrow(s), 2L)
  p <- ptm_pair(c("EGFR_N603", "TFRC_N727", "EGFR_N603"), c("EGFR_T693", "COIL_T303", "NOPE_S1"))
  expect_equal(nrow(p$row), 2L); expect_equal(nrow(p$missing), 1L); expect_equal(unname(p$n_total), c(75L, 130L, NA))
  expect_equal(names(p$cohorts)[1:2], c("glyco_site", "phospho_site")); expect_equal(nrow(p$cohorts), 4L)
})
test_that("every returned row carries its stable identifiers; a partner side comes from the shipped map, its shard or the release", {
  d <- make_fixture(); ptm_set_data_dir(d)
  g <- ptm_gene("EGFR", "glyco")
  expect_true(all(c("glyco_row_id", "phospho_row_id", "pair_id") %in% names(g)))
  expect_true(all(g$glyco_row_id == "PTMG:P00533:N603" | grepl("^PTMG:P00533:N", g$glyco_row_id)))
  # partner phosphosites live in shards not read; the package's site map (inst/extdata/site_ids.tsv.gz) knows them
  expect_equal(g$phospho_row_id[g$phospho_site == "VDAC2_S115"], "PTMP:P45880:S115")
  expect_false(anyNA(g$pair_id))
  p <- ptm_pair("TFRC_N727", "TFRC_S34"); expect_equal(p$row$pair_id, "PTMG:P02786:N727|PTMP:P02786:S34")
  ptm_set_release_dir(file.path(d, "release"), verify = FALSE)   # the feature tables complete the map
  g2 <- ptm_gene("EGFR", "glyco"); expect_equal(g2$pair_id[g2$phospho_site == "ALG3_S13"], "PTMG:P00533:N603|PTMP:Q92685:S13")
})
test_that("ptm_pair draws each cohort's interval at n - 3 - k_cov from the dictionary, and n - 4 where it states none", {
  d <- ptmlink:::.ptm_write_fixture(file.path(tempfile("kcov_"), "fx"), parquet = FALSE); ptm_set_data_dir(d)
  co <- ptm_pair("EGFR_N603", "EGFR_T693")$cohorts
  iv <- function(r, w) round(tanh(atanh(r) + c(-1, 1) * 1.959964 / sqrt(w)), 4)   # to four decimals, the precision rho is released at
  expect_equal(co$k_cov, c(1, 2)); expect_false("k" %in% names(co))                # the covariate count, never mistaken for the glycoform table's cohort count k
  expect_equal(c(co$lo[1], co$hi[1]), iv(0.40, 45 - 4), tolerance = 1e-12)          # BRCA, one covariate
  expect_equal(c(co$lo[2], co$hi[2]), iv(0.30, 30 - 3 - 2), tolerance = 1e-12)      # OV, two covariates in the fixture
  # a dictionary without k_cov and role: every interval at n - 4, every cohort a tumor cohort
  writeLines(jsonlite::toJSON(data.frame(cohort = c("BRCA", "OV"), tumour_type = c("breast", "ovarian")), dataframe = "rows"), file.path(d, "cohorts.json"))
  ptm_set_data_dir(d)
  co11 <- ptm_pair("EGFR_N603", "EGFR_T693")$cohorts
  expect_equal(co11$k_cov, c(1, 1)); expect_equal(c(co11$lo[2], co11$hi[2]), iv(0.30, 30 - 4), tolerance = 1e-12)
  expect_equal(ptm_cohorts()$role, c("tumour", "tumour"))
  # a normal-tissue record of a code never sets that code's tumor k
  writeLines(jsonlite::toJSON(data.frame(cohort = c("BRCA", "OV", "OV"), k_cov = c(1, 2, 0), role = c("tumour", "tumour", "normal")), dataframe = "rows"), file.path(d, "cohorts.json"))
  ptm_set_data_dir(d)
  expect_equal(ptm_pair("EGFR_N603", "EGFR_T693")$cohorts$k, c(1, 2)); expect_equal(sum(ptm_cohorts()$role == "tumour"), 2L)
})

# Where the data come from, and what a failure is called. A file the source does not carry is an absence (a gene
# with no reported pair); a host that cannot be reached is an error that says so, never "not a reported pair". The
# structure model is found through structure/bundles.json as the website finds it, else structure/<GENE>.json.gz.
# Every source here is a local directory or a file:// address; nothing leaves the machine.

# the fixture with its one structure record moved into a bundle, as the live website packs them
bundled_fixture <- function() {
  d <- ptmlink:::.ptm_write_fixture(file.path(tempfile("bundled_"), "fixture"), parquet = FALSE)
  raw <- paste(readLines(gzfile(file.path(d, "structure", "EGFR.json.gz")), warn = FALSE), collapse = "\n")
  dir.create(file.path(d, "structure", "bundles"))
  con <- gzfile(file.path(d, "structure", "bundles", "0000.json.gz"), "wb"); writeLines(paste0('{"EGFR":', raw, "}"), con); close(con)
  writeLines('{"format":1,"genes":{"EGFR":"structure/bundles/0000.json.gz"}}', file.path(d, "structure", "bundles.json"))
  unlink(file.path(d, "structure", "EGFR.json.gz"))
  d
}
# a directory published as a remote source: every file listed with its checksum, as the website's SHA256SUMS.json
publish <- function(d) {
  files <- list.files(d, recursive = TRUE); files <- files[!startsWith(files, "release/") & files != "SHA256SUMS.json"]
  sums <- setNames(as.list(unname(tools::sha256sum(file.path(d, files)))), files)
  writeLines(jsonlite::toJSON(sums, auto_unbox = TRUE), file.path(d, "SHA256SUMS.json"))
  paste0("file://", normalizePath(d))
}

test_that("structure is read through structure/bundles.json, locally and from a remote source", {
  d <- bundled_fixture(); ptm_set_data_dir(d)
  s <- ptm_sites("EGFR")
  expect_equal(s$plddt[s$site == "EGFR_N603"], 90); expect_equal(s$plddt[s$site == "EGFR_T693"], 40); expect_true(all(s$structure_status == "model read"))
  expect_message(t <- ptm_sites("TFRC"), "no AlphaFold model")
  expect_true(all(is.na(t$plddt))); expect_true(all(t$structure_status == "no AlphaFold model in this data source"))
  withr::local_envvar(R_USER_CACHE_DIR = tempfile("cache_"))
  ptm_set_data_url(publish(d))
  r <- ptm_sites("EGFR"); expect_equal(r$plddt, s$plddt); expect_equal(r$x, s$x)
  # the per-gene layout still reads where no bundle index exists (a data directory written before bundles)
  d2 <- ptmlink:::.ptm_write_fixture(file.path(tempfile("plain_"), "fixture"), parquet = FALSE); ptm_set_data_dir(d2)
  expect_equal(ptm_sites("EGFR")$plddt, s$plddt)
})

test_that("a structure model that cannot be fetched is said plainly and does not fail the call", {
  withr::local_envvar(R_USER_CACHE_DIR = tempfile("cache_"))
  d <- bundled_fixture(); ptm_set_data_url(publish(d))
  unlink(file.path(d, "structure", "bundles", "0000.json.gz"))       # listed, but the host cannot serve it
  expect_message(s <- ptm_sites("EGFR"), "could not be read")
  expect_equal(nrow(s), 4L); expect_true(all(is.na(s$plddt))); expect_match(s$structure_status[1], "^could not be read: could not reach")
})

test_that("an unreachable host is an error that says so, never an absent pair or gene", {
  withr::local_envvar(R_USER_CACHE_DIR = tempfile("cache_"))
  ptm_set_data_url(paste0("file://", tempfile("nowhere_")))
  expect_error(ptm_pair("IL6ST_N131", "TMX1_S270"), class = "ptm_unreachable")
  expect_error(ptm_pair("IL6ST_N131", "TMX1_S270"), "could not reach")
  expect_error(ptm_status("EGFR"), class = "ptm_unreachable")
  expect_error(ptm_gene("EGFR"), class = "ptm_unreachable")
})

test_that("offline with a cache, the cached checksum list is used with a warning", {
  withr::local_envvar(R_USER_CACHE_DIR = tempfile("cache_"))
  d <- bundled_fixture(); url <- publish(d)
  ptm_set_data_url(url); expect_false(is.null(ptm_pair("EGFR_N603", "ALG3_S13")$row))
  unlink(d, recursive = TRUE)                                         # the host is gone
  ptm_set_data_url(url)
  expect_warning(p <- ptm_pair("EGFR_N603", "ALG3_S13"), "checksum list cached")
  expect_false(is.null(p$row)); expect_equal(p$row$k_cohorts, 5L)
  expect_error(ptm_pair("TFRC_N727", "COIL_T303"), class = "ptm_unreachable")   # never cached: unreachable, not absent
})

test_that("an absent pair has one wording, and a gene without a data file is absent, not an error", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d)
  a <- ptm_pair("NOTAGENE_N1", "TMX1_S270"); b <- ptm_pair("EGFR_N603", "NOPE_S1")
  expect_null(a$row); expect_null(b$row)
  expect_equal(a$note, "not among the reported pairs of this release; the complete scored pairs are released with the article"); expect_equal(b$note, a$note)
  expect_output(print(a), "Not a reported pair: not among the reported pairs")
})

test_that("a pair is named by its pair_id, by site identifiers, or with the phosphosite first", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d)
  ref <- ptm_pair("EGFR_N603", "ALG3_S13")
  byid <- ptm_pair("PTMG:P00533:N603|PTMP:Q92685:S13")
  expect_equal(byid$row$pair_id, "PTMG:P00533:N603|PTMP:Q92685:S13"); expect_equal(byid$row$evidence, ref$row$evidence)
  expect_equal(ptm_pair("ptmg:p00533:n603|ptmp:q92685:s13")$row$pair_id, ref$row$pair_id)          # any case
  expect_equal(ptm_pair("PTMG:P00533:N603", "ALG3_S13")$row$pair_id, ref$row$pair_id)                # a site identifier for a key
  expect_message(sw <- ptm_pair("ALG3_S13", "EGFR_N603"), "the glycosite comes first")
  expect_equal(sw$row$pair_id, ref$row$pair_id)
  expect_error(ptm_pair("EGFR_N603", "EGFR_N196"), "one glycosite")
  expect_error(ptm_pair("EGFR_N603"), "pair_id")
  expect_error(ptm_pair("PTMG:P00000:N1|PTMP:P00000:S1"), "not a site identifier")
  e <- ptm_explain("PTMG:P00533:N603|PTMP:Q92685:S13"); expect_equal(unique(e$pair_id), ref$row$pair_id)
  two <- ptm_pair(c("PTMG:P00533:N603|PTMP:Q92685:S13", "PTMG:P02786:N727|PTMP:P38432:T303"))
  expect_equal(nrow(two$row), 2L)
})

test_that("site keys work in any case, and a site key given as a gene points to ptm_site()", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d)
  expect_equal(nrow(ptm_site("egfr_n603")), nrow(ptm_site("EGFR_N603")))
  expect_error(ptm_gene("EGFR_N603"), "is a site key, not a gene; ptm_site")
  expect_error(ptm_gene("NOTAGENE"), "not a recognized symbol")
  expect_error(ptm_site("NOTAGENE_N1"), "in no reported pair")
})
test_that("a site key is read in the release's own letter case, so a gene with lowercase letters is found", {
  d <- make_fixture(); ptm_set_data_dir(d)
  # a known gene keeps the index's spelling whatever was typed; an unknown gene and an identifier come back upper-cased
  expect_equal(ptmlink:::.ptm_canon_site(c("egfr_n603", " EGFR_s1064 ", "nope_s1", "ptmg:p00533:n603", "EGFR")),
               c("EGFR_N603", "EGFR_S1064", "NOPE_S1", "PTMG:P00533:N603", "EGFR"))
  # C10orf71 and its like: put one in the cached index and type it in other cases
  idx <- ptmlink:::.ptm_index(); extra <- idx$genes[1, , drop = FALSE]; extra$gene <- "C10orf71"; idx$genes <- rbind(idx$genes, extra)
  env <- ptmlink:::.ptm_env; env$cache[["index.json.gz"]] <- idx
  expect_equal(ptmlink:::.ptm_canon_site(c("c10ORF71_s677", "C10orf71_S677", "EGFR_N603")), c("C10orf71_S677", "C10orf71_S677", "EGFR_N603"))
  expect_equal(ptmlink:::.ptm_canon_gene(c("c10orf71", "nope")), c("C10orf71", "NOPE"))
  k <- ptmlink:::.ptm_site_key("c10orf71_s677"); expect_equal(k$gene, "C10orf71"); expect_equal(k$key, "C10orf71_S677"); expect_equal(k$side, "phospho")
  ptm_set_data_dir(d)   # drop the edited index
})


test_that("a partial filter list is completed from the defaults; an unknown name or a wrong kind of value is refused", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d)
  full <- ptm_gene("EGFR", "glyco", modifyList(ptm_default_filters(), list(include_cis = TRUE)))
  part <- ptm_gene("EGFR", "glyco", list(include_cis = TRUE))
  expect_equal(part$pair_id, full$pair_id); expect_equal(attr(part, "filters")$min_k, 2L); expect_true(attr(part, "filters")$include_cis)
  expect_equal(ptm_gene("EGFR", "glyco", NULL)$pair_id, ptm_gene("EGFR", "glyco")$pair_id)
  expect_error(ptm_gene("EGFR", "glyco", list(min_kk = 1)), "unknown filter \"min_kk\"; the filters are include_cis, min_k")
  expect_error(ptm_gene("EGFR", "glyco", list(min_k = "two")), "min_k is a whole number")
  expect_error(ptm_gene("EGFR", "glyco", list(significant_only = "no")), "significant_only is TRUE or FALSE")
  expect_error(ptm_gene("EGFR", "glyco", list(max_i2 = 5)), "max_i2 is a number between 0 and 1")
  expect_error(ptm_site("EGFR_N603", list(TRUE)), "every filter needs its name")
  expect_error(ptm_network("EGFR", "glyco"), "filters is a list of settings")
  expect_true(is.na(ptmlink:::.ptm_filters(list(max_i2 = NULL))$max_i2))
})

test_that("names are checked before anything is read, and an unknown name among several is skipped with a message", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d)
  expect_error(ptm_gene(NA), "gene is NA"); expect_error(ptm_gene(123), "gene is numeric, not text"); expect_error(ptm_gene(character(0)), "gene is empty")
  expect_error(ptm_gene(""), "missing or empty entry"); expect_error(ptm_status(NULL), "gene is empty")
  expect_equal(ptm_status(" egfr ")$query, "egfr")
  expect_error(ptm_gene("NOTAGENE"), "not a recognized symbol")
  expect_message(g <- ptm_gene(c("EGFR", "NOTAGENE", "EGFR_T693"), "glyco"), "NOTAGENE skipped: not a recognized symbol")
  expect_equal(unique(g$query_gene), "EGFR"); expect_equal(nrow(g), nrow(ptm_gene("EGFR", "glyco")))
  expect_error(ptm_network(c("EGFR", "TFRC")), "takes one gene, got 2")
  # a site in no reported pair is said as such, never shown as hidden by the filters
  expect_error(ptm_site("EGFR_N9999"), "EGFR_N9999 is in no reported pair of this release; ptm_sites\\(\"EGFR\"\\) lists")
  expect_message(s <- ptm_site(c("EGFR_T693", "EGFR_N9999"), list(include_cis = TRUE)), "EGFR_N9999 is in no reported pair")
  expect_equal(unique(s$query_site), "EGFR_T693")
  expect_error(ptm_set_data_dir(file.path(tempdir(), "no_such_dir_here")), "no such directory"); expect_error(ptm_set_data_dir(tempdir()), "is not a PTMlink data directory")
  expect_error(ptm_set_data_url("not an address"), "url is one web address")
  ptm_set_data_dir(d)
})

test_that("a printed header states only what is true of the rows in hand", {
  d <- ptmlink:::.ptm_write_fixture(parquet = FALSE); ptm_set_data_dir(d)
  old <- options(width = 80); on.exit(options(old))
  g <- ptm_gene("EGFR", "glyco", list(include_cis = TRUE)); n <- nrow(g)
  expect_match(utils::capture.output(print(g))[1], paste0("^EGFR glycosites: ", n, " pairs, sorted by evidence "))
  expect_match(utils::capture.output(print(g[1:2, ]))[1], paste0("^EGFR glycosites: 2 of ", n, " pairs, sorted by evidence "))
  expect_false(grepl("sorted", utils::capture.output(print(g[order(g$evidence), ]))[1]))          # re-sorted by the user: no claim
  two <- ptm_gene(c("EGFR", "TFRC"), "glyco", list(include_cis = TRUE)); out <- utils::capture.output(print(two, n = 2))   # EGFR's rows first, though TFRC holds the strongest pair
  expect_true(any(grepl("sorted by evidence within each gene", out[1:2]))); expect_true(any(grepl("^  TFRC ", out)) && any(grepl("^  EGFR ", out)))   # the top row of each gene
  expect_true(any(grepl("^Showing the top 1 of each gene: 2 of ", out)))
  p <- ptm_pair(c("EGFR_N603", "TFRC_N727"), c("ALG3_S13", "COIL_T303")); out <- utils::capture.output(print(p))
  expect_match(out[1], "^The pairs asked: 2 pairs \\(PTMlink ")                                   # in the order asked, 5.3 before 21.6
  k1 <- ptm_gene("EGFR", "glyco", list(min_k = 6)); expect_true(any(grepl("in fewer than 6 cohorts", utils::capture.output(print(k1)))))
})
