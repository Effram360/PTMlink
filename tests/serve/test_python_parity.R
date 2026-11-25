#!/usr/bin/env Rscript
# Checks that the Python reader gives the same rows, counts and verdicts as the R package on its fixture.
# tests/serve/test_python_parity.R -- the Python reader (python/ptmlink.py) against the R package, on the
# fixture the package writes. Both read the same files and must give the same rows in the same order, the
# same withheld counts, the same edges, the same rule verdicts and the same site positions. Skips with a
# message, and exit 77 so a test runner records a skip rather than a pass, where python3 or pkgload
# is absent or the served-tier arm lacks arrow/pyarrow. Pure, synthetic, read-only on the repo: the fixture
# lives in tempfile(). Run from the repo root.
ok <- TRUE; skipped <- FALSE
chk <- function(cond, msg) { if (!isTRUE(cond)) { cat("FAIL:", msg, "\n"); ok <<- FALSE } else cat("  ok:", msg, "\n") }
if (!nzchar(Sys.which("python3"))) { cat("SKIP: no python3 on PATH\n"); quit(status = 77L) }
if (!requireNamespace("pkgload", quietly = TRUE)) { cat("SKIP: pkgload not installed\n"); quit(status = 77L) }
suppressMessages(pkgload::load_all("ptmlink", quiet = TRUE))
have_arrow <- requireNamespace("arrow", quietly = TRUE)
have_pyarrow <- identical(suppressWarnings(system2("python3", c("-c", shQuote("import pyarrow")), stdout = FALSE, stderr = FALSE)), 0L)
d <- ptmlink:::.ptm_write_fixture(file.path(tempfile("parity_"), "fixture"), parquet = have_arrow)
rd <- file.path(d, "release")
ptm_set_data_dir(d); ptm_set_release_dir(rd, verify = FALSE)
py <- function(...) {
  out <- system2("python3", c("python/ptmlink.py", "--json", "--data-dir", shQuote(d), "--release-dir", shQuote(rd), "--no-verify", vapply(c(...), shQuote, character(1))), stdout = TRUE, stderr = FALSE)   # stderr carries the reader's messages, as R's message()
  st <- attr(out, "status"); if (!is.null(st) && st != 0) stop("python failed: ", paste(out, collapse = "\n"))
  jsonlite::fromJSON(paste(out, collapse = "\n"), simplifyVector = TRUE)
}
same_num <- function(a, b, tol = 1e-9) { a <- as.numeric(a); b <- as.numeric(b); length(a) == length(b) && all((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol)) }
same_chr <- function(a, b) { a <- as.character(a); b <- as.character(b); length(a) == length(b) && all((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)) }

# 1) status: the same resolution and status words
rs <- ptm_status(c("EGFR", "HER2", "P00533", "NOTAGENE")); ps <- py("status", "EGFR", "HER2", "P00533", "NOTAGENE")
chk(same_chr(rs$gene, ps$gene) && same_chr(rs$resolved_via, ps$resolved_via) && same_chr(rs$status_glyco, ps$status_glyco), "status: gene, resolution and status words agree")

# 2) gene rows: the same rows in the same order, the same numbers, the same identifiers, the same withheld counts
rg <- ptm_gene("EGFR", "glyco"); pg <- py("gene", "EGFR", "--side", "glyco")
chk(same_chr(rg$phospho_site, pg$rows$phospho_site) && same_num(rg$evidence, pg$rows$evidence), "gene EGFR glyco: rows and order agree")
chk(same_num(rg$n_total, pg$rows$n_total) && same_chr(rg$pair_id, pg$rows$pair_id), "gene EGFR glyco: n_total and pair_id agree")
chk(all(unlist(pg$withheld[names(attr(rg, "withheld"))]) == unname(attr(rg, "withheld"))), "gene EGFR glyco: withheld counts agree")
f <- modifyList(ptm_default_filters(), list(include_cis = TRUE, hide_self_coupled = TRUE, max_i2 = 0.5))
rg2 <- ptm_gene("EGFR", "glyco", f); pg2 <- py("gene", "EGFR", "--side", "glyco", "--include-cis", "--hide-self-coupled", "--max-i2", "0.5")
chk(same_chr(rg2$phospho_site, pg2$rows$phospho_site) && all(unlist(pg2$withheld[names(attr(rg2, "withheld"))]) == unname(attr(rg2, "withheld"))), "gene EGFR with three filters: rows and withheld agree")
rt <- ptm_gene("TFRC", "glyco", modifyList(ptm_default_filters(), list(include_cis = TRUE, hide_composition = TRUE))); pt <- py("gene", "TFRC", "--include-cis", "--hide-composition")
chk(same_chr(rt$phospho_site, pt$rows$phospho_site) && pt$withheld$composition == attr(rt, "withheld")[["composition"]], "the composition filter reads the used glycoform's label the same way")

# 3) a site
rsite <- ptm_site("EGFR_N603", modifyList(ptm_default_filters(), list(include_cis = TRUE))); psite <- py("site", "EGFR_N603", "--include-cis")
chk(same_chr(rsite$phospho_site, psite$rows$phospho_site) && same_num(rsite$evidence, psite$rows$evidence), "site EGFR_N603: rows agree")

# 4) a pair: the cohorts with their intervals, the representatives, the glycoform layer, the identifier
rp <- ptm_pair("EGFR_N603", "EGFR_T693"); pp <- py("pair", "EGFR_N603", "EGFR_T693")
chk(same_chr(rp$cohorts$cohort, pp$cohorts$cohort) && same_num(rp$cohorts$lo, pp$cohorts$lo) && same_num(rp$cohorts$hi, pp$cohorts$hi) && same_num(rp$cohorts$k, pp$cohorts$k) && any(rp$cohorts$k != 1), "pair: cohorts and intervals at n - 3 - k_cov agree")
chk(rp$n_total == pp$n_total && identical(rp$row$pair_id, pp$row$pair_id) && same_chr(rp$rep_glycans$rep_glycan, pp$rep_glycans$rep_glycan), "pair: n_total, pair_id and representatives agree")
chk(same_chr(rp$glycoforms$glycan, pp$glycoforms$glycan) && same_num(rp$glycoforms$q, pp$glycoforms$q) && identical(rp$glycoforms$anchored, pp$glycoforms$anchored), "pair: the glycoform layer agrees")
chk(is.null(ptm_pair("EGFR_N603", "NOPE_S1")$row) && is.null(py("pair", "EGFR_N603", "NOPE_S1")$row), "pair: an unserved pair is NULL in both")

# 5) the network: every column of every edge
rn <- ptm_network("EGFR", top_n = 2L); pn <- py("network", "EGFR", "--top-n", "2")
chk(identical(names(rn), names(pn$edges)) || setequal(names(rn), names(pn$edges)), "network: the 33 columns")
for (cn in names(rn)) {
  a <- rn[[cn]]; b <- pn$edges[[cn]]
  same <- if (is.numeric(a)) same_num(a, b) else if (is.logical(a)) identical(as.logical(a), as.logical(b)) else same_chr(a, b)
  if (!same) chk(FALSE, paste("network column", cn, "differs"))
}
chk(same_chr(rn$partner, pn$edges$partner) && identical(rn$drawn, pn$edges$drawn), "network: order, rank and drawn agree")
chk(attr(rn, "n_cis_edges") == pn$n_cis_edges && all(unlist(pn$withheld[names(attr(rn, "withheld"))]) == unname(attr(rn, "withheld"))), "network: cis count and withheld agree")

# 6) the reading rules: the same verdict and the same sentence on every rule
re_ <- ptm_explain("TFRC_N727", "COIL_T303"); pe <- py("explain", "TFRC_N727", "COIL_T303")
chk(same_chr(re_$rule, pe$rule) && identical(as.logical(re_$applies), as.logical(pe$applies)), "explain: rules and verdicts agree")
chk(same_chr(re_$statement, pe$statement), "explain: every statement is the same text")
chk(same_chr(re_$pair_id, pe$pair_id) && same_chr(re_$threshold, pe$threshold), "explain: pair_id and thresholds agree")
re2 <- ptm_explain("EGFR_N352", "SLC20A2_S256"); pe2 <- py("explain", "EGFR_N352", "SLC20A2_S256")
chk(identical(as.logical(re2$applies), as.logical(pe2$applies)) && same_chr(re2$statement, pe2$statement), "explain: a pair with missing values and a sequon flag agrees")

# 7) sites with their modeled position
rss <- ptm_sites(c("EGFR", "TFRC")); pss <- py("sites", "EGFR", "TFRC")
chk(same_chr(rss$site, pss$site) && same_chr(rss$row_id, pss$row_id) && same_num(rss$position_canonical, pss$position_canonical), "sites: keys, identifiers and positions agree")
chk(same_num(rss$plddt, pss$plddt) && same_num(rss$x, pss$x) && identical(rss$residue_matches_model, pss$residue_matches_model) && same_num(rss$frame_offset, pss$frame_offset), "sites: pLDDT, coordinates, residue check and frame offset agree")
chk(same_chr(rss$structure_status, pss$structure_status), "sites: the structure status agrees")
# 7b) the structure model packed in bundles, as the live website packs it: both readers find it through structure/bundles.json
db_ <- file.path(dirname(d), "fixture_bundled"); dir.create(db_); invisible(file.copy(list.files(d, full.names = TRUE), db_, recursive = TRUE))
raw <- paste(readLines(gzfile(file.path(db_, "structure", "EGFR.json.gz")), warn = FALSE), collapse = "\n"); dir.create(file.path(db_, "structure", "bundles"))
con <- gzfile(file.path(db_, "structure", "bundles", "0000.json.gz"), "wb"); writeLines(paste0('{"EGFR":', raw, "}"), con); close(con)
writeLines('{"format":1,"genes":{"EGFR":"structure/bundles/0000.json.gz"}}', file.path(db_, "structure", "bundles.json")); unlink(file.path(db_, "structure", "EGFR.json.gz"))
ptm_set_data_dir(db_); rsb <- suppressMessages(ptm_sites(c("EGFR", "TFRC"))); ptm_set_data_dir(d)
psb <- jsonlite::fromJSON(paste(system2("python3", c("python/ptmlink.py", "--json", "--data-dir", shQuote(db_), "sites", "EGFR", "TFRC"), stdout = TRUE, stderr = FALSE), collapse = "\n"))
chk(same_num(rsb$plddt, rss$plddt) && same_num(psb$plddt, rss$plddt) && same_chr(rsb$structure_status, psb$structure_status), "sites: a bundled structure model reads the same in both readers")

# 7c) a pair named by its pair_id, or phosphosite first; one wording for an absent pair
pid <- rp$row$pair_id
chk(identical(ptm_pair(pid)$row$pair_id, pid) && identical(py("pair", pid)$row$pair_id, pid), "pair: a pair_id names the pair in both readers")
sw_r <- suppressMessages(ptm_pair("EGFR_T693", "EGFR_N603")); sw_p <- py("pair", "EGFR_T693", "EGFR_N603")
chk(identical(sw_r$row$pair_id, pid) && identical(sw_p$row$pair_id, pid) && grepl("glycosite comes first", sw_p$input_note), "pair: phosphosite first is read in the right order in both readers")
chk(identical(ptm_pair("NOTAGENE_N1", "COIL_T303")$note, py("pair", "NOTAGENE_N1", "COIL_T303")$note) && identical(ptm_pair("EGFR_N603", "NOPE_S1")$note, py("pair", "EGFR_N603", "NOPE_S1")$note) &&
    identical(ptm_pair("EGFR_N603", "NOPE_S1")$note, "not among the reported pairs of this release; the complete scored pairs are released with the article"), "pair: one wording for an absent pair")
chk(identical(ptm_explain(pid)$statement, py("explain", pid)$statement), "explain: a pair_id gives the same rules in both readers")
# 7d) a user error on the command line is one line on stderr and a non-zero exit, never a traceback
cli <- function(...) { o <- suppressWarnings(system2("python3", c("python/ptmlink.py", "--data-dir", shQuote(d), ...), stdout = TRUE, stderr = TRUE)); list(out = o, status = attr(o, "status")) }
e1 <- cli("gene", "EGFR_N603"); e2 <- cli("gene", "EGFR", "--side", "phosphosite"); e3 <- cli("gene", "NOTAGENE"); e4 <- cli("pair", "EGFR_N603", "EGFR_N196")
chk(all(vapply(list(e1, e2, e3, e4), function(e) identical(e$status, 1L) && length(e$out) == 1L && !any(grepl("Traceback", e$out)), logical(1))) && grepl("is a site key", e1$out),
    "cli: user errors are one line with exit status 1")

# 8) the served tier, where both readers have arrow
if (have_arrow && have_pyarrow) {
  rq <- ptm_query(list(list("q_value", "<", 1e-5)), filters = modifyList(ptm_default_filters(), list(include_cis = TRUE)))
  pq <- py("query", "--where", "q_value", "<", "1e-5", "--include-cis")
  chk(same_chr(rq$pair_id, pq$rows$pair_id) && same_num(rq$evidence, pq$rows$evidence) && attr(rq, "n_matched") == pq$n_matched, "query: rows, order and n_matched agree")
  rq2 <- ptm_query(list(list("phospho_gene", "in", c("VDAC2", "COIL"))), columns = c("glyco_site", "phospho_site", "evidence"))
  pq2 <- py("query", "--where", "phospho_gene", "in", "VDAC2,COIL", "--columns", "glyco_site", "phospho_site", "evidence")
  chk(same_chr(rq2$phospho_site, pq2$rows$phospho_site) && setequal(names(rq2), names(pq2$rows)), "query: a set condition and a column projection agree")
  # 8b) a list: members, the three row tables, the partner-site and convergent tests, the connectivity expectation (the permutation is random in each reader and not compared)
  rl <- ptm_list(c("EGFR", "VDAC2", "HER2", "NOTAGENE", "TFRC_N727"), n_perm = 0); pl <- py("list", "EGFR", "VDAC2", "HER2", "NOTAGENE", "TFRC_N727", "--n-perm", "0")
  chk(same_chr(rl$members$gene, pl$members$gene) && same_chr(rl$members$resolved_via, pl$members$resolved_via) && same_chr(rl$members$site, pl$members$site) && same_num(rl$members$n_internal, pl$members$n_internal), "list: members resolve and count the same")
  chk(same_chr(rl$internal$pair_id, pl$internal$pair_id) && same_chr(rl$neighbourhood$pair_id, pl$neighbourhood$pair_id) && identical(rl$neighbourhood$site_named_glyco, pl$neighbourhood$site_named_glyco), "list: internal and neighbourhood rows agree, with the named-site flags")
  rc <- rl$convergent[order(rl$convergent$partner, rl$convergent$partner_side), ]; pc_ <- pl$convergent[order(pl$convergent$partner, pl$convergent$partner_side), ]
  same_rel <- function(a, b) { a <- as.numeric(a); b <- as.numeric(b); length(a) == length(b) && all(abs(a - b) <= 1e-9 * pmax(abs(a), abs(b)) + 1e-12) }
  chk(same_chr(rc$partner, pc_$partner) && same_num(rc$members_reached, pc_$members_reached) && same_num(rc$tier_reach, pc_$tier_reach) && same_rel(rc$fold, pc_$fold) && same_rel(rc$p, pc_$p) && same_rel(rc$q, pc_$q) && same_chr(rc$best_pair_id, pc_$best_pair_id), "list: the convergent partners, their P, q, fold and best pair agree")
  chk(rl$connectivity$internal_pairs == pl$connectivity$internal_pairs && same_rel(rl$connectivity$expected_pairs, pl$connectivity$expected_pairs) && rl$connectivity$members_in_tier == pl$connectivity$members_in_tier, "list: the connectivity count and expectation agree")
  chk(same_num(rl$composition$share, pl$composition$share) && same_chr(rl$composition$table, pl$composition$table), "list: the composition share agrees")
  rl2 <- ptm_list("EGFR", n_perm = 0); pl2 <- py("list", "EGFR", "--n-perm", "0")
  chk(same_rel(sort(rl2$convergent$p), sort(pl2$convergent$p)) && all(abs(rl2$convergent$p - 0.5) < 1e-12), "list: a one-member hypergeometric P of 0.5 in both readers")
  # the partner sites at protein level, with the members reached and the best pair
  ps <- rl$partner_sites[order(rl$partner_sites$site), ]; pps <- pl$partner_sites[order(pl$partner_sites$site), ]
  chk(same_chr(ps$site, pps$site) && same_chr(ps$gene, pps$gene) && same_chr(ps$side, pps$side) && same_num(ps$members_reached, pps$members_reached) && same_num(ps$tier_reach, pps$tier_reach) && same_num(ps$tier_n, pps$tier_n) &&
      same_rel(ps$p, pps$p) && same_rel(ps$q, pps$q) && same_rel(ps$fold, pps$fold) && same_chr(ps$members, pps$members) && same_num(ps$evidence_best, pps$evidence_best) && same_chr(ps$best_pair_id, pps$best_pair_id), "list: the partner sites, their P, q, fold, members and best pair agree")
  chk(identical(rl$level, pl$level) && rl$level == "site" && identical(vapply(rl$partner_sites, function(c) class(c)[1], ""), vapply(rl$partner_sites, function(c) class(c)[1], "")) && identical(names(rl$partner_sites), names(pl$partner_sites)), "list: the level and the partner-site columns agree")
  # the list tests count significant pairs only, in both readers: a release table with one more cross-protein pair that
  # is reported for its perturbation support alone (q = 0.6) gives the same universe, the same partners and the same P
  chk(rl$universe$pairs == pl$universe$pairs && identical(rl$universe$rule, pl$universe$rule), "list: the pairs the tests count are the same")
  rd2 <- file.path(dirname(d), "release_nonsig"); dir.create(rd2); invisible(file.copy(list.files(rd, full.names = TRUE), rd2))
  served <- as.data.frame(arrow::read_parquet(file.path(rd, "site_served.parquet")))
  extra <- served[served$glyco_site == "EGFR_N196" & served$phospho_site == "VDAC2_S115", ]
  extra$phospho_site <- "COIL_T303"; extra$phospho_gene <- "COIL"; extra$q_value <- 0.6; extra$evidence <- 0.3; extra$anch_any <- TRUE
  extra$phospho_row_id <- "PTMP:P38432:T303"; extra$pair_id <- paste0(extra$glyco_row_id, "|", extra$phospho_row_id)
  arrow::write_parquet(rbind(served, extra), file.path(rd2, "site_served.parquet"))
  ptm_set_release_dir(rd2, verify = FALSE)
  rs <- ptm_list("TFRC", n_perm = 0)
  psig <- jsonlite::fromJSON(paste(system2("python3", c("python/ptmlink.py", "--json", "--data-dir", shQuote(d), "--release-dir", shQuote(rd2), "--no-verify", "list", "TFRC", "--n-perm", "0"), stdout = TRUE, stderr = FALSE), collapse = "\n"), simplifyVector = TRUE)
  cs <- rs$partner_sites[rs$partner_sites$site == "COIL_T303", ]; cp <- psig$partner_sites[psig$partner_sites$site == "COIL_T303", ]
  chk(rs$universe$pairs == 4L && psig$universe$pairs == 4L && cs$tier_reach == 1L && cp$tier_reach == 1L && same_rel(cs$p, cp$p) && same_rel(rl$connectivity$expected_pairs, pl$connectivity$expected_pairs),
      "list: a pair reported only for its perturbation support is counted by neither reader")
  ptm_set_release_dir(rd, verify = FALSE)
  for (input in list(c("EGFR", "TFRC_N999"), c("EGFR", "EGFR_N196"), c("NOTAGENE", "TFRC_N727"))) {
    a <- ptm_list(input, n_perm = 0); b <- do.call(py, as.list(c("list", input, "--n-perm", "0")))
    chk(a$level == "site" && b$level == "site" && setequal(a$sites, b$sites) &&
        same_chr(a$neighbourhood$pair_id, b$neighbourhood$pair_id) &&
        same_num(a$members$n_neighbourhood, b$members$n_neighbourhood), paste("mixed literal scope and dedup:", paste(input, collapse = ",")))
  }
  # a site-level list: the level, the member counts, the site test against the tier's sites, the internal pair and the connectivity
  rs <- ptm_list(c("EGFR_N603", "EGFR_N196", "TFRC_N727", "VDAC2_S115"), n_perm = 0); pss <- py("list", "EGFR_N603", "EGFR_N196", "TFRC_N727", "VDAC2_S115", "--n-perm", "0")
  chk(rs$level == "site" && pss$level == "site" && same_chr(rs$members$site, pss$members$site) && same_num(rs$members$n_internal, pss$members$n_internal) && same_num(rs$members$n_neighbourhood, pss$members$n_neighbourhood), "site list: the level and the member counts agree")
  chk(same_chr(rs$partner_sites$site, pss$partner_sites$site) && same_rel(rs$partner_sites$p, pss$partner_sites$p) && same_num(rs$partner_sites$members_testable, pss$partner_sites$members_testable) && same_num(rs$partner_sites$tier_n, pss$partner_sites$tier_n) &&
      same_num(rs$partner_sites$tier_reach, pss$partner_sites$tier_reach) && same_chr(rs$partner_sites$members, pss$partner_sites$members) && identical(names(rs$partner_sites), names(pss$partner_sites)), "site list: the partner sites, P, universe and columns agree")
  chk(same_chr(rs$convergent$partner, pss$convergent$partner) && same_rel(rs$convergent$p, pss$convergent$p) && same_num(rs$convergent$tier_reach, pss$convergent$tier_reach), "site list: the convergent partners agree")
  chk(same_chr(rs$internal$pair_id, pss$internal$pair_id) && rs$connectivity$internal_pairs == pss$connectivity$internal_pairs && same_rel(rs$connectivity$expected_pairs, pss$connectivity$expected_pairs) &&
      rs$connectivity$members_in_tier == pss$connectivity$members_in_tier, "site list: the internal pair and the site-level connectivity agree")
  chk(same_chr(ptm_list("EGFR TFRC", n_perm = 0)$members$query, py("list", "EGFR TFRC", "--n-perm", "0")$members$query), "list: a pasted line splits into the same entries in both readers")
} else { skipped <- TRUE; cat("  SKIP: query parity needs arrow in R and pyarrow in Python (arrow:", have_arrow, "pyarrow:", have_pyarrow, ")\n") }

# 8c) the tool descriptions: the same names, arguments and required lists from the same reference
py_raw <- function(...) { out <- system2("python3", c("python/ptmlink.py", "--json", "--data-dir", shQuote(d), "--release-dir", shQuote(rd), "--no-verify", vapply(c(...), shQuote, character(1))), stdout = TRUE, stderr = FALSE); jsonlite::fromJSON(paste(out, collapse = "\n"), simplifyVector = FALSE) }
rt <- jsonlite::fromJSON(ptm_tools(as = "json"), simplifyVector = FALSE); pt <- py_raw("tools")
chk(identical(vapply(rt$tools, function(t) t$name, ""), vapply(pt$tools, function(t) t$name, "")), "tools: the same tool names in the same order")
chk(identical(rt$instructions, pt$instructions), "tools: the same instructions")
same_tool <- vapply(seq_along(rt$tools), function(i) { a <- rt$tools[[i]]$inputSchema; b <- pt$tools[[i]]$inputSchema
  identical(sort(names(a$properties)), sort(names(b$properties))) && identical(unlist(a$required), unlist(b$required)) && identical(rt$tools[[i]]$description, pt$tools[[i]]$description) &&
    all(vapply(names(a$properties), function(p) identical(a$properties[[p]]$type, b$properties[[p]]$type) && identical(a$properties[[p]]$anyOf, b$properties[[p]]$anyOf) && identical(a$properties[[p]]$description, b$properties[[p]]$description), logical(1))) }, logical(1))
chk(all(same_tool), paste("tools: every schema agrees (", paste(vapply(rt$tools, function(t) t$name, "")[!same_tool], collapse = ", "), ")"))
# 8d) the MCP server answers the handshake, lists the tools and runs one, capped at the row limit
mcp <- function(lines) { out <- system2("python3", c("python/ptmlink_mcp.py", "--data-dir", shQuote(d), "--release-dir", shQuote(rd), "--no-verify"), input = lines, stdout = TRUE, stderr = FALSE); lapply(out, jsonlite::fromJSON, simplifyVector = FALSE) }
r <- mcp(c('{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18","capabilities":{},"clientInfo":{"name":"test","version":"0"}}}',
           '{"jsonrpc":"2.0","method":"notifications/initialized"}',
           '{"jsonrpc":"2.0","id":2,"method":"tools/list"}',
           '{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"ptm_gene","arguments":{"gene":["EGFR"],"side":"glyco","filters":{"include_cis":true},"limit":2}}}',
           '{"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"ptm_status","arguments":{"gene":["HER2"]}}}',
           '{"jsonrpc":"2.0","id":5,"method":"tools/call","params":{"name":"ptm_nope","arguments":{}}}',
           '{"jsonrpc":"2.0","id":6,"method":"nothing/here"}'))
chk(length(r) == 6 && r[[1]]$id == 1 && !is.null(r[[1]]$result$serverInfo) && grepl("symmetric", r[[1]]$result$instructions), "mcp: initialize answers with server info and the instructions; the notification is not answered")
# the instructions name the current method and take the release and its cohort count from the data, never a fixed ten
ins <- r[[1]]$result$instructions
chk(grepl("student_normal_bh", ins) && !grepl("student_normal_by|q_bh|ten CPTAC", ins) &&
    grepl(sprintf("This data source: release fixture; %d tumor cohorts pooled (BRCA, OV, PDAC)", nrow(ptm_cohorts())), ins, fixed = TRUE),
    "mcp: the instructions say student_normal_bh, BH only, and count the cohorts from cohorts.json")
# a dictionary with a normal-tissue record (role "normal"): counted apart, never as a pooled cohort, in both readers
dn <- file.path(dirname(d), "fixture_normal"); dir.create(dn); invisible(file.copy(list.files(d, full.names = TRUE), dn, recursive = TRUE))
jsonlite::write_json(data.frame(cohort = c("BRCA", "OV", "PDAC", "OV"), k_cov = c(1L, 2L, 1L, 0L), role = c("tumour", "tumour", "tumour", "normal")), file.path(dn, "cohorts.json"))
rl <- system2("python3", c("-c", shQuote(sprintf("import sys; sys.path.insert(0, 'python'); from ptmlink import PTMlink; from ptmlink_mcp import release_line; print(release_line(PTMlink(data_dir=%s)))", deparse(dn)))), stdout = TRUE)
chk(identical(rl, " This data source: release fixture; 3 tumor cohorts pooled (BRCA, OV, PDAC); 1 normal adjacent tissue (NAT) cohorts shown beside them, view-only (OV)."),
    "mcp: a NAT cohort is named apart and never counted as pooled")
pyc <- jsonlite::fromJSON(paste(system2("python3", c("-c", shQuote(sprintf("import sys, json; sys.path.insert(0, 'python'); from ptmlink import PTMlink; print(json.dumps(PTMlink(data_dir=%s).cohorts()))", deparse(dn)))), stdout = TRUE), collapse = "\n"))
ptm_set_data_dir(dn); rc <- ptm_cohorts(); ptm_set_data_dir(d)
chk(identical(pyc$role, rc$role) && sum(rc$role == "tumour") == 3L, "cohorts: the same roles in both readers")
chk(length(r[[2]]$result$tools) == length(rt$tools) && r[[2]]$result$tools[[1]]$name == rt$tools[[1]]$name, "mcp: tools/list returns the same tools")
# a tool result is a readable summary followed by the JSON; rows keep their key columns
last <- function(x) x$result$content[[length(x$result$content)]]$text
g <- jsonlite::fromJSON(last(r[[3]]), simplifyVector = TRUE)
chk(isFALSE(r[[3]]$result$isError) && nrow(g$rows) == 2 && g$n_rows_total == 4 && g$rows$evidence[1] == 17.1, "mcp: tools/call runs ptm_gene with the view and the row cap, reporting the full count")
chk(length(r[[3]]$result$content) == 2 && grepl("sorted by evidence", r[[3]]$result$content[[1]]$text) && !"cohorts_rho" %in% names(g$rows) &&
    all(c("pair_id", "rho_lo", "rho_hi", "evidence", "q_value") %in% names(g$rows)), "mcp: ptm_gene answers with a summary, then rows of key columns with the pooled interval")
s <- jsonlite::fromJSON(last(r[[4]]), simplifyVector = TRUE)
chk(s$gene == "ERBB2" && s$resolved_via == "alias of ERBB2", "mcp: tools/call runs ptm_status")
chk(!is.null(r[[5]]$error) && r[[5]]$error$code == -32602 && grepl("no tool named", r[[5]]$error$message) && !is.null(r[[6]]$error) && r[[6]]$error$code == -32601, "mcp: an unknown tool is an invalid-params error (-32602), as the MCP specification gives it; an unknown method is -32601")
# every tool answers on the fixture without an error: the server end to end
ex <- list(ptm_release = list(), ptm_cohorts = list(), ptm_glycans = list(), ptm_status = list(gene = list("EGFR")),
           ptm_gene = list(gene = list("EGFR"), side = "glyco"), ptm_site = list(site = list(ptm_gene("EGFR", "glyco")$glyco_site[1])),
           ptm_pair = list(glyco_site = ptm_gene("EGFR", "glyco")$glyco_site[1], phospho_site = ptm_gene("EGFR", "glyco")$phospho_site[1]),
           ptm_network = list(gene = "EGFR"), ptm_explain = list(glyco_site = ptm_gene("EGFR", "glyco")$glyco_site[1], phospho_site = ptm_gene("EGFR", "glyco")$phospho_site[1]),
           ptm_sites = list(x = list("EGFR")), ptm_columns = list(), ptm_reference = list(), ptm_cite = list(),
           ptm_query = list(genes = list("EGFR"), limit = 3), ptm_list = list(x = list("EGFR", "TFRC")))
calls <- vapply(seq_along(ex), function(i) as.character(jsonlite::toJSON(list(jsonrpc = "2.0", id = 100 + i, method = "tools/call",
  params = list(name = names(ex)[i], arguments = if (length(ex[[i]])) ex[[i]] else structure(list(), names = character(0)))), auto_unbox = TRUE)), "")
names_listed <- vapply(rt$tools, function(t) t$name, "")
ra <- mcp(c('{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18","capabilities":{},"clientInfo":{"name":"test","version":"0"}}}', calls))
errs <- vapply(ra[-1], function(x) isTRUE(x$result$isError) || !is.null(x$error), logical(1))
chk(setequal(names(ex), names_listed) && !any(errs), paste("mcp: every tool answers on the fixture (", paste(names(ex)[errs], collapse = ", "), ")"))

# 8e) the MCP server keeps the columns of a capped table that is not pair rows, rejects what the schema does not
# allow, reads a pair_id, and says what a limit of 0 did
rm_ <- mcp(c('{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18","capabilities":{},"clientInfo":{"name":"test","version":"0"}}}',
             '{"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"ptm_sites","arguments":{"x":"EGFR","limit":1}}}',
             '{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"ptm_gene","arguments":{"gene":"EGFR","filters":{"bogus":1}}}}',
             sprintf('{"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"ptm_pair","arguments":{"glyco_site":"%s"}}}', pid),
             '{"jsonrpc":"2.0","id":5,"method":"tools/call","params":{"name":"ptm_gene","arguments":{"gene":"EGFR","limit":0}}}',
             '{"jsonrpc":"2.0","id":6,"method":"tools/call","params":{"name":"ptm_gene","arguments":{"gene":"EGFR","nonsense":true}}}',
             '{"jsonrpc":"2.0","id":7,"method":"tools/call","params":{"name":"ptm_pair","arguments":{"glyco_site":"EGFR_T693","phospho_site":"EGFR_N603"}}}'))
st_ <- jsonlite::fromJSON(last(rm_[[2]]), simplifyVector = TRUE)
chk(isFALSE(rm_[[2]]$result$isError) && st_$n_rows_total == 4 && nrow(st_$rows) == 1 && all(c("site", "row_id", "plddt", "structure_status") %in% names(st_$rows)), "mcp: a capped site table keeps its columns")
chk(isTRUE(rm_[[3]]$result$isError) && grepl("unknown filter \"bogus\"", rm_[[3]]$result$content[[1]]$text, fixed = TRUE) && isTRUE(rm_[[6]]$result$isError) && grepl("has no argument \"nonsense\"", rm_[[6]]$result$content[[1]]$text, fixed = TRUE), "mcp: unknown filters and arguments are errors")
chk(isFALSE(rm_[[4]]$result$isError) && identical(jsonlite::fromJSON(last(rm_[[4]]))$row$pair_id, pid), "mcp: ptm_pair reads a pair_id")
chk(isFALSE(rm_[[5]]$result$isError) && grepl("the limit was 0", rm_[[5]]$result$content[[1]]$text), "mcp: a limit of 0 is said, not shown as no pair")
chk(isFALSE(rm_[[7]]$result$isError) && identical(jsonlite::fromJSON(last(rm_[[7]]))$row$pair_id, pid) && any(grepl("glycosite comes first", vapply(rm_[[7]]$result$content, function(c) c$text, ""))), "mcp: a phosphosite-first pair is read in order, with the note")

# 9) the citation and the reference come from the same files
pc <- py("cite")
chk(identical(ptm_cite()$title, pc$title) && identical(ptm_cite()$text, pc$text) && is.null(pc$authors) && is.null(ptm_cite()$authors), "cite: the same title and citation line, with no author list")
# the shipped citation file (no release directory set) names the preprint, and both readers print the same line from it
local({
  rd <- ptmlink:::.ptm_env$release_dir; assign("release_dir", NULL, envir = ptmlink:::.ptm_env); on.exit(assign("release_dir", rd, envir = ptmlink:::.ptm_env))
  old <- Sys.getenv("PTMLINK_RELEASE_DIR", NA); Sys.unsetenv("PTMLINK_RELEASE_DIR"); on.exit(if (!is.na(old)) Sys.setenv(PTMLINK_RELEASE_DIR = old), add = TRUE)
  rc <- ptm_cite()
  pj <- jsonlite::fromJSON(paste(system2("python3", c("python/ptmlink.py", "--json", "--data-dir", shQuote(d), "cite"), stdout = TRUE, stderr = FALSE), collapse = "\n"), simplifyVector = FALSE)
  chk(identical(rc$text, pj$text) && identical(rc$paper$doi, pj$paper$doi) && identical(rc$paper$journal, pj$paper$journal) && grepl("bioRxiv, 2026. doi:10.", rc$text, fixed = TRUE),
      "cite: the shipped citation names the preprint, the same line in both readers")
})
# the package's citation file is the repository's CITATION.cff, byte for byte (the export writes one from the other)
cff <- if (file.exists("public/CITATION.cff")) "public/CITATION.cff" else "CITATION.cff"
chk(identical(unname(tools::md5sum(cff)), unname(tools::md5sum("ptmlink/inst/extdata/ptmlink_citation.cff"))), "cite: the package's citation file is the repository's CITATION.cff")
chk(length(ptm_reference()$rules) == nrow(py("reference")$rules), "reference: the same rule count")

if (!ok) quit(status = 1)
# Failure dominates; a skipped served-tier arm is unverified, not passed.
if (skipped) { cat("test_ptmlink_python_parity: no failures, but the served-tier arm was SKIPPED\n"); quit(status = 77L) }
cat("test_ptmlink_python_parity: all checks passed\n")
