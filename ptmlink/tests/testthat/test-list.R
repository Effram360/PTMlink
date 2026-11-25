# ptm_list on the fixture: resolution and coverage, the three row tables under the view, the tests on the tier as
# released (partner sites, convergent partners, connectivity), the composition share, the printed
# summary, at protein level and at site level. The fixture's tier holds four rows once same-protein rows are
# dropped: EGFR_N603-ALG3_S13, EGFR_N196-VDAC2_S115, EGFR_N352-SLC20A2_S256, TFRC_N727-COIL_T303.
fixture <- function() { d <- ptmlink:::.ptm_write_fixture(parquet = TRUE); ptm_set_data_dir(d); ptm_set_release_dir(file.path(d, "release")); d }

test_that("entries are parsed from elements, pasted text and a data.frame: every token is an entry", {
  e <- ptmlink:::.ptm_entries("EGFR, MET\nTFRC; CD44\nITGB1\nAKT1 MTOR\nEGFR")
  expect_equal(e$query, c("EGFR", "MET", "TFRC", "CD44", "ITGB1", "AKT1", "MTOR"))   # a repeated entry is kept once, where it first appears
  expect_equal(names(e), "query")                                                     # an entry is a token and nothing else
  expect_equal(ptmlink:::.ptm_entries(c("EGFR", "MET"))$query, c("EGFR", "MET"))
  expect_equal(ptmlink:::.ptm_entries(data.frame(id = c("EGFR", "TFRC"), lfc = c(0, -2)))$query, c("EGFR", "TFRC"))   # only the first column is read
  expect_equal(ptmlink:::.ptm_entries("EGFR up")$query, c("EGFR", "up"))              # a trailing word is an entry of its own, never stripped
  expect_equal(nrow(ptmlink:::.ptm_entries(c("", " ", NA))), 0L)
})

test_that("the tier is read once into codes, degrees and reaches, and is the same as the website's sidecar", {
  testthat::skip_if_not_installed("arrow")
  d <- fixture(); tier <- ptmlink:::.ptm_tier()
  expect_equal(tier$M, 4L); expect_equal(tier$genes_glyco, c("EGFR", "TFRC")); expect_equal(tier$genes_phospho, c("ALG3", "COIL", "SLC20A2", "VDAC2"))
  expect_equal(tier$sites_glyco, c("EGFR_N196", "EGFR_N352", "EGFR_N603", "TFRC_N727")); expect_equal(tier$site_gene_glyco, c("EGFR", "EGFR", "EGFR", "TFRC"))
  expect_equal(tier$gene_glyco["EGFR", ], data.frame(deg = 3L, reach_genes = 3L, reach_sites = 3L, row.names = "EGFR"))
  expect_equal(unname(tier$site_phospho["COIL_T303", "reach_genes"]), 1L); expect_equal(sum(tier$rows$gg == 1L), 3L)
  expect_identical(ptmlink:::.ptm_tier(), tier)   # cached
  js <- jsonlite::fromJSON(gzfile(file.path(d, "tier_significant.json.gz")))
  expect_true(js$significant_only); expect_match(js$rule, "q_value < 0.05", fixed = TRUE)
  expect_equal(js$M, 4L); expect_equal(js$n$sites_phospho, 4L); expect_equal(js$genes$EGFR, c(3L, 3L, 3L, 0L, 0L, 0L)); expect_equal(js$genes$VDAC2, c(0L, 0L, 0L, 1L, 1L, 1L)); expect_equal(js$sites$TFRC_N727, c(1L, 1L))
  expect_equal(js$columns$genes, c("deg_glyco", "reach_genes_glyco", "reach_sites_glyco", "deg_phospho", "reach_genes_phospho", "reach_sites_phospho"))
})

test_that("a protein-level list resolves the entries, splits the rows and tests on the tier", {
  testthat::skip_if_not_installed("arrow")
  fixture()
  l <- ptm_list(c("EGFR", "VDAC2", "HER2", "NOTAGENE", "TFRC_N727", "EGFR_T693"), level = "protein", n_perm = 20)
  expect_s3_class(l, "ptm_list"); expect_equal(l$level, "protein")
  m <- l$members
  expect_equal(m$kind, c("gene", "gene", "gene", "gene", "site", "site"))
  expect_equal(m$gene, c("EGFR", "VDAC2", "ERBB2", NA, "TFRC", "EGFR"))
  expect_equal(m$resolved_via[3:4], c("alias of ERBB2", "not a recognized symbol"))
  expect_equal(m$site[5], "TFRC_N727"); expect_equal(m$side_named[5:6], c("glyco", "phospho"))
  expect_equal(names(m), c("query", "kind", "gene", "resolved_via", "site", "side_named", "status_glyco", "status_phospho", "n_served_as_glyco", "n_served_as_phospho", "n_reported", "n_internal", "n_neighbourhood"))
  expect_equal(m$n_reported, c(6L, 1L, 0L, NA, 2L, 1L))   # a gene entry: its protein's pairs on both sides; a site entry: that site's own
  expect_equal(m$n_served_as_glyco[1], 4L); expect_equal(m$status_glyco[2], "protein measured, no site")
  # the rows: the default view hides the two same-protein rows; the internal pair is EGFR_N196 with VDAC2_S115
  expect_equal(l$internal$pair_id, "PTMG:P00533:N196|PTMP:P45880:S115"); expect_equal(nrow(l$same_protein), 0L); expect_equal(names(l$internal), names(l$neighbourhood))
  expect_equal(l$neighbourhood$phospho_site, c("COIL_T303", "SLC20A2_S256", "ALG3_S13"))   # sorted on evidence
  expect_true(l$neighbourhood$site_named_glyco[1]); expect_false(any(l$internal$site_named_glyco)); expect_false(any(l$neighbourhood$site_named_phospho))
  expect_equal(unname(l$withheld["cis"]), 2L); expect_equal(l$provenance$n_matched, 6L); expect_equal(l$provenance$genes, c("EGFR", "VDAC2", "ERBB2", "TFRC")); expect_equal(l$provenance$level, "protein")
  expect_equal(m$n_internal[1:2], c(1L, 1L)); expect_equal(m$n_neighbourhood[c(1, 5)], c(2L, 1L)); expect_true(is.na(m$n_internal[4]))
  # connectivity for {EGFR, VDAC2, ERBB2, TFRC}: one internal pair on the tier; expected (3 + 1) * 1 / 4 = 1
  cn <- l$connectivity
  expect_equal(cn$members_in_tier, 3L); expect_equal(cn$internal_pairs, 1L); expect_equal(cn$expected_pairs, 1); expect_equal(cn$ratio, 1)
  expect_equal(cn$n_perm, 20L); expect_true(cn$p_perm > 0 && cn$p_perm <= 1); expect_true(is.finite(cn$perm_median))
  expect_equal(names(cn), c("members_in_tier", "internal_pairs", "expected_pairs", "ratio", "p_perm", "n_perm", "perm_median"))
  # partner sites: the three phospho-side partner sites each reach one of the two testable members; the tier has two glyco genes
  ps <- l$partner_sites
  expect_equal(ps$site, c("COIL_T303", "SLC20A2_S256", "ALG3_S13")); expect_equal(ps$gene, c("COIL", "SLC20A2", "ALG3")); expect_equal(unique(ps$side), "phospho")   # equal P: ordered on evidence_best
  expect_equal(ps$members_testable, rep(2L, 3)); expect_equal(ps$tier_n, rep(2L, 3)); expect_equal(ps$tier_reach, rep(1L, 3)); expect_equal(ps$members, c("TFRC", "EGFR", "EGFR"))
  expect_equal(ps$p, rep(1, 3)); expect_equal(ps$q, rep(1, 3)); expect_equal(ps$fold, rep(1, 3))
  expect_equal(ps$evidence_best[1], 21.6); expect_equal(ps$best_pair_id[1], "PTMG:P02786:N727|PTMP:P38432:T303")
  # convergent partners: the same numbers per protein
  cv <- l$convergent
  expect_equal(cv$partner, c("COIL", "SLC20A2", "ALG3")); expect_equal(unique(cv$partner_side), "phospho"); expect_equal(cv$p, rep(1, 3)); expect_equal(cv$evidence_best[cv$partner == "COIL"], 21.6)
  # composition: nothing labeled on both sides in these tables
  expect_equal(l$composition$table, c("internal", "same_protein", "neighbourhood")); expect_equal(l$composition$share, c(0, NA, 0))
  out <- utils::capture.output(print(l))
  expect_match(out[1], "^List of 4 proteins against the reported pairs \\(PTMlink "); expect_true(any(grepl("^Pairs in view: 1 within the list, 3 with a partner outside it\\.$", out)))
  expect_true(any(grepl("^Not recognized: NOTAGENE\\.$", out))); expect_true(any(grepl("^  None is reached by two or more members; 3 reached by one member each", out)))   # single pairs are counted, never ranked as findings
  for (h in c("^Members$", "^Connectivity among the members$", "^Partner sites reached by more members than expected$", "^Partner proteins reached by more members than expected$", "^  partner sites: +x\\$partner_sites$"))
    expect_true(any(grepl(h, out)), info = h)
  expect_true(all(nchar(out) <= 80L)); expect_false(any(grepl(" $", out)))
  # the row tables and the partner tables print as summaries too
  expect_match(utils::capture.output(print(l$neighbourhood))[1], "^Pairs with a partner outside the list: 3 pairs, sorted by evidence")
  expect_match(utils::capture.output(print(l$partner_sites))[1], "^3 partner sites, sorted by P; 0 at q < 0.05$")
  expect_equal(l$universe$pairs, 4L); expect_match(l$universe$rule, "q_value < 0.05")
})

test_that("a site-level list tests named sites against the tier's sites", {
  testthat::skip_if_not_installed("arrow")
  fixture()
  l <- ptm_list(c("EGFR_N603", "EGFR_N196", "TFRC_N727", "VDAC2_S115"), n_perm = 20)
  expect_equal(l$level, "site"); expect_equal(l$members$n_internal, c(0L, 1L, 0L, 1L)); expect_equal(l$members$n_neighbourhood, c(1L, 0L, 1L, 0L))
  # internal: the named glycosite with the named phosphosite; neighborhood: the two other named glycosites' rows
  expect_equal(l$internal$pair_id, "PTMG:P00533:N196|PTMP:P45880:S115"); expect_equal(l$neighbourhood$glyco_site, c("TFRC_N727", "EGFR_N603"))
  expect_equal(nrow(l$same_protein), 0L); expect_equal(l$provenance$sites, c("EGFR_N603", "EGFR_N196", "TFRC_N727", "VDAC2_S115")); expect_equal(l$provenance$genes, c("EGFR", "TFRC", "VDAC2")); expect_equal(l$provenance$level, "site")
  # partner sites: three named glycosites of the four on the tier; COIL_T303 reaches one, P = 1 - C(3,3)/C(4,3) = 0.75
  ps <- l$partner_sites
  expect_equal(ps$site, c("COIL_T303", "ALG3_S13")); expect_equal(ps$members_testable, c(3L, 3L)); expect_equal(ps$tier_n, c(4L, 4L)); expect_equal(ps$tier_reach, c(1L, 1L))
  expect_equal(ps$p, c(0.75, 0.75)); expect_equal(ps$fold, c(4 / 3, 4 / 3)); expect_equal(ps$members, c("TFRC_N727", "EGFR_N603"))
  expect_equal(names(ps), c("site", "gene", "side", "members_reached", "members_testable", "tier_reach", "tier_n", "fold", "p", "q", "n_pairs", "members", "evidence_best", "best_pair_id"))
  cv <- l$convergent; expect_equal(cv$partner, c("COIL", "ALG3")); expect_equal(cv$tier_reach, c(1L, 1L))
  expect_equal(names(cv), c("partner", "partner_side", "members_reached", "members_testable", "tier_reach", "tier_n", "fold", "p", "q", "n_pairs", "members", "evidence_best", "best_pair_id"))
  # connectivity at site level: expected = (1 + 1 + 1) * 1 / 4
  cn <- l$connectivity; expect_equal(cn$members_in_tier, 4L); expect_equal(cn$internal_pairs, 1L); expect_equal(cn$expected_pairs, 0.75)
  expect_output(print(l), "^List of 4 sites against the reported pairs")
  # a one-sided site list has no internal pair and no expectation; a gene entry at site level is refused
  l2 <- ptm_list(c("EGFR_N603", "TFRC_N727"), n_perm = 0)
  expect_equal(nrow(l2$internal), 0L); expect_equal(l2$connectivity$internal_pairs, 0L); expect_true(is.na(l2$connectivity$ratio)); expect_equal(l2$partner_sites$site, c("COIL_T303", "ALG3_S13"))
  expect_error(ptm_list(c("EGFR", "EGFR_T693"), level = "site"), "every entry must be a site key")
  expect_equal(ptm_list(c("EGFR_N603", "TFRC_N727"), level = "protein", n_perm = 0)$level, "protein")
})

test_that("a data.frame input reads its first column; the cis view adds the same-protein rows", {
  testthat::skip_if_not_installed("arrow")
  fixture()
  l <- ptm_list(data.frame(id = c("EGFR", "TFRC"), lfc = c(1.5, -2)), n_perm = 0, filters = modifyList(ptm_default_filters(), list(include_cis = TRUE)))
  expect_equal(l$members$query, c("EGFR", "TFRC")); expect_equal(l$members$gene, c("EGFR", "TFRC"))
  ps <- l$partner_sites
  expect_equal(ps$site, c("COIL_T303", "VDAC2_S115", "SLC20A2_S256", "ALG3_S13"))   # equal P and fold: ordered on evidence_best
  expect_equal(sort(l$same_protein$phospho_site), c("EGFR_T693", "TFRC_S34")); expect_equal(l$composition$share[l$composition$table == "same_protein"], 0.5)
  # a site key at protein level contributes its protein, and a pasted line splits on whitespace too
  l2 <- ptm_list(c("EGFR", "EGFR_T693", "TFRC"), level = "protein", n_perm = 0)
  expect_equal(l2$members$gene, c("EGFR", "EGFR", "TFRC")); expect_equal(l2$level, "protein")
  expect_equal(ptm_list("EGFR TFRC", n_perm = 0)$members$query, c("EGFR", "TFRC"))
})

test_that("a one-protein list gives a real hypergeometric P and no permutation; empty and unresolved lists are errors", {
  testthat::skip_if_not_installed("arrow")
  fixture()
  l <- ptm_list("EGFR", n_perm = 0)
  cv <- l$convergent; expect_equal(sort(cv$partner), c("ALG3", "SLC20A2", "VDAC2"))
  expect_equal(cv$p, rep(0.5, 3)); expect_equal(cv$fold, rep(2, 3)); expect_equal(cv$q, rep(0.5, 3))   # one testable member of two glyco genes; each partner reaches one
  expect_equal(l$partner_sites$p, rep(0.5, 3)); expect_equal(l$partner_sites$members, rep("EGFR", 3))
  expect_equal(l$connectivity$members_in_tier, 1L); expect_equal(l$connectivity$internal_pairs, 0L); expect_true(is.na(l$connectivity$ratio)); expect_true(is.na(l$connectivity$p_perm))
  expect_error(ptm_list(character(0)), "empty"); expect_error(ptm_list("NOTAGENE"), "no entry")
})


test_that("mixed lists preserve invalid and literal sites without expanding their proteins", {
  testthat::skip_if_not_installed("arrow")
  fixture()
  a <- ptm_list(c("EGFR", "TFRC_N999"), n_perm = 0)
  expect_equal(a$level, "site"); expect_true(a$mixed)
  expect_true(all(a$neighbourhood$glyco_gene != "TFRC"))
  b <- ptm_list(c("EGFR", "EGFR_N196"), n_perm = 0)
  expect_equal(length(b$sites), length(unique(b$sites)))
  expect_equal(b$members$n_neighbourhood[1], 3L)
  c <- ptm_list(c("EGFR", "TFRC_N727"), n_perm = 0)
  expect_true("TFRC_N727" %in% c$neighbourhood$glyco_site)
  expect_false("TFRC_S34" %in% c$sites)
})


test_that("site connectivity removes cross-site pairs on the same protein", {
 tab <- data.frame(glyco_site=c("EGFR_N1","TFRC_N1"),phospho_site=c("TFRC_S1","EGFR_S1"),glyco_gene=c("EGFR","TFRC"),phospho_gene=c("TFRC","EGFR"))
 tier <- ptmlink:::.ptm_tier_from(tab)
 cn <- ptmlink:::.ptm_connectivity(tier,"site",c("EGFR_N1","EGFR_S1"),n_perm=0)
 expect_equal(cn$internal_pairs,0L); expect_equal(cn$expected_pairs,0); expect_true(is.na(cn$ratio))
})


test_that("the tests count significant pairs only: a pair reported for its perturbation support alone is not reached", {
  tab <- data.frame(glyco_site = c("A_N1", "A_N1", "B_N1"), phospho_site = c("C_S1", "D_S1", "C_S1"), glyco_gene = c("A", "A", "B"), phospho_gene = c("C", "D", "C"),
                    q_value = c(0.001, 0.60, 0.04), stringsAsFactors = FALSE)
  all_rows <- ptmlink:::.ptm_tier_from(tab); sig <- ptmlink:::.ptm_tier_from(tab, significant_only = TRUE)
  expect_equal(all_rows$M, 3L); expect_equal(sig$M, 2L); expect_true(sig$significant_only)
  expect_equal(sig$sites_phospho, "C_S1")                                    # D_S1 pairs with A_N1 only at q = 0.60
  expect_equal(sig$gene_glyco["A", "deg"], 1L); expect_equal(all_rows$gene_glyco["A", "deg"], 2L)
  expect_error(ptmlink:::.ptm_tier_from(tab[, 1:4], significant_only = TRUE), "no q_value column")
  # a non-finite q is never significant
  tab$q_value[3] <- NA; expect_equal(ptmlink:::.ptm_tier_from(tab, significant_only = TRUE)$M, 1L)
})
