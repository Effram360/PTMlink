# The fixture: a two-gene, six-row release in a temporary directory, in the shapes the real files
# have. The package's own tests read it, and so does the repository's Python parity test, so the two
# readers are held to the same rows. It is internal; nothing in a user's session calls it.

.ptm_write_fixture <- function(d = file.path(tempdir(), "ptmlink_fixture"), parquet = requireNamespace("arrow", quietly = TRUE)) {
  dir.create(file.path(d, "genes"), recursive = TRUE, showWarnings = FALSE)
  wjson <- function(x, path, gz = TRUE) { js <- jsonlite::toJSON(x, dataframe = "rows", auto_unbox = TRUE, na = "null", digits = 6); if (gz) { con <- gzfile(path, "wb"); writeLines(js, con); close(con) } else writeLines(js, path) }
  # VDAC2 and COIL are partner genes with a phospho side only, so a list can hold an internal pair (EGFR with VDAC2)
  genes <- data.frame(gene = c("EGFR", "TFRC", "CD274", "VDAC2", "COIL"), status_glyco = c("usable", "usable", "protein measured, no site", "protein measured, no site", "protein measured, no site"), status_phospho = c("usable", "usable", "usable", "usable", "usable"),
                      n_cohorts_glyco_detected = c(10, 10, 0, 0, 0), n_cohorts_phospho_detected = c(10, 10, 6, 6, 9), n_scored_rows_glyco = c(7418, 3000, 0, 0, 0), n_scored_rows_phospho = c(3000, 100, 782, 500, 400),
                      n_served_as_glyco = c(4, 2, 0, 0, 0), n_served_as_phospho = c(2, 0, 0, 1, 1), stringsAsFactors = FALSE)
  # CD71 is listed twice for one symbol (an alias that is also a previous symbol): it names TFRC. G2 is an alias of two
  # symbols: it is ambiguous.
  aliases <- data.frame(query = c("ERBB1", "HER2", "P00533", "P02786", "CD71", "CD71", "G2", "G2"), resolves_to = c("EGFR", "ERBB2", "EGFR", "TFRC", "TFRC", "TFRC", "EGFR", "TFRC"),
                        kind = c("alias", "alias", "accession", "accession", "alias", "previous symbol", "alias", "alias"), stringsAsFactors = FALSE)
  wjson(list(release = "fixture", built_utc = "2026-09-02T00:00:00Z", genes = genes, aliases = aliases), file.path(d, "index.json.gz"))
  rows <- data.frame(glyco_site = c("EGFR_N603", "EGFR_N603", "EGFR_N196", "EGFR_N352", "TFRC_N727", "TFRC_N727"),
                     phospho_site = c("EGFR_T693", "ALG3_S13", "VDAC2_S115", "SLC20A2_S256", "COIL_T303", "TFRC_S34"),
                     glyco_gene = c("EGFR", "EGFR", "EGFR", "EGFR", "TFRC", "TFRC"), phospho_gene = c("EGFR", "ALG3", "VDAC2", "SLC20A2", "COIL", "TFRC"),
                     is_cis = c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE), rho_combined = c(0.40, 0.22, 0.35, 0.35, 0.38, 0.52), evidence = c(17.1, 5.3, 12.5, 10.9, 21.6, 51.0),
                     k_cohorts = c(5L, 5L, 6L, 6L, 9L, 9L), I2 = c(0.57, 0.3, 0.2, 0.1, 0.44, 0.9), q_value = c(1e-14, 8e-4, 1e-10, 1e-9, 1e-18, 1e-40),
                     anch_dose = FALSE, anch_gene_editing = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE), anch_enzymatic = FALSE, anch_any = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE),
                     q_value_anchored = c(1e-13, 1e-3, NA, NA, 1e-17, NA), evidence_pctile_k = c(0.999, 0.6, 0.95, 0.9, 0.999, 0.9999), rho_prot = c(NA, 0.12, 0.05, 0.1, 0.10, NA), evidence_prot = c(NA, 2, 1, 1, 3, NA),
                     shadow_clean = c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE), glyco_programs = c("", "", "", "", "N2H8F0S0G0=proliferation(+0.55)", "N2H8F0S0G0=proliferation(+0.55)"), phospho_programs = c("", "", "", "", "", "proliferation(+0.52)"),
                     cohorts_rho = c("BRCA:0.40:45;OV:0.30:30", "BRCA:0.20:40", "BRCA:0.35:50", "BRCA:0.35:50", "BRCA:0.40:60;PDAC:0.36:70", "BRCA:0.50:60"),
                     rep_glycans = c("BRCA:N2H8F0S0G0;OV:N2H9F0S0G0", "BRCA:N4H5F1S1G0", "BRCA:N2H8F0S0G0", "BRCA:N2H8F0S0G0", "BRCA:N2H8F0S0G0;PDAC:N2H8F0S0G0", "BRCA:N2H8F0S0G0"),
                     n_rep_glycans = c(2L, 1L, 1L, 1L, 1L, 1L),
                     rho_adjusted_lo = c(NA, NA, NA, NA, 0.09, NA), rho_adjusted_hi = c(NA, NA, NA, NA, 0.25, NA),
                     self_glyco = c(0.30, 0.30, 0.20, 0.10, 0.21, 0.21), self_phospho = c(0.90, 0.20, 0.15, 0.10, 0.77, 0.85), self_coupling = c(0.90, 0.30, 0.20, NA, 0.77, 0.85), n_cohorts_self = c(5L, 5L, 6L, NA, 9L, 9L),
                     loco_pctile_k = c(0.4, NA, 0.8, 0.7, 0.95, 0.99), cohort_max_influence = c("BRCA", NA, "BRCA", "BRCA", "PDAC", "BRCA"),
                     stringsAsFactors = FALSE)
  # site records with the identifiers the real shards carry
  gsites <- data.frame(site = c("EGFR_N603", "EGFR_N196", "EGFR_N352", "TFRC_N727"), gene = c("EGFR", "EGFR", "EGFR", "TFRC"), uniprot_acc = c("P00533", "P00533", "P00533", "P02786"),
                       position = c(603L, 196L, 352L, 727L), residue = "N", residue_ok = TRUE, window_15 = c("NCTSISGDLHILPVA", "VSSDFLSNMSMDFQN", "CKKCEGPCRKVCNGI", "VSSDFLSNMSMDFQN"),
                       topology_class = c("PASS_EXTRACELLULAR_LUMENAL", "PASS_EXTRACELLULAR_LUMENAL", "UNDECIDABLE", "PASS_EXTRACELLULAR_LUMENAL"), n_cohorts_served = c(5L, 6L, 6L, 9L), self_coupling_median = c(0.3, 0.2, 0.1, 0.21),
                       row_id = c("PTMG:P00533:N603", "PTMG:P00533:N196", "PTMG:P00533:N352", "PTMG:P02786:N727"), glyco_site_ambiguous = c(FALSE, FALSE, TRUE, FALSE), stringsAsFactors = FALSE)
  psites <- data.frame(site = c("EGFR_T693", "TFRC_S34"), gene = c("EGFR", "TFRC"), uniprot_acc = c("P00533", "P02786"), residue_shipped = c("T", "S"), position_gencode34 = c(693L, 34L), position_canonical = c(693L, 34L),
                       frame_match = "window found; same position", residue_canonical = c("T", "S"), residue_ok = TRUE, window_15 = c("RRLLQERELVEPLTP", "MMDQARSAFSNLFGG"), n_cohorts_served = c(5L, 9L), self_coupling_median = c(0.9, 0.85),
                       row_id = c("PTMP:P00533:T693", "PTMP:P02786:S34"), stringsAsFactors = FALSE)
  shard <- function(g) list(gene = g, release = "fixture", status = list(glyco = "usable", phospho = "usable"), glyco_sites = gsites[gsites$gene == g, , drop = FALSE],
                            phospho_sites = psites[psites$gene == g, , drop = FALSE],
                            n_rows_as_glyco = sum(rows$glyco_gene == g), n_rows_as_phospho = sum(rows$phospho_gene == g), rows = rows[rows$glyco_gene == g | rows$phospho_gene == g, ])
  for (g in c("EGFR", "TFRC", "VDAC2", "COIL")) wjson(shard(g), file.path(d, "genes", paste0(g, ".json.gz")))
  # k_cov as the combined release's dictionary carries it (OV given two covariates here so the interval's k is exercised)
  wjson(data.frame(cohort = c("BRCA", "OV", "PDAC"), tumour_type = c("breast", "ovarian", "pancreatic"), k_cov = c(1L, 2L, 1L), role = "tumour", stringsAsFactors = FALSE), file.path(d, "cohorts.json"), gz = FALSE)
  wjson(data.frame(composition = "N2H8F0S0G0", glycan_class = "oligomannose", stringsAsFactors = FALSE), file.path(d, "glycans.json.gz"))
  wjson(list(release = "fixture", built_utc = "2026-09-02T00:00:00Z", n_genes = 2, n_served_rows = 6, has_adjusted = FALSE, has_structure = TRUE), file.path(d, "build.json"), gz = FALSE)
  # the tier counts the website's List page reads, from the same rows of reported pairs (list.R::.ptm_tier_json)
  wjson(.ptm_tier_json(.ptm_tier_from(rows[, c("glyco_site", "phospho_site", "glyco_gene", "phospho_gene", "q_value")], significant_only = TRUE), release = "fixture", built_utc = "2026-09-02T00:00:00Z"), file.path(d, "tier_significant.json.gz"))
  # the glycoform tier for one glycosite only (TFRC_N727 has none, so ptm_pair returns NULL there)
  dir.create(file.path(d, "glycoforms"), showWarnings = FALSE)
  gf <- '{"site":"EGFR_N603","release":"fixture","columns":["glycan","rho","evidence","k","I2","q","anch"],"pairs":{"EGFR_T693":[["N2H8F0S0G0",0.4012,17.1,5,0.57,1e-14,true],["N2H9F0S0G0",0.3100,9.2,4,0.2,1e-06,false]],"ALG3_S13":[["N4H5F1S1G0",0.2200,5.3,5,0.3,8e-04,true]]}}'
  con <- gzfile(file.path(d, "glycoforms", "EGFR_N603.json.gz"), "wb"); writeLines(gf, con); close(con)
  # a structure sidecar for EGFR only: 700 modeled residues of a 700-residue chain, the site residues where
  # the sites claim them, a straight-line trace, confident until residue 650 and disordered after
  dir.create(file.path(d, "structure"), showWarnings = FALSE)
  seq <- rep("A", 700); seq[196] <- "N"; seq[352] <- "N"; seq[603] <- "N"; seq[693] <- "T"
  st <- list(gene = "EGFR", acc = "P00533", entity = "AF-P00533-F1", v = 6L, first = 1L, last = 700L, n = 700L, seq_len = 700L, plddt_global = 80,
             seq = paste(seq, collapse = ""), xyz = lapply(seq_len(700), function(i) c(round(i * 3.8, 1), 0, 0)), plddt = ifelse(seq_len(700) <= 650, 90, 40),
             tm = list(c(646L, 668L)), signal = list(c(1L, 24L)), source = "fixture", licence = "CC BY 4.0", url = "https://alphafold.ebi.ac.uk/entry/P00533", built_utc = "2026-09-12T00:00:00Z")
  con <- gzfile(file.path(d, "structure", "EGFR.json.gz"), "wb"); writeLines(jsonlite::toJSON(st, auto_unbox = TRUE, digits = 6), con); close(con)
  # a release directory beside the data directory: the reported pairs as Parquet with the identifiers, the
  # feature tables, the checksum list, the versions and the citation
  rd <- file.path(d, "release"); dir.create(rd, showWarnings = FALSE)
  ids <- function(site, tab) unname(setNames(tab$row_id, tab$site)[site])
  served <- rows; served$glyco_row_id <- ids(served$glyco_site, gsites); served$phospho_row_id <- ids(served$phospho_site, psites)
  served$phospho_row_id[is.na(served$phospho_row_id)] <- c(ALG3_S13 = "PTMP:Q92685:S13", VDAC2_S115 = "PTMP:P45880:S115", SLC20A2_S256 = "PTMP:Q08357:S256", COIL_T303 = "PTMP:P38432:T303")[served$phospho_site[is.na(served$phospho_row_id)]]
  served$pair_id <- paste0(served$glyco_row_id, "|", served$phospho_row_id)
  served$glyco_site_ambiguous <- unname(setNames(gsites$glyco_site_ambiguous, gsites$site)[served$glyco_site]) %in% TRUE
  if (isTRUE(parquet)) arrow::write_parquet(served, file.path(rd, "site_served.parquet"))
  utils::write.table(gsites, file.path(rd, "features_glyco.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(rbind(psites, data.frame(site = c("ALG3_S13", "VDAC2_S115", "SLC20A2_S256", "COIL_T303"), gene = c("ALG3", "VDAC2", "SLC20A2", "COIL"), uniprot_acc = c("Q92685", "P45880", "Q08357", "P38432"), residue_shipped = c("S", "S", "S", "T"),
                                          position_gencode34 = c(13L, 115L, 256L, 303L), position_canonical = c(13L, 115L, 256L, 303L), frame_match = "window found; same position", residue_canonical = c("S", "S", "S", "T"), residue_ok = TRUE, window_15 = "x",
                                          n_cohorts_served = c(5L, 6L, 6L, 9L), self_coupling_median = NA, row_id = c("PTMP:Q92685:S13", "PTMP:P45880:S115", "PTMP:Q08357:S256", "PTMP:P38432:T303"), stringsAsFactors = FALSE)),
                     file.path(rd, "features_phospho.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(component = c("release", "built_utc"), value = c("fixture", "2026-09-02T00:00:00Z")), file.path(rd, "versions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  sums <- list(); for (f in list.files(rd)) sums[[f]] <- unname(tools::sha256sum(file.path(rd, f)))
  writeLines(jsonlite::toJSON(sums, auto_unbox = TRUE, pretty = TRUE), file.path(rd, "SHA256SUMS.json"))
  writeLines(c("cff-version: 1.2.0", "message: fixture", 'title: "PTMlink fixture"', "type: dataset", "version: 0.0.0", "date-released: 2026-09-02", "license:", "  - MIT", "  - CC-BY-4.0", "authors:", "  - family-names: Wei", "    given-names: Effram"), file.path(rd, "CITATION.cff"))
  invisible(d)
}
