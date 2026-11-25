# Reproduce the manuscript's candidate screen from the full released site_served.parquet.
# The website's standard toggles use different thresholds and do not implement this conjunction.
# source(system.file("examples", "candidate_screen.R", package = "ptmlink"))
# rows <- arrow::read_parquet("path/to/release/site_served.parquet", col_select = tidyselect::all_of(candidate_columns))
# candidates <- candidate_screen(rows)
# perturbation_candidates <- candidates[candidates$anch_any %in% TRUE, ]
# Rows remain ranked by primary evidence; adjusted intervals are not multiplicity-corrected.
candidate_columns <- c("pair_id", "glyco_site", "phospho_site", "glyco_gene", "phospho_gene",
  "evidence", "q_value", "k_cohorts", "rho_prot", "self_coupling", "glyco_programs",
  "phospho_programs", "rho_adjusted_lo", "rho_adjusted_hi", "loco_pctile_k", "anch_any", "is_cis")

candidate_screen <- function(rows) {
  missing <- setdiff(candidate_columns, names(rows))
  if (length(missing)) stop("Missing candidate fields: ", paste(missing, collapse = ", "))
  rows <- as.data.frame(rows)
  # Missing label means no label; missing numerical evidence never passes a criterion.
  labeled <- function(x) !is.na(x) & nzchar(trimws(x))
  keep <- !(rows$is_cis %in% TRUE) &
    !is.na(rows$q_value) & rows$q_value < 0.05 &
    !is.na(rows$k_cohorts) & rows$k_cohorts >= 6 &
    !is.na(rows$rho_prot) & abs(rows$rho_prot) < 0.15 &
    !is.na(rows$self_coupling) & rows$self_coupling < 0.6 &
    !labeled(rows$glyco_programs) & !labeled(rows$phospho_programs) &
    !is.na(rows$rho_adjusted_lo) & !is.na(rows$rho_adjusted_hi) &
    ((rows$rho_adjusted_lo > 0 & rows$rho_adjusted_hi > 0) |
     (rows$rho_adjusted_lo < 0 & rows$rho_adjusted_hi < 0)) &
    !is.na(rows$loco_pctile_k) & rows$loco_pctile_k > 0.5
  out <- rows[which(keep), , drop = FALSE]
  out[order(-out$evidence, out$glyco_site, out$phospho_site), , drop = FALSE]
}
