# ptm_explain -- the reading rules on a reported pair, as typed fields.
#
# The website states every caveat in prose. An agent needs the same facts as data: which rule applies
# to this pair, on which column, against which threshold, and the sentence the site would show. The
# rules live in one file, inst/reference.json, and this function reads them from there, so the
# package, its reference and the Python reader cannot drift apart.

.ptm_reference_raw <- function() {
  if (!is.null(.ptm_env$reference)) return(.ptm_env$reference)
  p <- system.file("reference.json", package = "ptmlink")
  if (!nzchar(p)) p <- file.path("inst", "reference.json")
  .ptm_env$reference <- jsonlite::fromJSON(p, simplifyVector = FALSE); .ptm_env$reference
}
.ptm_rules <- function() .ptm_reference_raw()$rules

# the value a rule reads from a row, as the site reads it
.ptm_rule_value <- function(rule, row, glyco_ambiguous = NA) {
  g <- function(nm) if (nm %in% names(row)) row[[nm]][1] else NA
  switch(rule$id,
    direction = NA,
    significance = g("q_value"),
    perturbation = g("anch_any"),
    same_protein = g("is_cis"),
    protein_shadow = g("rho_prot"),
    self_coupling = g("self_coupling"),
    tissue_composition = { gl <- .ptm_glyco_labelled(as.character(g("glyco_programs")), as.character(g("rep_glycans"))); pl <- nzchar(ifelse(is.na(g("phospho_programs")), "", g("phospho_programs"))); gl && pl },
    depth = g("evidence_pctile_k"),
    heterogeneity = g("I2"),
    one_cohort = g("k_cohorts"),
    without_one_cohort = g("loco_pctile_k"),
    adjusted = { lo <- g("rho_adjusted_lo"); hi <- g("rho_adjusted_hi"); if (is.na(lo) || is.na(hi)) NA else (lo > 0 | hi < 0) },
    glycoform_pooling = g("n_rep_glycans"),
    sequon_ambiguity = { v <- g("glyco_site_ambiguous"); if (is.na(v)) glyco_ambiguous else v },
    NA)
}
.ptm_rule_applies <- function(rule, value) {
  if (is.null(rule$comparator) || rule$comparator == "always") return(TRUE)
  if (rule$comparator == "informational") return(NA)
  if (length(value) != 1 || is.na(value)) return(NA)
  t <- rule$threshold
  switch(rule$comparator, "<" = value < t, "<=" = value <= t, ">" = value > t, ">=" = value >= t, "abs>=" = abs(value) >= t, "==" = isTRUE(value == t), "is_true" = isTRUE(as.logical(value)), NA)
}
.ptm_rule_statement <- function(rule, value, applies, row) {
  g <- function(nm) if (nm %in% names(row)) row[[nm]][1] else NA
  fmt <- function(x, d = 3) if (length(x) != 1 || is.na(x)) "not measured" else if (is.logical(x)) tolower(as.character(x)) else if (abs(x) < 1e-3 && x != 0) formatC(x, format = "e", digits = 1) else formatC(x, format = "f", digits = d)
  # a rule that tests a value says so when the pair has none: its sentence would otherwise read as a finding
  if (is.na(applies) && !is.null(rule$comparator) && !rule$comparator %in% c("always", "informational"))
    return(paste0("Not assessed: this pair has no value for ", paste(unlist(rule$field), collapse = ", "), "."))
  s <- rule$statement
  cohorts <- function(k) if (length(k) != 1 || is.na(k)) "an unstated number of cohorts" else if (k == 1) "1 cohort" else paste(fmt(k, 0), "cohorts")
  # depth: evidence_pctile_k ranks the pair among ALL scored site pairs with the same k, not only reported ones; the
  # share from the top is floored at the column's resolution (0.1 %) so a top pair never reads "top 0.0 %", as on the site
  extra <- switch(rule$id,
    perturbation = if (isTRUE(applies)) { ch <- c(if (isTRUE(g("anch_dose"))) "dose-response (the pair was measured together)", if (isTRUE(g("anch_gene_editing"))) "gene editing (each site on its own)", if (isTRUE(g("anch_enzymatic"))) "sialidase (each site on its own)"); paste0(" Channels: ", paste(ch, collapse = "; "), ".") } else "",
    self_coupling = if (isTRUE(applies)) {
      # the side with the larger absolute correlation with its own protein, and that correlation with its sign
      sg <- g("self_glyco"); sp <- g("self_phospho"); known <- !is.na(sg) && !is.na(sp)
      side <- if (known) (if (abs(sg) >= abs(sp)) "glycosite" else "phosphosite") else "stronger side"
      own <- if (known) (if (abs(sg) >= abs(sp)) sg else sp) else value
      if (own < 0) paste0(" The ", side, " runs against its protein, at ", fmt(own, 2), ".") else paste0(" The ", side, " tracks its protein at ", fmt(own, 2), ".")
    } else "",
    protein_shadow = if (!is.na(value)) paste0(" Protein rho ", fmt(value, 2), ".") else "",
    depth = if (!is.na(value)) paste0(" This pair is in the top ", fmt(max(0.1, 100 * (1 - value)), 1), "% of all scored site pairs measured in ", cohorts(g("k_cohorts")), ".") else "",
    heterogeneity = if (!is.na(value)) paste0(" I2 = ", fmt(value, 2), " over ", cohorts(g("k_cohorts")), ".") else "",
    without_one_cohort = if (!is.na(value)) paste0(" Without its most influential cohort", if (is.na(g("cohort_max_influence"))) "" else paste0(" (", .ptm_cohort_label(g("cohort_max_influence")), ")"), ", the pair keeps a larger share of its evidence than ", fmt(100 * value, 0), "% of reported pairs measured in ", cohorts(g("k_cohorts")), ".") else "",
    adjusted = if (!is.na(value)) paste0(" Adjusted interval [", fmt(g("rho_adjusted_lo"), 2), ", ", fmt(g("rho_adjusted_hi"), 2), "].") else "",
    glycoform_pooling = if (!is.na(value)) paste0(" ", fmt(value, 0), if (isTRUE(value == 1)) " representative glycoform" else " distinct representative glycoforms", " across the cohorts.") else "",
    significance = if (!is.na(value)) paste0(" q = ", fmt(value), "; evidence ", fmt(g("evidence"), 1), " over ", cohorts(g("k_cohorts")), " and ", if (is.na(g("n_total"))) "an unstated number of" else formatC(g("n_total"), format = "d", big.mark = ","), " tumors.") else "",
    "")
  paste0(s, extra)
}

#' The reading rules on one or more pairs, as typed fields
#'
#' One row per rule per pair: which reading rule applies, on which column, against which threshold, and
#' the sentence the website would show. The rules are those stated on the website's Help page and in the
#' package reference ([ptm_reference()]); every threshold is read from that one file.
#' @param glyco_site,phospho_site site keys; vectors ask several pairs (recycled to equal length). A pair's stable
#'   identifier (`pair_id`) may be given as `glyco_site` alone, as in [ptm_pair()].
#' @param rows optionally, a data.frame of rows of reported pairs already in hand (from [ptm_gene()], [ptm_query()] or
#'   [ptm_site()]); then the pairs are taken from it and nothing is fetched. `glyco_site` and `phospho_site`
#'   may be omitted.
#' @return a data.frame with columns `glyco_site`, `phospho_site`, `pair_id`, `rule`, `field`, `value`
#'   (as character), `threshold`, `applies` (TRUE when the rule's condition holds, FALSE when it does not, NA
#'   when the rule is informational or the value is missing), `statement`, and `filter` (the website filter
#'   or package setting that acts on the rule, if any). Attribute `note` restates the symmetry rule.
#' @examples
#' \dontrun{
#' # reads the website's data by default; optionally a local copy: ptm_set_data_dir("path/to/data")
#' e <- ptm_explain("TFRC_N727", "COIL_T303")
#' e[e$applies %in% TRUE, c("rule", "statement")]
#' ptm_explain(rows = head(ptm_gene("EGFR", "glyco"), 3))
#' }
#' @export
ptm_explain <- function(glyco_site = NULL, phospho_site = NULL, rows = NULL) {
  rules <- .ptm_rules()
  if (is.null(rows)) {
    if (is.null(glyco_site)) stop("ptm_explain() takes a glycosite and a phosphosite (ptm_explain(\"IL6ST_N131\", \"TMX1_S270\")), one pair_id, or rows = a result of ptm_gene(), ptm_site() or ptm_query()", call. = FALSE)
    glyco_site <- .ptm_chr(glyco_site, "glyco_site", "a glycosite key such as \"IL6ST_N131\", or a pair_id")
    if (!is.null(phospho_site)) phospho_site <- .ptm_chr(phospho_site, "phospho_site", "a phosphosite key such as \"TMX1_S270\"")
    if (is.null(phospho_site)) { k <- .ptm_split_pair_id(glyco_site, "ptm_explain"); glyco_site <- k$g; phospho_site <- k$p }
    if (length(glyco_site) != length(phospho_site) && length(glyco_site) != 1L && length(phospho_site) != 1L)
      stop("glyco_site and phospho_site pair up element by element: give the same number of each, or one of either (got ", length(glyco_site), " and ", length(phospho_site), ")", call. = FALSE)
    n <- max(length(glyco_site), length(phospho_site)); glyco_site <- rep_len(.ptm_canon_site(glyco_site), n); phospho_site <- rep_len(.ptm_canon_site(phospho_site), n)
    parts <- lapply(seq_len(n), function(i) { p <- ptm_pair(glyco_site[i], phospho_site[i]); if (is.null(p$row)) NULL else p$row })
    missing <- vapply(parts, is.null, logical(1))
    if (any(missing)) message(paste(paste(glyco_site[missing], phospho_site[missing], sep = " x "), collapse = ", "), if (sum(missing) == 1) " is " else " are ", .PTM_ABSENT_PAIR)
    rows <- if (any(!missing)) do.call(rbind, parts[!missing]) else NULL
    if (is.null(rows)) return(invisible(NULL))
  } else {
    if (!is.data.frame(rows) || !all(c("glyco_site", "phospho_site") %in% names(rows))) stop("rows is a result of ptm_gene(), ptm_site() or ptm_query(), with its glyco_site and phospho_site columns", call. = FALSE)
    rows <- as.data.frame(rows, stringsAsFactors = FALSE)
    if (!nrow(rows)) stop("rows is empty: there is no pair to explain", call. = FALSE)
    if (!is.null(glyco_site) && !is.null(phospho_site)) {
      rows <- rows[paste(rows$glyco_site, rows$phospho_site) %in% paste(.ptm_canon_site(glyco_site), .ptm_canon_site(phospho_site)), , drop = FALSE]
      if (!nrow(rows)) stop("rows holds none of the pairs asked for (", paste(utils::head(paste(glyco_site, phospho_site, sep = " x "), 3), collapse = ", "), ")", call. = FALSE)
    }
    if (!"n_total" %in% names(rows) && "cohorts_rho" %in% names(rows)) rows$n_total <- .ptm_n_total(rows$cohorts_rho)
    rows <- .ptm_add_ids(rows)
  }
  amb <- function(r) {
    if ("glyco_site_ambiguous" %in% names(rows)) return(NA)
    sh <- .ptm_or_absent(.ptm_shard(r$glyco_gene[1]))
    if (is.null(sh) || !is.data.frame(sh$glyco_sites) || !"glyco_site_ambiguous" %in% names(sh$glyco_sites)) return(NA)
    v <- sh$glyco_sites$glyco_site_ambiguous[sh$glyco_sites$site == r$glyco_site[1]]; if (length(v)) v[1] else NA
  }
  out <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i, , drop = FALSE]; a <- amb(r)
    do.call(rbind, lapply(rules, function(rule) {
      value <- .ptm_rule_value(rule, r, a); applies <- .ptm_rule_applies(rule, value)
      data.frame(glyco_site = r$glyco_site, phospho_site = r$phospho_site, pair_id = if ("pair_id" %in% names(r)) r$pair_id else NA_character_,
                 rule = rule$id, field = if (is.null(rule$field)) NA_character_ else paste(unlist(rule$field), collapse = ","),
                 value = if (length(value) != 1 || is.na(value)) NA_character_ else as.character(value),
                 threshold = if (is.null(rule$threshold)) NA_character_ else paste0(rule$comparator, " ", rule$threshold),
                 applies = as.logical(applies), statement = .ptm_rule_statement(rule, value, applies, r),
                 filter = if (is.null(rule$filter)) NA_character_ else rule$filter, stringsAsFactors = FALSE)
    }))
  }))
  rownames(out) <- NULL
  attr(out, "note") <- "the score is symmetric and licenses no direction; a rule that applies is a reading instruction, not a verdict on the pair"
  class(out) <- c("ptm_explain", "data.frame")
  out
}
