# ptm_query -- whole-database questions over the released table of reported pairs (site_served.parquet).
#
# The shards answer one gene at a time; the Parquet answers "every reported pair with q below 1e-10 in
# at least five cohorts" in one call. It needs the arrow package and a local release directory; the
# file is never downloaded piecemeal. The rules are the package's rules: the total score sorts,
# always and only; the default view applies after the conditions; no column is recomputed.

.ptm_release_dir <- function() if (!is.null(.ptm_env$release_dir)) .ptm_env$release_dir else Sys.getenv("PTMLINK_RELEASE_DIR", "")

#' Point the package at a local release directory (the archive deposit, unpacked)
#' @param path directory holding `site_served.parquet` and, beside it, `SHA256SUMS.json`, `versions.tsv`,
#'   `features_glyco.tsv`, `features_phospho.tsv` and `CITATION.cff`
#' @param verify check `site_served.parquet` against `SHA256SUMS.json` now (the default; about ten seconds
#'   for the released 2 GB file). A mismatch is an error, never a silent pass. With `verify = FALSE` the directory
#'   is accepted as it is, and a missing Parquet is reported by [ptm_query()] when it is called.
#' @return the normalized path, invisibly
#' @export
ptm_set_release_dir <- function(path, verify = TRUE) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) stop("path is one directory name", call. = FALSE)
  if (!dir.exists(path)) stop("no such directory: ", path, call. = FALSE)
  pq <- file.path(path, "site_served.parquet")
  if (isTRUE(verify)) {
    if (!file.exists(pq)) stop("no site_served.parquet in ", path, call. = FALSE)
    sums <- file.path(path, "SHA256SUMS.json")
    if (!file.exists(sums)) stop("no SHA256SUMS.json in ", path, "; pass verify = FALSE to use the directory unverified", call. = FALSE)
    want <- jsonlite::fromJSON(sums)[["site_served.parquet"]]
    if (is.null(want)) stop("SHA256SUMS.json in ", path, " carries no entry for site_served.parquet", call. = FALSE)
    got <- unname(tools::sha256sum(pq))
    if (!identical(got, want)) stop("site_served.parquet in ", path, " does not match SHA256SUMS.json (", got, " against ", want, ")", call. = FALSE)
    .ptm_env$release_sha <- got
  } else .ptm_env$release_sha <- NA_character_
  .ptm_env$release_dir <- normalizePath(path); .ptm_env$ids_release_read <- FALSE; .ptm_env$tier <- NULL
  invisible(.ptm_env$release_dir)
}
.ptm_versions <- function() {
  d <- .ptm_release_dir(); p <- file.path(d, "versions.tsv")
  if (!nzchar(d) || !file.exists(p)) return(list())
  t <- utils::read.delim(p, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (all(c("component", "value") %in% names(t))) setNames(as.list(t$value), t$component) else list()
}
.ptm_dataset <- function() {
  d <- .ptm_release_dir()
  if (!nzchar(d)) stop("ptm_query() and ptm_list() read the full table of reported pairs (site_served.parquet), which is released with the published article; once you have it, give its release directory to ptm_set_release_dir().", call. = FALSE)
  if (!file.exists(file.path(d, "site_served.parquet"))) stop("no site_served.parquet in ", d, call. = FALSE)
  if (!requireNamespace("arrow", quietly = TRUE)) stop("ptm_query() and ptm_list() need the arrow package: install.packages(\"arrow\")", call. = FALSE)
  arrow::open_dataset(file.path(d, "site_served.parquet"))
}
.PTM_OPS <- c("<", "<=", ">", ">=", "==", "!=", "in", "not_in", "is_na", "not_na")
# A query whose conditions match more rows than this is refused unless allow_full = TRUE: every matching row is read
# into memory before the view and the limit apply.
.PTM_QUERY_MAX_ROWS <- 1000000L
# the value of a condition, as its column's type: checked, never bent to fit (11.5 is not read as 11, "yes" not as NA)
.ptm_condition_value <- function(col, ty, v) {
  shown <- function() paste(utils::head(format(v), 3), collapse = ", ")
  if (is.factor(v)) v <- as.character(v)
  if (is.null(v) || !length(v) || anyNA(v)) stop(col, ": a condition compares with a value, got none; use the operator \"is_na\" or \"not_na\" to test for a missing value", call. = FALSE)
  if (grepl("^u?int", ty) || ty %in% c("double", "float")) {
    if (!is.numeric(v)) stop(col, " holds numbers; got ", if (is.character(v)) paste0("\"", v[1], "\"") else shown(), call. = FALSE)
    if (ty %in% c("double", "float")) return(as.numeric(v))
    if (any(v != round(v))) stop(col, " holds whole numbers; got ", shown(), call. = FALSE)
    return(as.integer(v))
  }
  if (ty == "bool") { if (!is.logical(v)) stop(col, " holds TRUE or FALSE; got ", if (is.character(v)) paste0("\"", v[1], "\"") else shown(), call. = FALSE); return(v) }
  if (!is.character(v)) stop(col, " holds text; got ", shown(), call. = FALSE)
  v
}
# the genes of a query as the release spells them: an alias or accession resolved as ptm_gene() resolves it; a name
# the index does not know is said and kept as typed (it then matches no pair)
.ptm_query_genes <- function(genes) {
  out <- vapply(genes, function(g) {
    r <- tryCatch(.ptm_resolve(g), error = function(e) if (inherits(e, "ptm_ambiguous")) stop(e) else NULL)   # an ambiguous alias stays an error
    if (is.null(r)) return(.ptm_canon_gene(g))   # no index to read: the name as typed, in the release's letter case
    if (is.na(r$gene)) { message(g, " is not a recognized symbol, alias or UniProt accession: it matches no pair"); return(.ptm_canon_gene(g)) }
    r$gene
  }, character(1))
  unique(unname(out))
}
# one condition -> one arrow Expression, with the value checked against the column's type
.ptm_condition <- function(cond, schema) {
  if (!is.list(cond) || length(cond) < 2 || length(cond) > 3 || !is.character(cond[[1]]) || !is.character(cond[[2]]) || length(cond[[1]]) != 1 || length(cond[[2]]) != 1)
    stop("a condition is list(column, operator, value), e.g. list(\"q_value\", \"<\", 1e-10)", call. = FALSE)
  col <- as.character(cond[[1]]); op <- as.character(cond[[2]])
  if (!col %in% names(schema)) stop("no column ", col, " in site_served.parquet; ptm_columns() lists them", call. = FALSE)
  if (!op %in% .PTM_OPS) stop("operator ", op, " is not one of ", paste(.PTM_OPS, collapse = ", "), call. = FALSE)
  f <- arrow::Expression$field_ref(col)
  if (op == "is_na") return(arrow::Expression$create("is_null", f))
  if (op == "not_na") return(!arrow::Expression$create("is_null", f))
  if (length(cond) < 3) stop("operator ", op, " needs a value", call. = FALSE)
  ty <- schema[[col]]$type$ToString(); v <- cond[[3]]
  v <- .ptm_condition_value(col, ty, if (is.list(v)) unlist(v) else v)
  if (op %in% c("in", "not_in")) {
    e <- arrow::Expression$create("is_in", f, options = list(value_set = arrow::Array$create(v), skip_nulls = TRUE))
    return(if (op == "in") e else !e)
  }
  if (length(v) != 1) stop("operator ", op, " takes one value; use \"in\" for a set", call. = FALSE)
  switch(op, "<" = f < v, "<=" = f <= v, ">" = f > v, ">=" = f >= v, "==" = f == v, "!=" = f != v)
}
# rows with one of `values` in column a or column b (two text columns); an empty set matches no row
.ptm_either_side <- function(a, b, values, schema) {
  set <- arrow::Array$create(as.character(values[!is.na(values)]), type = arrow::utf8())
  one <- function(col) arrow::Expression$create("is_in", arrow::Expression$field_ref(col), options = list(value_set = set, skip_nulls = TRUE))
  one(a) | one(b)
}

#' Query the release: every reported pair that meets typed conditions
#'
#' Reads `site_served.parquet` through arrow and keeps the rows that meet the conditions. The result is sorted
#' on `evidence`, always: the total score is the one ranking, and no argument changes it. The default view
#' ([ptm_default_filters()]) is applied after the conditions, so `include_cis = TRUE` is needed to see
#' same-protein pairs, exactly as on the website.
#'
#' For the v1.0.0 release (7.1 million reported pairs, a 2 GB file) a query with narrow conditions takes about 10
#' seconds and a few GB of free memory. Every pair the conditions match is read into memory before the view and
#' `limit` apply, so `columns` and `limit` shorten the result, not the read: conditions that match more than one
#' million pairs take minutes and tens of GB, and are refused unless `allow_full = TRUE`.
#' @param where a list of conditions, each `list(column, operator, value)`; operators `<`, `<=`, `>`, `>=`,
#'   `==`, `!=`, `in`, `not_in` (a vector of values), `is_na`, `not_na` (no value). Conditions are combined
#'   with AND. Column names are those of [ptm_columns()]. A value has its column's type: a number for a numeric
#'   column (a whole number for a count such as `k_cohorts`), `TRUE` or `FALSE` for a logical one, text otherwise.
#' @param genes symbols, aliases or UniProt accessions: keep pairs with a gene from this set on either side
#' @param sites site keys: keep pairs with a site from this set on either side
#' @param columns the columns to return (default every column); the three identifier columns and `n_total` are
#'   always added
#' @param filters the default view, see [ptm_default_filters()]
#' @param limit return at most this many rows, from the top of the ranking; the attribute `n_view` keeps the
#'   number of rows that pass the conditions and the view
#' @param allow_full permit a query with no `where`, `genes` or `sites`, or one whose conditions match more than one
#'   million pairs (minutes, and tens of GB in memory)
#' @return a data.frame sorted on `evidence` descending, with `n_total` added, the identifier columns
#'   (`glyco_row_id`, `phospho_row_id`, `pair_id`) as the release ships them (NA where a site has no identifier: its protein has no UniProt accession, or the site does not map onto the canonical sequence; about 3 % of pairs), and attributes: `withheld`
#'   (rows the default view hid, by reason), `n_matched` (rows the conditions matched before the view), and
#'   `provenance` (release, build time, file and checksum, the conditions, the filters, the sort key, the
#'   license, the time of the query), so a result can be cited or reproduced.
#' @examples
#' \dontrun{
#' ptm_set_release_dir("path/to/release")
#' top <- ptm_query(list(list("q_value", "<", 1e-10), list("k_cohorts", ">=", 5)), limit = 100)
#' attr(top, "provenance")
#' ptm_query(genes = c("EGFR", "TFRC"), columns = c("glyco_site", "phospho_site", "evidence"))
#' ptm_query(list(list("anch_any", "==", TRUE), list("shadow_clean", "==", TRUE)))
#' }
#' @export
ptm_query <- function(where = list(), genes = NULL, sites = NULL, columns = NULL, filters = ptm_default_filters(), limit = NULL, allow_full = FALSE) {
  filters <- .ptm_filters(filters)
  if (!is.list(where)) stop("where is a list of conditions, each list(column, operator, value), e.g. list(list(\"q_value\", \"<\", 1e-10))", call. = FALSE)
  if (!is.null(genes)) genes <- .ptm_chr(genes, "genes", .PTM_GENE_EXAMPLE)
  if (!is.null(sites)) sites <- .ptm_chr(sites, "sites", .PTM_SITE_EXAMPLE)
  if (!is.null(columns)) columns <- .ptm_chr(columns, "columns", "column names; ptm_columns() lists them")
  if (!is.null(limit) && !(is.numeric(limit) && length(limit) == 1L && is.finite(limit) && limit >= 0)) stop("limit is a whole number, 0 or more, or NULL for every row", call. = FALSE)
  ds <- .ptm_dataset(); schema <- ds$schema
  if (length(where) && !is.list(where[[1]])) where <- list(where)   # one bare condition
  if (!length(where) && is.null(genes) && is.null(sites) && !isTRUE(allow_full)) stop("no condition: pass `where`, `genes` or `sites`, or allow_full = TRUE to read every reported pair", call. = FALSE)
  exprs <- lapply(where, .ptm_condition, schema = schema)
  if (!is.null(genes)) exprs <- c(exprs, list(.ptm_either_side("glyco_gene", "phospho_gene", .ptm_query_genes(genes), schema)))
  if (!is.null(sites)) exprs <- c(exprs, list(.ptm_either_side("glyco_site", "phospho_site", .ptm_canon_site(sites), schema)))
  filt <- if (length(exprs)) Reduce(`&`, exprs) else NULL
  if (!isTRUE(allow_full) && !is.null(filt)) {
    n <- arrow::Scanner$create(ds, filter = filt, projection = "glyco_site")$CountRows()
    if (n > .PTM_QUERY_MAX_ROWS) stop("the conditions match ", format(n, big.mark = ","), " reported pairs, more than ", format(.PTM_QUERY_MAX_ROWS, big.mark = ","),
      ": reading them all takes minutes and tens of GB of memory (limit shortens the answer, not the read). Add a condition, or pass allow_full = TRUE to read them anyway", call. = FALSE)
  }
  out <- .ptm_scan(ds, filt, columns, filters, limit, list(where = where, genes = genes, sites = sites))
  .ptm_as_pairs(out, "whole-release query")
}
# the columns the default view, the identifiers and n_total need; always read whatever the projection asks
.PTM_NEED <- c("glyco_site", "phospho_site", "glyco_gene", "phospho_gene", "evidence", "is_cis", "k_cohorts", "q_value", "I2", "anch_any", "shadow_clean", "self_coupling",
               "glyco_programs", "phospho_programs", "rep_glycans", "cohorts_rho", "rho_adjusted_lo", "rho_adjusted_hi", "loco_pctile_k", "glyco_row_id", "phospho_row_id", "pair_id")
# one scan of the reported pairs: the arrow filter pushed down, the default view after it, the one sort key,
# provenance on the result. ptm_query() and ptm_list() both come through here, so there is one code path
# from the Parquet to a data.frame of rows of reported pairs.
.ptm_scan <- function(ds, filt, columns, filters, limit, asked) {
  schema <- ds$schema
  proj <- if (is.null(columns)) names(schema) else { bad <- setdiff(columns, names(schema)); if (length(bad)) stop("no column(s) ", paste(bad, collapse = ", "), "; ptm_columns() lists them", call. = FALSE); unique(c(columns, intersect(.PTM_NEED, names(schema)))) }
  sc <- if (is.null(filt)) arrow::Scanner$create(ds, projection = proj) else arrow::Scanner$create(ds, filter = filt, projection = proj)
  rows <- as.data.frame(sc$ToTable(), stringsAsFactors = FALSE)
  n_matched <- nrow(rows)
  res <- .ptm_filter(rows, filters); out <- res$rows; n_view <- nrow(out)
  if (!is.null(limit) && nrow(out) > limit) out <- out[seq_len(limit), , drop = FALSE]
  if (!is.null(columns)) { keep <- unique(c(columns, "glyco_row_id", "phospho_row_id", "pair_id", if ("n_total" %in% names(out)) "n_total")); out <- out[, intersect(keep, names(out)), drop = FALSE] }
  rownames(out) <- NULL
  attr(out, "withheld") <- res$withheld; attr(out, "n_matched") <- n_matched; attr(out, "n_view") <- n_view
  attr(out, "provenance") <- c(.ptm_provenance(), asked, list(filters = filters, limit = limit, sorted_on = "evidence (descending); the total score is the only ranking",
    n_matched = n_matched, n_in_view = n_view, n_returned = nrow(out), licence = "data CC BY 4.0", queried_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    note = "the score is symmetric and licenses no direction; adjusted, parent and self-coupling columns are annotations, never a ranking"))
  out
}
# the release in use: stamp, build time, the file and its checksum
.ptm_provenance <- function() {
  v <- .ptm_versions()
  list(inference = .ptm_inference(v), release = .ptm_public_label(v$release, v$public_release), built_utc = if (!is.null(v$built_utc)) v$built_utc else NA_character_,
       file = file.path(.ptm_release_dir(), "site_served.parquet"), sha256 = if (!is.null(.ptm_env$release_sha)) .ptm_env$release_sha else NA_character_)
}
