# ptmlink -- the R client for the PTMlink release.
#
# One database, three ways to reach it. The per-gene data files the website reads (its data/ directory:
# index.json.gz, genes/<GENE>.json.gz, glycoforms/, structure/, cohorts.json, glycans.json.gz,
# build.json) come from a local directory (ptm_set_data_dir()) or by checksum-verified download
# on first use (ptm_set_data_url()); the released table of reported pairs (release/<v>/site_served.parquet)
# comes from a local release directory (ptm_set_release_dir()) and answers whole-database
# questions through ptm_query(). Nothing in this package recomputes a score: every number is
# read from the release. The reading rules every function applies are stated once, in
# inst/reference.json, and ptm_explain() reads them from there.
#
# The data policy in one sentence: the package installs and checks with no release present;
# ptm_status() works on the small index alone and is the documented first call.
#
# Files: ptmlink.R (data sources, the shard readers, the default view), query.R (the Parquet
# tier), explain.R (the reading rules on a pair), sites.R (reported sites with their modeled
# position), reference.R (the machine-readable reference, the columns, the citation),
# fixture.R (the fixture the tests and the Python parity test share).

#' ptmlink: query the PTMlink database from R
#'
#' PTMlink ranks pairs of one N-glycosylation site and one phosphorylation site by how strongly the two co-vary
#' across CPTAC tumor cohorts. The package reads the same data as the PTMlink website and needs no setup:
#' `library(ptmlink); ptm_gene("EGFR")`.
#'
#' @section Where to start:
#' * [ptm_status()] says how a gene name resolves and what was measured for it.
#' * [ptm_gene()] and [ptm_site()] list the ranked partners of a gene or a site.
#' * [ptm_pair()] gives one pair with its estimate in every cohort; [ptm_explain()] gives the cautions that apply.
#' * [ptm_network()] and [ptm_plot_network()] give the partners of a gene as a network.
#' * [ptm_query()] and [ptm_list()] read the full table of reported pairs, which is released with the article.
#' * [ptm_columns()] defines every column; [ptm_default_filters()] describes the default view; [ptm_cite()] says how to cite.
#'
#' @section Reading a result:
#' Pairs are ranked by `evidence`, the -log10 P of the estimate pooled across cohorts. A correlation has no
#' direction. Perturbation support can add a pair to the reported pairs and never raises its rank. The default
#' view shows pairs between two different proteins that are significant at q < 0.05 in at least two cohorts.
#'
#' `vignette("ptmlink")` walks through every function.
#' @keywords internal
"_PACKAGE"

#' @importFrom stats setNames median
#' @importFrom utils download.file read.delim write.table citation
#' @importFrom tools sha256sum R_user_dir
NULL

.ptm_env <- new.env(parent = emptyenv())
.ptm_env$data_dir <- NULL
.ptm_env$base_url <- ""   # set by ptm_set_data_url(); PTMLINK_DATA_URL is read at call time, never cached at load
.ptm_env$release_dir <- NULL
.ptm_env$cache <- list()
.ptm_env$remote <- NULL
.ptm_env$ids <- character(0)   # site key -> row_id, learned from every shard read and from the release's feature tables

#' Point the package at the release's data address for on-demand, checksum-verified download
#' @param url the address that serves `index.json.gz`, `genes/`, `SHA256SUMS.json`
#' @return the address, invisibly
#' @export
ptm_set_data_url <- function(url) {
  if (!is.character(url) || length(url) != 1L || is.na(url) || !grepl("^[A-Za-z][A-Za-z0-9+.-]*://", url)) stop("url is one web address, such as \"https://ptmlink.biomarkercenter.org/data\"", call. = FALSE)
  .ptm_env$data_dir <- NULL; .ptm_env$base_url <- sub("/$", "", url); .ptm_env$remote <- NULL; .ptm_env$cache <- list(); .ptm_env$ids <- character(0); .ptm_env$ids_package_read <- FALSE; invisible(url)
}
# With no source chosen, the package reads the public website, so library(ptmlink); ptm_gene("EGFR") works as is.
.PTM_DEFAULT_URL <- "https://ptmlink.biomarkercenter.org/data"
.ptm_base_url <- function() if (nzchar(.ptm_env$base_url)) .ptm_env$base_url else Sys.getenv("PTMLINK_DATA_URL", .PTM_DEFAULT_URL)

#' Point the package at a local release data directory (the website's `data/`)
#' @param path directory holding `index.json.gz`, `genes/`, `cohorts.json`, `build.json`
#' @return the normalized path, invisibly
#' @export
ptm_set_data_dir <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) stop("path is one directory name", call. = FALSE)
  if (!dir.exists(path)) stop("no such directory: ", path, call. = FALSE)
  if (!file.exists(file.path(path, "index.json.gz"))) stop(path, " is not a PTMlink data directory: it has no index.json.gz. The right directory also holds build.json and genes/.", call. = FALSE)
  .ptm_env$data_dir <- normalizePath(path); .ptm_env$cache <- list(); .ptm_env$ids <- character(0); .ptm_env$ids_package_read <- FALSE; invisible(.ptm_env$data_dir)
}

# THREE KINDS OF FAILURE, KEPT APART. A file the data source does not list (`ptm_absent`) is a fact about the
# release: a gene with no reported pair has no data file. A host that cannot be reached (`ptm_unreachable`) or a
# download that does not match the checksum list (`ptm_checksum`) is a fact about the connection, and is never
# reported as the absence of a gene or a pair. Callers catch `ptm_absent` alone; the other two reach the user.
.ptm_error <- function(msg, class) stop(structure(class = c(class, "ptm_error", "error", "condition"), list(message = msg, call = NULL)))
.ptm_absent <- function(rel) .ptm_error(paste0(rel, " is not in this data source"), "ptm_absent")
.ptm_host <- function(base) { h <- sub("^[A-Za-z]+://([^/]*).*$", "\\1", base); if (identical(h, base) || !nzchar(h)) base else h }
.ptm_unreachable <- function(base, detail, answered = FALSE) .ptm_error(if (answered)
  paste0(.ptm_host(base), " answered ", detail, ". If the address is wrong, set it with ptm_set_data_url(): it is the folder that holds SHA256SUMS.json. Otherwise try again in a few minutes, or use a local copy of the data with ptm_set_data_dir().")
  else paste0("could not reach ", .ptm_host(base), if (nzchar(detail)) paste0(" (", detail, ")"),
  ". Check the internet connection and try again, or use a local copy of the data with ptm_set_data_dir()."), "ptm_unreachable")

# INPUT CHECKS, worded once. A name argument is one or more non-empty strings; anything else is refused with the
# argument's name and an example, before any file is read.
.PTM_GENE_EXAMPLE <- "a gene symbol, alias or UniProt accession such as \"EGFR\""
.PTM_SITE_EXAMPLE <- "a site key such as \"EGFR_N603\" or \"EGFR_T693\""
.ptm_chr <- function(x, arg, example, one = FALSE) {
  if (is.factor(x)) x <- as.character(x)
  if (!length(x)) stop(arg, " is empty: give ", example, call. = FALSE)
  if (all(is.na(x))) stop(arg, " is NA: give ", example, call. = FALSE)
  if (!is.character(x)) stop(arg, " is ", class(x)[1], ", not text: give ", example, call. = FALSE)
  x <- trimws(x)
  if (anyNA(x) || any(!nzchar(x))) stop(arg, " holds a missing or empty entry: give ", example, call. = FALSE)
  if (one && length(x) != 1L) stop(arg, " takes one value here, got ", length(x), call. = FALSE)
  x
}
# The view's settings, checked: a partial list is completed from the defaults, an unknown name or a value of the
# wrong kind is refused. Every function that takes `filters` comes through here.
.ptm_filters <- function(f) {
  def <- unclass(ptm_default_filters())
  if (is.null(f)) return(structure(def, class = "ptm_filters"))
  if (!is.list(f)) stop("filters is a list of settings, such as list(min_k = 1); ?ptm_default_filters describes them", call. = FALSE)
  f <- unclass(f); nms <- names(f)
  if (length(f) && (is.null(nms) || any(!nzchar(nms)))) stop("every filter needs its name, as in list(min_k = 1)", call. = FALSE)
  bad <- setdiff(nms, names(def))
  if (length(bad)) stop("unknown filter ", paste0("\"", bad, "\"", collapse = ", "), "; the filters are ", paste(names(def), collapse = ", "), call. = FALSE)
  flags <- names(def)[vapply(def, is.logical, logical(1))]
  for (nm in nms) {
    v <- f[[nm]]
    if (nm %in% flags) {
      if (!(is.logical(v) && length(v) == 1L && !is.na(v))) stop(nm, " is TRUE or FALSE", call. = FALSE)
    } else if (nm == "min_k") {
      if (!(is.numeric(v) && length(v) == 1L && is.finite(v) && v >= 1 && v == round(v))) stop("min_k is a whole number of cohorts, 1 or more", call. = FALSE)
      v <- as.integer(v)
    } else {
      if (is.null(v) || (length(v) == 1L && is.na(v))) v <- NA_real_
      if (!(is.numeric(v) && length(v) == 1L) || (is.finite(v) && (v < 0 || v > 1))) stop(nm, " is a number between 0 and 1, or NA for no limit", call. = FALSE)
      v <- as.numeric(v)
    }
    def[[nm]] <- v
  }
  structure(def, class = "ptm_filters")
}
# download.file reports the reason as a warning and then fails with "cannot open URL"; keep the reason
.ptm_fetch <- function(url, dest, base) {
  why <- character(0)
  ok <- tryCatch(withCallingHandlers({ utils::download.file(url, dest, mode = "wb", quiet = TRUE); TRUE },
                                     warning = function(w) { why <<- c(why, conditionMessage(w)); invokeRestart("muffleWarning") }),
                 error = function(e) { why <<- c(why, conditionMessage(e)); FALSE })
  if (!isTRUE(ok)) {
    st <- grep("status was '", why, value = TRUE, fixed = TRUE)   # the reason may itself hold an apostrophe ("Couldn't ...")
    detail <- if (length(st)) sub("^.*status was '(.*)'.*$", "\\1", st[1]) else if (length(why)) why[1] else ""
    unlink(dest); .ptm_unreachable(base, detail, answered = grepl("^[0-9]{3}( |$)", detail))   # an HTTP status, not a failed connection
  }
  invisible(TRUE)
}

.ptm_remote <- function(base) {
  if (!is.null(.ptm_env$remote) && identical(.ptm_env$remote$base, base)) return(.ptm_env$remote)
  # One fresh manifest per source selection/session; namespace disk files by URL.
  stamp <- tempfile(); writeLines(base, stamp, useBytes = TRUE)
  id <- unname(tools::sha256sum(stamp)); unlink(stamp)
  root <- file.path(tools::R_user_dir("ptmlink", "cache"), "data", id)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  manifest <- file.path(root, "SHA256SUMS.json")
  tmp <- tempfile(tmpdir = root); on.exit(unlink(tmp), add = TRUE)
  got <- tryCatch(.ptm_fetch(paste0(base, "/SHA256SUMS.json"), tmp, base), ptm_unreachable = function(e) e)
  if (inherits(got, "ptm_unreachable")) {
    # offline with a cache: the checksum list saved by the last session still verifies every cached file
    old <- if (file.exists(manifest)) tryCatch(jsonlite::fromJSON(manifest), error = function(e) NULL) else NULL
    if (!is.list(old) || !length(old)) stop(got)
    warning("could not reach ", .ptm_host(base), "; using the checksum list cached on ", format(file.mtime(manifest), "%Y-%m-%d"),
            " (files not cached yet cannot be read)", call. = FALSE)
    .ptm_env$remote <- list(base = base, root = root, sums = old, manifest = manifest)
    return(.ptm_env$remote)
  }
  sums <- tryCatch(jsonlite::fromJSON(tmp), error = function(e) NULL)
  if (!is.list(sums) || !length(sums)) stop("invalid checksum manifest at ", base, call. = FALSE)
  if (!file.rename(tmp, manifest)) stop("cannot save checksum manifest", call. = FALSE)
  .ptm_env$remote <- list(base = base, root = root, sums = sums, manifest = manifest)
  .ptm_env$remote
}
# Does the data source carry this file? Asks the checksum list (remote) or the directory (local); fetches nothing else.
.ptm_has <- function(rel) {
  if (!is.null(.ptm_env$data_dir)) return(file.exists(file.path(.ptm_env$data_dir, rel)))
  base <- .ptm_base_url(); if (!nzchar(base)) return(FALSE)
  !is.null(.ptm_remote(base)$sums[[rel]])
}

.ptm_path <- function(rel) {
  # a name built from what the caller typed (a gene, a site) never leaves the data source
  if (grepl("^[/\\\\]|\\\\", rel) || ".." %in% strsplit(rel, "/", fixed = TRUE)[[1]]) .ptm_absent(rel)
  if (!is.null(.ptm_env$data_dir)) { p <- file.path(.ptm_env$data_dir, rel); if (!file.exists(p)) .ptm_absent(rel); return(p) }
  base <- .ptm_base_url()
  if (!nzchar(base)) stop("no data source: call ptm_set_data_dir(<data directory>) or ptm_set_data_url(<address>)", call. = FALSE)
  remote <- .ptm_remote(base)
  if (identical(rel, "SHA256SUMS.json")) return(remote$manifest)
  want <- remote$sums[[rel]]
  if (is.null(want)) .ptm_absent(rel)
  local <- file.path(remote$root, rel)
  dir.create(dirname(local), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(local) || !identical(unname(tools::sha256sum(local)), want)) {
    tmp <- tempfile(tmpdir = dirname(local)); on.exit(unlink(tmp), add = TRUE)
    .ptm_fetch(paste0(base, "/", rel), tmp, base)
    tryCatch(.ptm_verify(tmp, rel, base), ptm_checksum = function(e) {
      # the site may have been updated since this session read its checksum list: read the list again, once, and
      # accept the file only if the new list names it
      .ptm_env$remote <- NULL
      if (identical(.ptm_remote(base)$sums[[rel]], want)) stop(e)
      .ptm_env$cache <- list()
      .ptm_fetch(paste0(base, "/", rel), tmp, base)
      .ptm_verify(tmp, rel, base)
    })
    if (!file.rename(tmp, local)) stop("cannot save verified download ", rel, call. = FALSE)
  }
  local
}
# Cached and newly downloaded files are checked against the current session manifest.
.ptm_verify <- function(local, rel, base) {
  sums <- .ptm_remote(base)$sums
  if (is.null(sums[[rel]])) { unlink(local); .ptm_absent(rel) }
  if (!identical(unname(tools::sha256sum(local)), sums[[rel]])) {
    unlink(local)
    .ptm_error(paste0("checksum mismatch for ", rel, " -- download removed. The site may be in the middle of an update; try again in a few minutes."), "ptm_checksum")
  }
  invisible(TRUE)
}
# Run `expr`; a file the data source does not carry gives `absent` instead. Every other failure propagates.
.ptm_or_absent <- function(expr, absent = NULL) tryCatch(expr, ptm_absent = function(e) absent)
.ptm_json <- function(rel) {
  key <- if (is.null(.ptm_env$data_dir)) paste(.ptm_base_url(), rel, sep = "|") else rel
  if (!is.null(.ptm_env$cache[[key]])) return(.ptm_env$cache[[key]])
  p <- .ptm_path(rel)
  txt <- if (grepl("\\.gz$", p)) paste(readLines(gzfile(p), warn = FALSE), collapse = "\n") else paste(readLines(p, warn = FALSE), collapse = "\n")
  v <- jsonlite::fromJSON(txt, simplifyVector = TRUE)
  .ptm_env$cache[[key]] <- v
  if (grepl("^genes/", rel)) .ptm_learn_ids(v)
  v
}
# Release metadata, never the release name, says the method. Two revisions are read: student_normal_bh (Student tails,
# signed normal pooling, q_value BH, the only correction) and the legacy Fisher-normal/BH of the releases written before
# the marker existed. Any other marker reads as legacy, as the release build itself reads it.
.ptm_inference <- function(meta = list()) {
  if (!is.null(meta[["inference"]])) meta <- meta[["inference"]]
  modern <- identical(meta$inference_revision, "student_normal_bh")
  out <- list(inference_revision = if (modern) "student_normal_bh" else "legacy_fisher_bh",
              p_method = if (modern) "signed_student_normal_weighted" else "fisher_normal",
              q_method = "BH", q_family = "all_scored_site_pairs",
              anchored_q_family = "perturbation_anchored_rows_within_each_layer",
              served_rule = "q_value<0.05 OR anch_any")
  for (nm in c("q_family", "anchored_q_family", "glycoform_q_family", "cohort_q_family"))
    if (!is.null(meta[[nm]])) out[[nm]] <- meta[[nm]]
  out
}
.ptm_shard <- function(gene) {
  sh <- .ptm_json(paste0("genes/", gene, ".json.gz"))
  sh$inference <- .ptm_inference(sh); attr(sh$rows, "inference") <- sh$inference
  sh
}
# A gene's AlphaFold record, found as the website finds it: through structure/bundles.json (gene -> bundle file, each
# bundle an object keyed by gene) when the data source carries the index, else structure/<GENE>.json.gz (a data
# directory written before bundles). NULL when the source has no model for the gene; a network or checksum failure
# propagates.
.ptm_structure <- function(gene) {
  if (.ptm_has("structure/bundles.json")) {
    idx <- .ptm_json("structure/bundles.json"); path <- idx$genes[[gene]]
    if (is.null(path)) return(NULL)
    b <- .ptm_json(path); return(b[[gene]])
  }
  .ptm_or_absent(.ptm_json(paste0("structure/", gene, ".json.gz")))
}

# STABLE IDENTIFIERS. Every site record in a shard carries its row_id (PTMG:<accession>:N<position>,
# PTMP:<accession>:<residue><canonical position>); a pair's id is the two joined by "|", as
# site_served.parquet ships it. A shard knows the ids of ITS OWN sites only, so a partner site's id comes from
# the map of every reported site the package ships (inst/extdata/site_ids.tsv.gz, written when the release is built),
# from the partner's shard when that is read, or from the release's feature tables when ptm_set_release_dir() has
# been called. A site without a UniProt accession has no id, and its pairs keep pair_id NA. ptm_query() reads the
# ids from the Parquet directly.
.ptm_learn_ids <- function(sh) {
  for (tab in list(sh$glyco_sites, sh$phospho_sites)) {
    if (is.data.frame(tab) && all(c("site", "row_id") %in% names(tab)) && nrow(tab)) {
      ok <- !is.na(tab$site) & !is.na(tab$row_id); .ptm_env$ids[tab$site[ok]] <- tab$row_id[ok]
    }
  }
  invisible(NULL)
}
.ptm_ids_from_release <- function() {
  d <- .ptm_release_dir(); if (!nzchar(d) || isTRUE(.ptm_env$ids_release_read)) return(invisible(NULL))
  for (f in c("features_glyco.tsv", "features_phospho.tsv")) {
    p <- file.path(d, f)
    if (file.exists(p)) { t <- utils::read.delim(p, stringsAsFactors = FALSE, quote = "", comment.char = ""); if (all(c("site", "row_id") %in% names(t))) { ok <- !is.na(t$site) & !is.na(t$row_id) & nzchar(t$row_id); .ptm_env$ids[t$site[ok]] <- t$row_id[ok] } }
  }
  .ptm_env$ids_release_read <- TRUE; invisible(NULL)
}
.ptm_ids_from_package <- function() {
  if (isTRUE(.ptm_env$ids_package_read)) return(invisible(NULL))
  p <- system.file("extdata", "site_ids.tsv.gz", package = "ptmlink")
  if (nzchar(p)) {
    t <- utils::read.delim(gzfile(p), colClasses = "character", quote = "", comment.char = "")
    new <- !t$site %in% names(.ptm_env$ids); .ptm_env$ids[t$site[new]] <- t$row_id[new]
  }
  .ptm_env$ids_package_read <- TRUE; invisible(NULL)
}
.ptm_add_ids <- function(rows) {
  if (!nrow(rows)) { for (nm in c("glyco_row_id", "phospho_row_id", "pair_id")) if (!nm %in% names(rows)) rows[[nm]] <- character(0); return(rows) }
  .ptm_ids_from_package(); .ptm_ids_from_release()
  look <- function(keys) { v <- unname(.ptm_env$ids[keys]); v[is.na(keys)] <- NA_character_; v }
  if (!"glyco_row_id" %in% names(rows) || anyNA(rows$glyco_row_id)) { g <- look(rows$glyco_site); if ("glyco_row_id" %in% names(rows)) g[!is.na(rows$glyco_row_id)] <- rows$glyco_row_id[!is.na(rows$glyco_row_id)]; rows$glyco_row_id <- g }
  if (!"phospho_row_id" %in% names(rows) || anyNA(rows$phospho_row_id)) { p <- look(rows$phospho_site); if ("phospho_row_id" %in% names(rows)) p[!is.na(rows$phospho_row_id)] <- rows$phospho_row_id[!is.na(rows$phospho_row_id)]; rows$phospho_row_id <- p }
  rows$pair_id <- ifelse(is.na(rows$glyco_row_id) | is.na(rows$phospho_row_id), NA_character_, paste0(rows$glyco_row_id, "|", rows$phospho_row_id))
  rows
}

#' The release stamp
#' @return the release's build record as a list (`release`, `released`, the public release date from the package's
#'   citation file, `built_utc`, the date the data were built, the counts and the layers the build carries)
#' @export
ptm_release <- function() { b <- .ptm_json("build.json"); b$inference <- .ptm_inference(b); b <- .ptm_public(b); b$released <- .ptm_release_date(); class(b) <- "ptm_release"; b }
# the public release date, as the citation files state it
.ptm_release_date <- function() {
  p <- system.file("extdata", "ptmlink_citation.cff", package = "ptmlink"); if (!nzchar(p)) p <- file.path("inst", "extdata", "ptmlink_citation.cff")
  tryCatch(.ptm_read_cff(p)$date, error = function(e) NA_character_)
}

# The label a reader sees. The build record keeps the internal build label; the twelve-cohort build v1.2 is the first
# public release, v1.0.0, as on the website. A build record that states `public_release` wins. The internal label is
# kept as `build_release`.
.PTM_PUBLIC_LABEL <- c(v1.2 = "v1.0.0")
.ptm_public_label <- function(rel, public = NULL) {
  if (!is.null(public) && length(public) == 1L && !is.na(public) && nzchar(public)) return(public)
  if (is.null(rel) || length(rel) != 1L || is.na(rel)) return(NA_character_)
  if (rel %in% names(.PTM_PUBLIC_LABEL)) unname(.PTM_PUBLIC_LABEL[[rel]]) else rel
}
.ptm_public <- function(b) {
  if (is.list(b) && !is.null(b$release)) { b$build_release <- b$release; b$release <- .ptm_public_label(b$release, b$public_release) }
  b
}

#' The cohort dictionary and the glycan table
#' @return a data.frame read from the release. For the cohorts: one row per cohort with its `role`, "tumour" (a tumor cohort, pooled
#'   into every score) or "normal" (normal tissue of a tumor cohort, shown beside it, never pooled or counted), and
#'   `k_cov`, the covariate count of its pooling weight n - 3 - k_cov; a release whose file has no `role` is
#'   all tumor. `label` is the article's display name (GC for the analysis key STAD, ccRCC, Rare kidney cancers);
#'   arguments and the release files use `cohort`.
#' @export
ptm_cohorts <- function() {
  d <- .ptm_json("cohorts.json")
  if (is.data.frame(d) && nrow(d) && !"role" %in% names(d)) d$role <- "tumour"
  d <- .ptm_display_cohorts(d)
  if (is.data.frame(d)) class(d) <- c("ptm_cohorts", "data.frame")
  d
}
# The article's display names, as on the website: GC for the gastric cohort, ccRCC, rare kidney
# cancers, and NAT for normal adjacent tissue. `cohort` keeps the analysis key the release files and every argument use.
.PTM_COHORT_LABEL <- c(CCRCC = "ccRCC", RARE_KIDNEY = "Rare kidney cancers", STAD = "GC")
.PTM_TUMOUR_TYPE <- c(STAD = "gastric cancer", RARE_KIDNEY = "non-clear-cell renal cancer")
.ptm_cohort_label <- function(x) { x <- as.character(x); ifelse(!is.na(x) & x %in% names(.PTM_COHORT_LABEL), unname(.PTM_COHORT_LABEL[x]), x) }
.ptm_display_cohorts <- function(d) {
  if (!is.data.frame(d) || !nrow(d) || !"cohort" %in% names(d)) return(d)
  lab <- .ptm_cohort_label(d$cohort)
  d$label <- if ("label" %in% names(d)) ifelse(is.na(d$label) | !nzchar(d$label), lab, d$label) else lab
  if ("tumour_type" %in% names(d)) {
    base <- ifelse(d$cohort %in% names(.PTM_TUMOUR_TYPE), unname(.PTM_TUMOUR_TYPE[d$cohort]),
                   sub("\\s*\\(rare kidney histologies\\)", "", sub(",\\s*normal adjacent tissue( \\(NAT\\))?$", "", d$tumour_type)))
    d$tumour_type <- ifelse(d$role %in% "normal", paste0(base, ", normal adjacent tissue (NAT)"), base)
  }
  d[, c("cohort", "label", setdiff(names(d), c("cohort", "label"))), drop = FALSE]
}
# cohort -> k_cov of the tumor cohorts the dictionary states; a dictionary without k_cov, or no dictionary,
# gives every cohort k = 1 (the n - 4 weight)
.ptm_kcov <- function() {
  d <- .ptm_or_absent(.ptm_json("cohorts.json"))
  if (!is.data.frame(d) || !all(c("cohort", "k_cov") %in% names(d))) return(setNames(numeric(0), character(0)))
  if ("role" %in% names(d)) d <- d[is.na(d$role) | d$role != "normal", , drop = FALSE]
  ok <- !is.na(d$k_cov)
  setNames(as.numeric(d$k_cov[ok]), d$cohort[ok])
}
#' @rdname ptm_cohorts
#' @export
ptm_glycans <- function() { d <- .ptm_json("glycans.json.gz"); if (is.data.frame(d)) class(d) <- c("ptm_glycans", "data.frame"); d }

#' The default view, as the website applies it
#'
#' Every function that takes `filters` accepts this list, a changed copy of it
#' (`modifyList(ptm_default_filters(), list(min_k = 1))`), or a list naming only the settings to change
#' (`list(min_k = 1)`). An unknown name or a value of the wrong kind is an error.
#' @return a named list of filter settings; `significant_only` is TRUE, as on the website:
#'   only pairs at q < 0.05 on the total score are returned, and a caller who
#'   wants every reported pair (the perturbation arm included) sets it to FALSE itself -- the
#'   filters combine as a plain conjunction and nothing overrides it; `max_i2` is an optional ceiling on the heterogeneity
#'   statistic (NA, the default, applies none; the website offers 0.5); `adjusted_excludes_zero`
#'   keeps only pairs whose composition-adjusted 95 % interval excludes zero (a filter, never a sort);
#'   `shadow_clean_only` keeps only pairs with weak parent-protein co-variation (`shadow_clean`, i.e.
#'   |rho_prot| < 0.15 and not same-protein; a pair with no shadow measured is withheld);
#'   `hide_self_coupled` drops pairs whose stronger side tracks its own protein at |rho| >= 0.6
#'   (`self_coupling`; a pair with no value is kept); `min_loco_pctile` is an optional floor on
#'   `loco_pctile_k`, the percentile among reported pairs of the same `k_cohorts` of the share of evidence a
#'   pair keeps without the cohort whose removal costs it most (NA, the default, applies none; single-cohort
#'   pairs have no value and are withheld when it is set)
#' @export
ptm_default_filters <- function() structure(list(include_cis = FALSE, min_k = 2L, hide_composition = FALSE, anchored_only = FALSE, significant_only = TRUE, max_i2 = NA_real_, adjusted_excludes_zero = FALSE, shadow_clean_only = FALSE, hide_self_coupled = FALSE, min_loco_pctile = NA_real_), class = "ptm_filters")

.ptm_index <- function() .ptm_json("index.json.gz")
# a query that resolves to no gene: a site key is pointed at ptm_site(), anything else at ptm_status()
.ptm_not_symbol <- function(q) {
  if (grepl("^.+_[NSTY][0-9]+$", toupper(trimws(q)))) stop(q, " is a site key, not a gene; ptm_site(\"", q, "\") lists its partners", call. = FALSE)
  stop(q, " is not a recognized symbol, alias or UniProt accession; ptm_status(\"", q, "\") says how symbols resolve", call. = FALSE)
}
.ptm_resolve <- function(gene) {
  idx <- .ptm_index(); g <- idx$genes
  up <- toupper(trimws(gene))
  hit <- g$gene[!is.na(g$gene) & toupper(g$gene) == up]
  if (length(hit) == 1) return(list(gene = hit, via = "current symbol"))
  a <- idx$aliases; a <- a[!is.na(a$query) & !is.na(a$resolves_to) & toupper(a$query) == up, , drop = FALSE]
  # an alias listed twice for one symbol (two letter cases, or an alias that is also a previous symbol) names that symbol
  to <- sort(unique(a$resolves_to))
  if (length(to) == 1) return(list(gene = to, via = paste(sub(" of$", "", a$kind[1]), "of", to)))
  if (length(to) > 1) .ptm_error(paste0(gene, " is ambiguous: it is an alias of ", paste(c(paste(to[-length(to)], collapse = ", "), to[length(to)]), collapse = " and "), "; give one of these symbols"), "ptm_ambiguous")
  list(gene = NA_character_, via = "unknown")
}
# A gene name as the index spells it (C10orf71, not C10ORF71); an unknown name comes back in upper case.
.ptm_canon_gene <- function(gene) {
  up <- toupper(trimws(as.character(gene)))
  idx <- tryCatch(.ptm_index(), error = function(e) NULL)
  if (is.null(idx)) return(up)
  g <- idx$genes$gene; g <- g[!is.na(g)]; hit <- g[match(up, toupper(g))]
  # an alias or accession that names one symbol is read as that symbol (HER2_N530 is ERBB2_N530)
  miss <- which(is.na(hit))
  if (length(miss)) {
    a <- idx$aliases; a <- a[!is.na(a$query) & !is.na(a$resolves_to), , drop = FALSE]; aq <- toupper(a$query)
    for (i in miss) { to <- unique(a$resolves_to[aq == up[i]]); if (length(to) == 1) hit[i] <- to }
  }
  ifelse(is.na(hit), up, hit)
}
# A typed site key in the release's own letter case: the gene part as the index spells it, the residue in upper
# case. A site identifier (PTMG:..., PTMP:...) and anything that is not a site key come back in upper case.
.ptm_canon_site <- function(site) {
  up <- toupper(trimws(as.character(site)))
  m <- regmatches(up, regexec("^(.+)_([NSTY][0-9]+)$", up))
  vapply(seq_along(up), function(i) {
    if (!length(m[[i]]) || grepl("^PTM[GP]:", up[i])) up[i] else paste0(.ptm_canon_gene(m[[i]][2]), "_", m[[i]][3])
  }, "")
}

#' Why a gene returns what it returns -- the documented first call
#' @param gene one or more symbols, aliases or UniProt accessions
#' @return one row per query: the resolved symbol, how it resolved, the status word on each
#'   side, and the counts behind it; `status_*` is one of "usable", "measured, no reported pair",
#'   "sites detected, not scored", "protein measured, no site", "absent from every matrix",
#'   or "not a recognized symbol"
#' @examples
#' \dontrun{
#' # reads the website's data by default; optionally a local copy: ptm_set_data_dir("path/to/data")
#' ptm_status("HER2")      # resolves the alias and says why ERBB2 returns what it returns
#' ptm_status(c("EGFR", "TFRC", "NOTAGENE"))
#' }
#' @export
ptm_status <- function(gene) {
  gene <- .ptm_chr(gene, "gene", .PTM_GENE_EXAMPLE)
  one <- function(q) {
    # an ambiguous alias is one row of the answer, never the loss of the other rows
    r <- tryCatch(.ptm_resolve(q), ptm_ambiguous = function(e) list(gene = NA_character_, via = sub("^[^:]*: ", "ambiguous: ", conditionMessage(e))))
    idx <- .ptm_index()
    if (is.na(r$gene)) {
      amb <- startsWith(r$via, "ambiguous"); word <- if (amb) "ambiguous alias" else "not a recognized symbol"
      via <- if (amb) r$via else if (grepl("^.+_[NSTY][0-9]+$", toupper(q))) "a site key, not a gene (ptm_site() lists its partners)" else "not a recognized symbol"
      return(data.frame(query = q, resolved_via = via, gene = NA_character_, status_glyco = word, status_phospho = word, stringsAsFactors = FALSE))
    }
    g <- idx$genes[!is.na(idx$genes$gene) & idx$genes$gene == r$gene, , drop = FALSE]
    if (!nrow(g)) return(data.frame(query = q, resolved_via = r$via, gene = r$gene, status_glyco = "not in the index", status_phospho = "not in the index", stringsAsFactors = FALSE))
    cbind(data.frame(query = q, resolved_via = r$via, stringsAsFactors = FALSE), g)
  }
  parts <- lapply(gene, one)
  cols <- unique(c("query", "resolved_via", "gene", unlist(lapply(parts, names))))
  out <- do.call(rbind, lapply(parts, function(p) { for (nm in setdiff(cols, names(p))) p[[nm]] <- NA; p[, cols, drop = FALSE] }))
  for (nm in intersect(c("status_glyco", "status_phospho"), names(out))) out[[nm]] <- sub("no served pair", "no reported pair", out[[nm]], fixed = TRUE)   # the index's word, in the readers' terms
  rownames(out) <- NULL; class(out) <- c("ptm_status", "data.frame"); out
}
# what was measured for a gene that has no reported pair, in one line
.ptm_no_pair <- function(st) paste0(st$gene, " has no reported pair. Glycosites: ", st$status_glyco, ". Phosphosites: ", st$status_phospho, ".")
# a result of rows of reported pairs: its class, what was asked, and how many rows the query returned (so a printed
# subset can say "3 of 8,579")
.ptm_as_pairs <- function(out, query, by = NULL) {
  class(out) <- c("ptm_pairs", "data.frame"); attr(out, "query") <- query; attr(out, "n_query") <- nrow(out); attr(out, "by") <- by
  out
}

# the pair's total sample count across the cohorts that measured it, from the packed per-cohort
# breakdown; NA when the row carries none. Shown beside every correlation.
.ptm_n_total <- function(packed) vapply(strsplit(ifelse(is.na(packed), "", packed), ";", fixed = TRUE), function(cells) {
  cells <- cells[nzchar(cells)]; if (!length(cells)) return(NA_integer_)
  sum(as.integer(sub(".*:", "", cells))) }, integer(1))
# THE LABELS THE SCORE USED (the site's programsForRow, web/js/data.js): a glycosite's labels are attributed to
# glycoforms ("N3H4F0S1G0=plasma(-0.50)|N2H9F0S0G0=myeloid(+0.51)"), and a row counts as labeled on the glyco
# side only when a label sits on a glycoform among the row's representatives (rep_glycans, "BRCA:N2H8F0S0G0;...").
# Unattributed labels (no "=") fall back to any label, as the site does.
.ptm_glyco_labelled <- function(glyco_programs, rep_glycans) {
  vapply(seq_along(glyco_programs), function(i) {
    lab <- glyco_programs[i]; if (is.na(lab) || !nzchar(lab)) return(FALSE)
    segs <- strsplit(lab, "|", fixed = TRUE)[[1]]; segs <- segs[nzchar(segs)]
    attributed <- any(grepl("=", segs, fixed = TRUE))
    if (!attributed || is.na(rep_glycans[i]) || !nzchar(rep_glycans[i])) return(TRUE)
    used <- sub("^[^:]*:", "", strsplit(rep_glycans[i], ";", fixed = TRUE)[[1]])
    gly <- sub("=.*$", "", segs); body <- sub("^[^=]*=", "", segs)
    any(gly %in% used & nzchar(body))
  }, logical(1))
}
.ptm_filter <- function(rows, f) {
  f <- .ptm_filters(f)
  empty <- c(cis = 0L, k = 0L, composition = 0L, anchored = 0L, significant = 0L, heterogeneous = 0L, loco = 0L, adjusted = 0L, shadow = 0L, self = 0L)
  if (!nrow(rows)) return(list(rows = rows, withheld = empty))
  if ("cohorts_rho" %in% names(rows) && !"n_total" %in% names(rows)) rows$n_total <- .ptm_n_total(rows$cohorts_rho)
  keep <- rep(TRUE, nrow(rows)); withheld <- empty
  if (!f$include_cis) { d <- rows$is_cis %in% TRUE; withheld["cis"] <- sum(keep & d); keep <- keep & !d }
  # the shadow-clean filter: |parent-protein rho| < 0.15 and the pair is not same-protein; NA shadows are withheld
  if (isTRUE(f$shadow_clean_only)) { d <- !(rows$shadow_clean %in% TRUE); withheld["shadow"] <- sum(keep & d); keep <- keep & !d }
  # the self-coupled filter: drop pairs whose stronger side tracks its own
  # protein at |rho| >= 0.6, the candidate screen's cutoff (an abundance story before a site story); a pair with no self-coupling
  # measured is KEPT, since the filter removes a known abundance story rather than requiring a clean one
  if (isTRUE(f$hide_self_coupled)) { d <- !is.na(rows$self_coupling) & rows$self_coupling >= 0.6; withheld["self"] <- sum(keep & d); keep <- keep & !d }
  d <- is.na(rows$k_cohorts) | rows$k_cohorts < f$min_k; withheld["k"] <- sum(keep & d); keep <- keep & !d
  # a composition PAIR carries a program label on BOTH sites; a label on one side is a property of that feature and would hide whole genes
  if (f$hide_composition) { d <- .ptm_glyco_labelled(rows$glyco_programs, rows$rep_glycans) & nzchar(ifelse(is.na(rows$phospho_programs), "", rows$phospho_programs)); withheld["composition"] <- sum(keep & d); keep <- keep & !d }
  if (f$anchored_only) { d <- !(rows$anch_any %in% TRUE); withheld["anchored"] <- sum(keep & d); keep <- keep & !d }
  if (f$significant_only) { d <- !(is.finite(rows$q_value) & rows$q_value < 0.05); withheld["significant"] <- sum(keep & d); keep <- keep & !d }
  # an I2 ceiling: a pooled rho is a common effect only where the cohorts agree; k = 1 rows carry NA and pass
  if (is.finite(f$max_i2)) { d <- is.finite(rows$I2) & rows$I2 > f$max_i2; withheld["heterogeneous"] <- sum(keep & d); keep <- keep & !d }
  # leave-one-cohort-out: a floor on the percentile, among pairs of the same k, of the share of evidence kept
  # without the cohort whose removal costs most; single-cohort pairs carry NA and are withheld, as the site does
  if (is.finite(f$min_loco_pctile)) {
    lp <- if ("loco_pctile_k" %in% names(rows)) rows$loco_pctile_k else rep(NA_real_, nrow(rows))
    d <- !(is.finite(lp) & lp >= f$min_loco_pctile); withheld["loco"] <- sum(keep & d); keep <- keep & !d
  }
  # the adjusted block's filter: keep pairs whose composition-adjusted 95 % interval excludes zero
  if (isTRUE(f$adjusted_excludes_zero)) {
    lo <- if ("rho_adjusted_lo" %in% names(rows)) rows$rho_adjusted_lo else rep(NA_real_, nrow(rows)); hi <- if ("rho_adjusted_hi" %in% names(rows)) rows$rho_adjusted_hi else rep(NA_real_, nrow(rows))
    d <- !(is.finite(lo) & is.finite(hi) & (lo > 0 | hi < 0)); withheld["adjusted"] <- sum(keep & d); keep <- keep & !d }
  out <- rows[keep, , drop = FALSE]; out <- out[order(-out$evidence), , drop = FALSE]
  out <- .ptm_add_ids(out)
  attr(out, "withheld") <- withheld; attr(out, "filters") <- f
  list(rows = out, withheld = withheld)
}
# stack per-query frames that may differ in columns; a `query` column names the query each row came from
.ptm_stack <- function(parts, key = "query", keys = names(parts)) {
  parts <- Map(function(p, k) { if (!is.null(p) && nrow(p)) { p[[key]] <- k; p[, c(key, setdiff(names(p), key)), drop = FALSE] } else NULL }, parts, keys)
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(NULL)
  cols <- unique(unlist(lapply(parts, names)))
  out <- do.call(rbind, lapply(parts, function(p) { for (nm in setdiff(cols, names(p))) p[[nm]] <- NA; p[, cols, drop = FALSE] }))
  rownames(out) <- NULL; out
}

#' Ranked partners of a gene, on one side
#' @param gene one or more symbols, aliases or accessions. With several, a name that is not recognized is said in
#'   a message and skipped; with one, it is an error.
#' @param side "glyco" (the gene's glycosites) or "phospho" (its phosphosites)
#' @param filters the view: [ptm_default_filters()], or a list naming only the settings to change, such as
#'   `list(significant_only = FALSE)`; an unknown name is an error
#' @return a data.frame of rows of reported pairs sorted on total evidence, with `glyco_row_id`, `phospho_row_id`
#'   and `pair_id` (NA where a site has no identifier: its protein has no UniProt accession, or the site does not map onto the canonical sequence);
#'   attribute `withheld` counts rows the view hides (a matrix, one row per gene, when several genes
#'   were asked); attribute `filters` is the view that was applied; attribute `sites` is the gene's site table on
#'   that side (a list by gene for several).
#'   Several genes are stacked with a `query_gene` column, each gene's rows sorted on evidence.
#' @examples
#' \dontrun{
#' # reads the website's data by default; optionally a local copy: ptm_set_data_dir("path/to/data")
#' g <- ptm_gene("EGFR", "glyco")
#' attr(g, "withheld")                                 # what the default view hid, by reason
#' ptm_gene("EGFR", "glyco", modifyList(ptm_default_filters(), list(shadow_clean_only = TRUE)))
#' ptm_gene(c("EGFR", "TFRC"), "glyco")                # stacked, with query_gene
#' }
#' @export
ptm_gene <- function(gene, side = c("glyco", "phospho"), filters = ptm_default_filters()) {
  side <- match.arg(side); gene <- .ptm_chr(gene, "gene", .PTM_GENE_EXAMPLE); filters <- .ptm_filters(filters)
  many <- length(gene) > 1L; seen <- character(0)
  one <- function(q) {
    # one gene asked: an unknown name is an error. Several asked: it is said and skipped, and the rest are returned.
    r <- if (many) tryCatch(.ptm_resolve(q), error = function(e) list(gene = NA_character_, via = sub("^[^:]*: ", "", conditionMessage(e)))) else .ptm_resolve(q)
    if (!is.na(r$gene)) { if (r$gene %in% seen) return(NULL); seen <<- c(seen, r$gene) }   # the same protein named twice: read once
    if (is.na(r$gene)) {
      if (!many) .ptm_not_symbol(q)
      message(q, " skipped: ", if (grepl("^.+_[NSTY][0-9]+$", toupper(q))) "a site key, not a gene (ptm_site() lists its partners)" else if (r$via == "unknown") "not a recognized symbol, alias or UniProt accession" else r$via)
      return(NULL)
    }
    st <- ptm_status(r$gene)
    if ((st$n_served_as_glyco + st$n_served_as_phospho) == 0) { message(.ptm_no_pair(st)); return(NULL) }
    sh <- .ptm_shard(r$gene)
    rows <- sh$rows; rows <- rows[if (side == "glyco") rows$glyco_gene == r$gene else rows$phospho_gene == r$gene, , drop = FALSE]
    res <- .ptm_filter(rows, filters)
    list(gene = r$gene, via = r$via, rows = res$rows, withheld = res$withheld, sites = if (side == "glyco") sh$glyco_sites else sh$phospho_sites)
  }
  parts <- Filter(Negate(is.null), lapply(gene, one))
  if (!length(parts)) return(invisible(NULL))
  sides <- if (side == "glyco") "glycosites" else "phosphosites"
  if (!many) {
    p <- parts[[1]]
    out <- p$rows; attr(out, "withheld") <- p$withheld; attr(out, "sites") <- p$sites; attr(out, "resolved_via") <- p$via
    return(.ptm_as_pairs(out, paste0(p$gene, " ", sides)))
  }
  genes <- vapply(parts, function(p) p$gene, character(1))
  out <- .ptm_stack(lapply(parts, function(p) p$rows), "query_gene", genes)
  if (is.null(out)) { out <- parts[[1]]$rows; out$query_gene <- character(0) }
  attr(out, "withheld") <- do.call(rbind, setNames(lapply(parts, function(p) p$withheld), genes))
  attr(out, "sites") <- setNames(lapply(parts, function(p) p$sites), genes)
  attr(out, "resolved_via") <- setNames(vapply(parts, function(p) p$via, character(1)), genes); attr(out, "filters") <- filters
  .ptm_as_pairs(out, paste0(paste(genes, collapse = ", "), " ", sides), by = "gene")
}

#' The one-hop network of a gene: one edge per partner gene and side of the pair
#'
#' Every reported row in which the gene is on either side, under the filters, rolled up to one edge per
#' (partner gene, role). `role` is "glyco" when the query gene carries the glycosite of the pair and
#' "phospho" when the partner does; a partner reached on both sides gives two rows and must not be
#' deduplicated. An edge's numbers are those of its best row (the largest `evidence`; ties broken on
#' `glyco_site` then `phospho_site` in byte order), never a maximum of one column beside a minimum of
#' another, plus the count of site pairs behind it and their median evidence. Sorted on `evidence_best`
#' descending: the total score sorts, everywhere. Depth is shown, not corrected: a best-of many site pairs
#' is higher than a best-of one, so read `n_site_pairs` beside `evidence_best`. The same table the
#' website's network view draws and exports; the procedure is the one in `web/js/data.js::geneNetwork`.
#' @param gene one gene symbol, alias or UniProt accession
#' @param filters the view: [ptm_default_filters()], or a list naming only the settings to change; `include_cis` admits same-protein rows,
#'   which are counted in the `n_cis_edges` attribute and never become a partner edge
#' @param top_n the website's node cap (default 25): it sets the `drawn` column and nothing else. The frame
#'   always holds every edge of the view. `NULL` marks every row drawn.
#' @param pin partner symbols drawn in addition to the cap, as a pair page does for the partner it came from
#' @return a data.frame with 33 columns, one row per (partner, role) edge, with attributes `withheld`,
#'   `gene`, `resolved_via`, `n_edges`, `n_partner_genes`, `n_cis_edges`, `n_both_roles`, `n_rows` (rows of reported pairs
#'   under the filters), `release`, `top_n` and `note`.
#'   The score is symmetric: an edge has no direction.
#' @examples
#' \dontrun{
#' net <- ptm_network("EGFR")
#' head(net[net$drawn, c("partner", "role", "n_site_pairs", "evidence_best", "rho_at_best")])
#' ptm_pair(net$best_glyco_site[1], net$best_phospho_site[1])
#' }
#' @export
ptm_network <- function(gene, filters = ptm_default_filters(), top_n = 25L, pin = character(0)) {
  if (length(gene) > 1L) stop("ptm_network() takes one gene, got ", length(gene), "; call it once per gene (the tables stack with rbind())", call. = FALSE)
  gene <- .ptm_chr(gene, "gene", .PTM_GENE_EXAMPLE, one = TRUE); filters <- .ptm_filters(filters)
  if (!is.null(top_n) && !(is.numeric(top_n) && length(top_n) == 1L && is.finite(top_n) && top_n >= 0)) stop("top_n is a whole number, or NULL to mark every partner as drawn", call. = FALSE)
  r <- .ptm_resolve(gene)
  if (is.na(r$gene)) .ptm_not_symbol(gene)
  st <- ptm_status(r$gene)
  cols <- c("gene", "partner", "role", "rank", "drawn", "n_site_pairs", "n_own_sites", "n_partner_sites",
            "n_site_pairs_positive", "n_site_pairs_negative", "sign_mixed", "evidence_best", "evidence_median", "pctile_at_best",
            "rho_at_best", "k_at_best", "n_at_best", "I2_at_best", "q_at_best", "q_anchored_at_best",
            "rho_prot", "shadow_clean_at_best", "self_coupling_at_best", "self_glyco_at_best", "self_phospho_at_best",
            "n_site_pairs_anchored", "anch_any", "anch_dose", "anch_gene_editing", "anch_enzymatic",
            "own_site_at_best", "best_glyco_site", "best_phospho_site")
  note <- "an edge aggregates every reported site pair between two proteins on one side of the pair; the statistic is symmetric and licenses no direction, and the number on an edge is the strongest of its site pairs, not a consensus"
  finish <- function(edges, withheld, n_rows, n_cis) {
    attr(edges, "withheld") <- withheld; attr(edges, "filters") <- filters; attr(edges, "gene") <- r$gene; attr(edges, "resolved_via") <- r$via
    attr(edges, "n_edges") <- nrow(edges); attr(edges, "n_partner_genes") <- length(unique(edges$partner))
    attr(edges, "n_cis_edges") <- n_cis; attr(edges, "n_both_roles") <- nrow(edges) - length(unique(edges$partner))
    attr(edges, "n_rows") <- n_rows; attr(edges, "top_n") <- top_n; attr(edges, "note") <- note
    meta <- tryCatch(ptm_release(), error = function(e) list())
    attr(edges, "release") <- if (is.null(meta$release)) NA_character_ else meta$release
    attr(edges, "inference") <- .ptm_inference(meta)
    class(edges) <- c("ptm_network", "data.frame")
    edges
  }
  empty <- as.data.frame(setNames(replicate(length(cols), logical(0), simplify = FALSE), cols), stringsAsFactors = FALSE)
  if ((st$n_served_as_glyco + st$n_served_as_phospho) == 0) { message(.ptm_no_pair(st)); return(invisible(finish(empty, NULL, 0L, 0L))) }
  sh <- .ptm_shard(r$gene)
  edges <- .ptm_network_edges(sh$rows, r$gene, filters, top_n, pin, cols)
  finish(edges$edges, edges$withheld, edges$n_rows, edges$n_cis)
}
# the roll-up itself, on a data.frame of shard rows: kept apart from the fetch so a test can run it on a fixture
.ptm_network_edges <- function(rows0, gene, filters, top_n, pin, cols) {
  res <- .ptm_filter(rows0, filters); rows <- res$rows
  col <- function(nm) if (nm %in% names(rows)) rows[[nm]] else rep(NA, nrow(rows))
  is_cis <- col("is_cis") %in% TRUE
  cis_n <- sum(is_cis); rows <- rows[!is_cis, , drop = FALSE]
  if (!nrow(rows)) return(list(edges = as.data.frame(setNames(replicate(length(cols), logical(0), simplify = FALSE), cols)), withheld = res$withheld, n_rows = nrow(res$rows), n_cis = cis_n))
  as_glyco <- rows$glyco_gene == gene
  role <- ifelse(as_glyco, "glyco", "phospho"); partner <- ifelse(as_glyco, rows$phospho_gene, rows$glyco_gene)
  key <- paste0(partner, "\t", role)
  n_tot <- if ("cohorts_rho" %in% names(rows)) .ptm_n_total(rows$cohorts_rho) else rep(NA_integer_, nrow(rows))
  pick <- function(nm, b) if (nm %in% names(rows)) rows[[nm]][b] else rep(NA, length(b))
  # the best row of each edge: the largest evidence, ties on glyco_site then phospho_site in byte order (radix = C locale)
  o <- order(key, -rows$evidence, rows$glyco_site, rows$phospho_site, method = "radix")
  b <- o[!duplicated(key[o])]
  g <- match(key, key[b]); nb <- length(b)
  count <- function(hit) tabulate(g[which(hit)], nb)
  own_site <- ifelse(as_glyco, rows$glyco_site, rows$phospho_site); other_site <- ifelse(as_glyco, rows$phospho_site, rows$glyco_site)
  pos <- count(rows$rho_combined >= 0); neg <- count(rows$rho_combined < 0)
  flag <- function(nm) if (nm %in% names(rows)) rows[[nm]] %in% TRUE else rep(FALSE, nrow(rows))
  edges <- data.frame(gene = rep(gene, nb), partner = partner[b], role = role[b], rank = 0L, drawn = FALSE,
                      n_site_pairs = tabulate(g, nb), n_own_sites = count(!duplicated(paste(g, own_site))), n_partner_sites = count(!duplicated(paste(g, other_site))),
                      n_site_pairs_positive = pos, n_site_pairs_negative = neg, sign_mixed = pos > 0 & neg > 0,
                      evidence_best = rows$evidence[b], evidence_median = vapply(split(rows$evidence, factor(g, levels = seq_len(nb))), stats::median, numeric(1), USE.NAMES = FALSE),
                      pctile_at_best = as.numeric(pick("evidence_pctile_k", b)),
                      rho_at_best = rows$rho_combined[b], k_at_best = as.integer(rows$k_cohorts[b]), n_at_best = as.integer(n_tot[b]), I2_at_best = as.numeric(pick("I2", b)),
                      q_at_best = as.numeric(pick("q_value", b)), q_anchored_at_best = as.numeric(pick("q_value_anchored", b)),
                      rho_prot = as.numeric(pick("rho_prot", b)), shadow_clean_at_best = as.logical(pick("shadow_clean", b)), self_coupling_at_best = as.numeric(pick("self_coupling", b)),
                      self_glyco_at_best = as.numeric(pick("self_glyco", b)), self_phospho_at_best = as.numeric(pick("self_phospho", b)),
                      n_site_pairs_anchored = count(flag("anch_any")), anch_any = count(flag("anch_any")) > 0, anch_dose = count(flag("anch_dose")) > 0,
                      anch_gene_editing = count(flag("anch_gene_editing")) > 0, anch_enzymatic = count(flag("anch_enzymatic")) > 0,
                      own_site_at_best = own_site[b], best_glyco_site = rows$glyco_site[b], best_phospho_site = rows$phospho_site[b], stringsAsFactors = FALSE)
  edges <- edges[order(-edges$evidence_best, edges$partner, ifelse(edges$role == "glyco", 0L, 1L), method = "radix"), cols, drop = FALSE]
  rownames(edges) <- NULL; edges$rank <- seq_len(nrow(edges))
  edges$drawn <- if (is.null(top_n)) rep(TRUE, nrow(edges)) else (edges$rank <= top_n | edges$partner %in% pin)
  list(edges = edges, withheld = res$withheld, n_rows = nrow(res$rows), n_cis = cis_n)
}

.ptm_site_key <- function(site) {
  site <- .ptm_canon_site(site)
  m <- regmatches(site, regexec("^(.+)_([NSTY])([0-9]+)$", site))[[1]]
  if (!length(m)) stop(site, " is not a site key (GENE_N123 or GENE_S123)", call. = FALSE)
  list(gene = m[2], key = paste0(m[2], "_", m[3], m[4]), side = if (m[3] == "N") "glyco" else "phospho")
}

#' Partners of one or more sites
#' @param site site keys such as EGFR_N603 or EGFR_T693. A site in no reported pair is an error; with several
#'   sites it is said in a message and skipped.
#' @param filters the view: [ptm_default_filters()], or a list naming only the settings to change
#' @return a data.frame of rows of reported pairs, sorted on total evidence, with identifier columns and a `withheld`
#'   attribute; several sites are stacked with a `query_site` column, each site's rows sorted on evidence
#' @examples
#' \dontrun{
#' # reads the website's data by default; optionally a local copy: ptm_set_data_dir("path/to/data")
#' ptm_site("EGFR_T693")   # glycosites that co-vary with one phosphosite
#' }
#' @export
ptm_site <- function(site, filters = ptm_default_filters()) {
  site <- .ptm_chr(site, "site", .PTM_SITE_EXAMPLE); filters <- .ptm_filters(filters)
  many <- length(site) > 1L; seen <- character(0)
  one <- function(s) {
    k <- .ptm_site_key(s)
    if (k$key %in% seen) return(NULL)   # the same site named twice: read once
    seen <<- c(seen, k$key)
    sh <- .ptm_or_absent(.ptm_shard(k$gene))
    rows <- if (is.null(sh)) NULL else sh$rows[sh$rows$glyco_site == k$key | sh$rows$phospho_site == k$key, , drop = FALSE]
    if (is.null(rows) || !nrow(rows)) {
      # a site with no reported pair is said as such, never shown as "no pair passes the filters"
      why <- paste0(k$key, " is in no reported pair of this release; ", if (is.null(sh)) paste0("ptm_status(\"", k$gene, "\") says what was measured") else paste0("ptm_sites(\"", k$gene, "\") lists the reported sites of ", k$gene))
      if (!many) stop(why, call. = FALSE)
      message(why); return(NULL)
    }
    c(.ptm_filter(rows, filters), list(key = k$key))
  }
  parts <- Filter(Negate(is.null), lapply(site, one))
  if (!length(parts)) return(invisible(NULL))
  keys <- vapply(parts, function(p) p$key, character(1))
  if (!many) {
    out <- parts[[1]]$rows; attr(out, "withheld") <- parts[[1]]$withheld
    return(.ptm_as_pairs(out, paste("partners of", keys)))
  }
  out <- .ptm_stack(lapply(parts, function(p) p$rows), "query_site", keys)
  if (is.null(out)) { out <- parts[[1]]$rows; out$query_site <- character(0) }
  attr(out, "withheld") <- do.call(rbind, setNames(lapply(parts, function(p) p$withheld), keys)); attr(out, "filters") <- filters
  .ptm_as_pairs(out, paste("partners of", paste(keys, collapse = ", ")), by = "site")
}

# ONE WORDING for a pair the release does not report, in R, Python and the MCP server
.PTM_ABSENT_PAIR <- "not among the reported pairs of this release; the complete scored pairs are released with the article"
# A pair named by its stable identifier, "PTMG:<accession>:N<position>|PTMP:<accession>:<residue><position>", as the
# one argument; split into its two site identifiers (read as sites by .ptm_pair_keys).
.ptm_split_pair_id <- function(x, fn) {
  x <- toupper(trimws(as.character(x)))
  ok <- grepl("^PTMG:[^|]+[|]PTMP:[^|]+$", x)
  if (!length(x) || !all(ok)) stop(fn, "() takes a glycosite and a phosphosite (", fn, "(\"IL6ST_N131\", \"TMX1_S270\")) or one pair_id (",
                                   fn, "(\"PTMG:P40189:N131|PTMP:Q9H3N1:S270\")); got ", if (length(x)) x[!ok][1] else "nothing", call. = FALSE)
  parts <- strsplit(x, "|", fixed = TRUE)
  list(g = vapply(parts, `[`, "", 1), p = vapply(parts, `[`, "", 2))
}
# site key -> its side, from the key's residue
.ptm_key_side <- function(k) ifelse(grepl("_N[0-9]+$", k), "glyco", ifelse(grepl("_[STY][0-9]+$", k), "phospho", NA_character_))
# The site keys a pair argument stands for: a site key as given, or every key a site identifier (row_id) maps to in
# the shipped site list (a phosphosite identifier can stand for two CPTAC labels of one canonical residue). A
# phosphosite given first is read in the right order, with a message.
.ptm_pair_keys <- function(gs, ps) {
  keys <- function(x) {
    if (!grepl("^PTM[GP]:", x)) return(x)
    .ptm_ids_from_package(); .ptm_ids_from_release()
    k <- names(.ptm_env$ids)[!is.na(.ptm_env$ids) & .ptm_env$ids == x]
    if (!length(k)) stop(x, " is not a site identifier of this release", call. = FALSE)
    sort(unique(k))
  }
  g <- keys(gs); p <- keys(ps)
  sg <- unique(.ptm_key_side(g)); sp <- unique(.ptm_key_side(p))
  if (anyNA(c(sg, sp))) stop(c(g, p)[is.na(.ptm_key_side(c(g, p)))][1], " is not a site key (GENE_N123 or GENE_S123)", call. = FALSE)
  if (identical(sg, "phospho") && identical(sp, "glyco")) {
    message("the glycosite comes first: reading the pair as ", p[1], " x ", g[1])
    tmp <- g; g <- p; p <- tmp
  } else if (!identical(sg, "glyco") || !identical(sp, "phospho")) {
    stop("a pair is one glycosite (GENE_N123) and one phosphosite (GENE_S123, GENE_T123 or GENE_Y123); got ", g[1], " and ", p[1], call. = FALSE)
  }
  list(g = g, p = p)
}

#' One or more pairs, each with its per-cohort breakdown
#' @param glyco_site,phospho_site the two site keys (any case); vectors of equal length (or one recycled) ask several
#'   pairs. A pair's stable identifier, `"PTMG:<accession>:N<position>|PTMP:<accession>:<residue><position>"`, may be
#'   given as `glyco_site` alone, and a site's identifier may stand in for its key. Given phosphosite first, the
#'   pair is read in the right order, with a message.
#' @return for one pair a list: `row` (the reported row with its identifiers, or NULL), `cohorts` (cohort, its display label, rho, n, a
#'   95 % interval `lo`, `hi` at the cohort's weight n - 3 - k_cov, and `k_cov`, its covariate count from the cohort dictionary), `n_total` (samples across them), `rep_glycans` (the glycoform each cohort's correlation was measured
#'   on), `glycoforms` (every glycoform of the glycosite scored against the phosphosite in the glycoform layer, with that
#'   layer's own q; NULL when the data source does not carry the tier), and `note` (the reading rule: the score is
#'   symmetric and licenses no direction). For several pairs the same list with every frame stacked and keyed by
#'   `glyco_site` and `phospho_site`, `n_total` a named vector, and `missing` the pairs not reported.
#' @examples
#' \dontrun{
#' # reads the website's data by default; optionally a local copy: ptm_set_data_dir("path/to/data")
#' p <- ptm_pair("TFRC_N727", "COIL_T303")
#' p$cohorts; p$glycoforms
#' }
#' @export
ptm_pair <- function(glyco_site, phospho_site = NULL) {
  glyco_site <- .ptm_chr(glyco_site, "glyco_site", "a glycosite key such as \"IL6ST_N131\", or a pair_id")
  if (!is.null(phospho_site)) phospho_site <- .ptm_chr(phospho_site, "phospho_site", "a phosphosite key such as \"TMX1_S270\"")
  if (is.null(phospho_site)) { k <- .ptm_split_pair_id(glyco_site, "ptm_pair"); glyco_site <- k$g; phospho_site <- k$p }
  if (length(glyco_site) != length(phospho_site) && length(glyco_site) != 1L && length(phospho_site) != 1L)
    stop("glyco_site and phospho_site pair up element by element: give the same number of each, or one of either (got ", length(glyco_site), " and ", length(phospho_site), ")", call. = FALSE)
  n <- max(length(glyco_site), length(phospho_site)); glyco_site <- rep_len(.ptm_canon_site(glyco_site), n); phospho_site <- rep_len(.ptm_canon_site(phospho_site), n)
  absent <- list(row = NULL, cohorts = NULL, n_total = NA_integer_, rep_glycans = NULL, glycoforms = NULL, note = .PTM_ABSENT_PAIR)
  one <- function(gs, ps) {
    k <- .ptm_pair_keys(gs, ps); gs <- k$g[1]; ps <- k$p[1]
    gene <- sub("_N[0-9]+$", "", gs)
    sh <- .ptm_or_absent(.ptm_shard(gene))
    if (is.null(sh)) return(absent)
    # the partner's data file teaches the partner site's identifier, so the pair id is complete
    pgene <- .ptm_site_key(ps)$gene
    if (pgene != gene) .ptm_or_absent(.ptm_shard(pgene))
    row <- sh$rows[sh$rows$glyco_site %in% k$g & sh$rows$phospho_site %in% k$p, , drop = FALSE]
    if (!nrow(row)) return(absent)
    row <- row[order(-row$evidence), , drop = FALSE][1, , drop = FALSE]; gs <- row$glyco_site; ps <- row$phospho_site
    row <- .ptm_add_ids(row); if (!"n_total" %in% names(row)) row$n_total <- .ptm_n_total(row$cohorts_rho)
    cells <- strsplit(row$cohorts_rho[1], ";", fixed = TRUE)[[1]]
    parts <- do.call(rbind, strsplit(cells, ":", fixed = TRUE))
    # the 95 % interval at each cohort's own weight n - 3 - k, k its covariate count (cohorts.json k_cov; 1 when unstated)
    kc <- .ptm_kcov(); k <- ifelse(parts[, 1] %in% names(kc), unname(kc[parts[, 1]]), 1)
    rho <- as.numeric(parts[, 2]); nn <- as.integer(parts[, 3]); z <- atanh(pmin(pmax(rho, -0.999999), 0.999999))
    w <- nn - 3 - k; se <- ifelse(w > 0, 1 / sqrt(pmax(w, 1e-12)), NA_real_)
    # the interval to four decimals, the precision rho is released at
    co <- data.frame(cohort = parts[, 1], label = .ptm_cohort_label(parts[, 1]), rho = rho, n = nn, lo = round(tanh(z - 1.959964 * se), 4), hi = round(tanh(z + 1.959964 * se), 4), k_cov = k, stringsAsFactors = FALSE)
    reps <- if (!is.na(row$rep_glycans[1]) && nzchar(row$rep_glycans[1])) { p <- do.call(rbind, strsplit(strsplit(row$rep_glycans[1], ";", fixed = TRUE)[[1]], ":", fixed = TRUE)); data.frame(cohort = p[, 1], rep_glycan = p[, 2], stringsAsFactors = FALSE) } else NULL
    # the glycoform tier is optional in a data source: absent file or 404 -> NULL, quietly
    gf <- .ptm_or_absent(.ptm_json(paste0("glycoforms/", gs, ".json.gz")))
    glycoforms <- NULL
    if (!is.null(gf) && !is.null(gf$pairs[[ps]])) {
      m <- gf$pairs[[ps]]
      if (is.list(m)) m <- do.call(rbind, lapply(m, function(a) as.data.frame(as.list(unlist(a)), stringsAsFactors = FALSE)))
      m <- as.data.frame(m, stringsAsFactors = FALSE); names(m) <- gf$columns
      glycoforms <- data.frame(glycan = as.character(m$glycan), rho = as.numeric(m$rho), evidence = as.numeric(m$evidence), k = as.integer(m$k),
                               I2 = as.numeric(m$I2), q = as.numeric(m$q), anchored = as.logical(m$anch), stringsAsFactors = FALSE)
      for (nm in intersect(c("sum_weight", "q_value_anchored"), names(m))) glycoforms[[nm]] <- as.numeric(m[[nm]])
      attr(glycoforms, "inference") <- .ptm_inference(gf)
    }
    list(inference = sh$inference, row = row, cohorts = co, n_total = sum(co$n), rep_glycans = reps, glycoforms = glycoforms, note = "the score is symmetric and licenses no direction")
  }
  parts <- Map(function(g, p) structure(one(g, p), class = "ptm_pair"), glyco_site, phospho_site)
  # the keys as reported (a pair_id or a swapped order is read as the two site keys)
  glyco_site <- vapply(seq_len(n), function(i) if (!is.null(parts[[i]]$row)) parts[[i]]$row$glyco_site[1] else glyco_site[i], character(1))
  phospho_site <- vapply(seq_len(n), function(i) if (!is.null(parts[[i]]$row)) parts[[i]]$row$phospho_site[1] else phospho_site[i], character(1))
  if (n == 1) return(parts[[1]])
  keyed <- function(field) { fr <- lapply(seq_len(n), function(i) { x <- parts[[i]][[field]]; if (is.null(x) || !nrow(x)) return(NULL); cbind(data.frame(glyco_site = glyco_site[i], phospho_site = phospho_site[i], stringsAsFactors = FALSE), x) }); fr <- Filter(Negate(is.null), fr); if (length(fr)) { out <- do.call(rbind, fr); rownames(out) <- NULL; out } else NULL }
  found <- vapply(parts, function(p) !is.null(p$row), logical(1))
  rows <- if (any(found)) { r <- do.call(rbind, lapply(parts[found], function(p) p$row)); rownames(r) <- NULL; r } else NULL
  structure(list(row = rows, cohorts = keyed("cohorts"), n_total = setNames(vapply(parts, function(p) as.integer(p$n_total), integer(1)), paste(glyco_site, phospho_site, sep = "|")),
       rep_glycans = keyed("rep_glycans"), glycoforms = keyed("glycoforms"),
       missing = data.frame(glyco_site = glyco_site[!found], phospho_site = phospho_site[!found], stringsAsFactors = FALSE),
       note = "the score is symmetric and licenses no direction; a pair absent from `row` is not among the reported pairs of this release"), class = "ptm_pair")
}
