# The machine-readable reference, the column list and the citation.
#
# inst/reference.json is the one file an agent framework loads: every exported function with its
# arguments and return columns, the reading rules with their thresholds, and the columns of a reported row
# with their meaning. ptm_explain() reads its rules from it; a test holds it to the NAMESPACE and to
# the functions' formals, so it cannot describe a function that does not exist or an argument that
# was renamed.

#' The package reference as data
#' @param as "list" (the default) or "json" (the file's text, for a tool-description loader)
#' @return the parsed reference: `package`, `version`, `purpose`, `release_contract`, `rules`, `functions`, `columns`,
#'   `filters` and `tool_conventions`
#' @examples
#' ref <- ptm_reference()
#' vapply(ref$functions, function(f) f$name, character(1))
#' @export
ptm_reference <- function(as = c("list", "json")) {
  as <- match.arg(as)
  if (as == "json") { p <- system.file("reference.json", package = "ptmlink"); if (!nzchar(p)) p <- file.path("inst", "reference.json"); return(paste(readLines(p, warn = FALSE), collapse = "\n")) }
  structure(.ptm_reference_raw(), class = "ptm_reference")
}

#' The columns of a reported row, with their meaning
#' @return a data.frame `column`, `type`, `meaning`, `in_shards` (carried by the per-gene files the website
#'   reads; the rest are in the released Parquet only). When a release directory is set and arrow is installed,
#'   the type is read from the Parquet schema.
#' @export
ptm_columns <- function() {
  cols <- .ptm_reference_raw()$columns
  out <- data.frame(column = vapply(cols, function(c) c$name, character(1)), type = vapply(cols, function(c) if (is.null(c$type)) NA_character_ else c$type, character(1)),
                    meaning = vapply(cols, function(c) c$meaning, character(1)), in_shards = vapply(cols, function(c) isTRUE(c$in_shards), logical(1)), stringsAsFactors = FALSE)
  d <- .ptm_release_dir()
  if (nzchar(d) && file.exists(file.path(d, "site_served.parquet")) && requireNamespace("arrow", quietly = TRUE)) {
    s <- arrow::open_dataset(file.path(d, "site_served.parquet"))$schema
    ty <- setNames(vapply(s$fields, function(f) f$type$ToString(), character(1)), names(s))
    extra <- setdiff(names(ty), out$column)
    if (length(extra)) out <- rbind(out, data.frame(column = extra, type = unname(ty[extra]), meaning = "see the release data dictionary", in_shards = FALSE, stringsAsFactors = FALSE))
    hit <- out$column %in% names(ty); out$type[hit] <- unname(ty[out$column[hit]])
  }
  class(out) <- c("ptm_columns", "data.frame")
  out
}

# a small CFF reader: top-level scalars, the authors list and the preferred citation's journal, year and DOI; enough
# for a citation, not a YAML parser
.ptm_read_cff <- function(path) {
  l <- readLines(path, warn = FALSE)
  scalar <- function(key) { i <- grep(paste0("^", key, ":"), l); if (!length(i)) return(NA_character_); v <- sub(paste0("^", key, ":\\s*"), "", l[i[1]]); v <- gsub('^"|"$', "", v); if (v == ">-" || v == "") { j <- i[1] + 1; buf <- character(0); while (j <= length(l) && grepl("^\\s+", l[j])) { buf <- c(buf, trimws(l[j])); j <- j + 1 }; v <- paste(buf, collapse = " ") }; v }
  a0 <- grep("^authors:", l); a1 <- grep("^[a-z-]+:", l); a1 <- a1[a1 > a0][1]; if (is.na(a1)) a1 <- length(l) + 1
  block <- if (length(a0)) l[(a0 + 1):(a1 - 1)] else character(0)
  fam <- trimws(sub("^\\s*-?\\s*family-names:\\s*", "", grep("family-names:", block, value = TRUE))); giv <- trimws(sub("^\\s*given-names:\\s*", "", grep("given-names:", block, value = TRUE)))
  authors <- if (length(fam)) paste(giv[seq_along(fam)], fam) else character(0)
  lic <- { i <- grep("^license:", l); if (!length(i)) NA_character_ else { v <- sub("^license:\\s*", "", l[i[1]]); if (nzchar(v)) v else { j <- i[1] + 1; buf <- character(0); while (j <= length(l) && grepl("^\\s+-", l[j])) { buf <- c(buf, trimws(sub("^\\s+-\\s*", "", l[j]))); j <- j + 1 }; paste(buf, collapse = ", ") } } }
  # the paper that describes the resource (its preprint, then its article), when the file names one
  p0 <- grep("^preferred-citation:", l)
  paper <- if (length(p0)) {
    nx <- grep("^[a-z-]+:", l); nx <- nx[nx > p0[1]][1]; if (is.na(nx)) nx <- length(l) + 1
    b <- l[(p0[1] + 1):(nx - 1)]
    field <- function(key) { v <- grep(paste0("^  ", key, ":"), b, value = TRUE); if (!length(v)) NA_character_ else gsub('^"|"$', "", trimws(sub(paste0("^  ", key, ":\\s*"), "", v[1]))) }
    list(journal = field("journal"), year = field("year"), doi = field("doi"), url = field("url"), status = field("status"))
  } else NULL
  list(title = scalar("title"), version = scalar("version"), date = scalar("date-released"), authors = authors, license = lic, doi = scalar("doi"), url = scalar("url"), paper = paper, message = scalar("message"))
}

# The citation line. With a paper: its title, journal, year and DOI, then the resource's release, date and address, in
# the words of the website's About page.
# Without one (a release written before the paper existed): the resource alone. inst/CITATION and the Python
# reader's cite() print the same line.
.ptm_cite_text <- function(c) {
  has <- function(v) !is.null(v) && !is.na(v) && nzchar(v)
  resource <- paste0("ersion ", c$version, if (has(c$date)) paste0(", released ", c$date) else "", if (has(c$doi)) paste0(", doi:", c$doi) else "")
  p <- c$paper
  if (!is.null(p) && has(p$doi))
    return(paste0(c$title, ". ", if (has(p$journal)) p$journal else "Preprint", if (has(p$year)) paste0(", ", p$year) else "", ". doi:", p$doi,
                  ". PTMlink release v", c$version, if (has(c$date)) paste0(", released ", c$date) else "", if (has(c$doi)) paste0(", doi:", c$doi) else "",
                  if (has(c$url)) paste0(", ", c$url) else "", ". Data CC BY 4.0, code MIT."))
  paste0(c$title, if (has(c$date)) paste0(" (", substr(c$date, 1, 4), ")") else "", ". V", resource, if (has(c$url)) paste0(". ", c$url) else "", ". Data CC BY 4.0, code MIT.")
}

#' How to cite the release in use
#' @return a list: `title`, `version`, `date`, `license`, `doi` (the data archive's, NA until the deposit exists) and
#'   `url`, `paper` (the paper that describes the resource, when the citation file names one: `journal`, `year`,
#'   `doi`, `url` and `status`, which is `"preprint"` until the article is published), `release` (the stamp of the
#'   data source in use, if one is set), `text` (one citation line: the title, the paper's journal, year and DOI,
#'   then the resource's version and address; it carries no author list, which the paper itself gives) and
#'   `r_citation` (the package's own citation object, for `print` or `toBibtex`). Read from the release directory's `CITATION.cff` when one is set, else from the copy
#'   the package ships.
#' @examples
#' ptm_cite()$text
#' @export
ptm_cite <- function() {
  d <- .ptm_release_dir(); p <- if (nzchar(d) && file.exists(file.path(d, "CITATION.cff"))) file.path(d, "CITATION.cff") else system.file("extdata", "ptmlink_citation.cff", package = "ptmlink")
  if (!nzchar(p)) p <- file.path("inst", "extdata", "ptmlink_citation.cff")
  c <- .ptm_read_cff(p); c$authors <- NULL; c$message <- NULL
  rel <- tryCatch(ptm_release()$release, error = function(e) NA_character_)
  c$release <- rel
  c$text <- .ptm_cite_text(c)
  c$r_citation <- tryCatch(utils::citation("ptmlink"), error = function(e) NULL)
  class(c) <- "ptm_cite"
  c
}
