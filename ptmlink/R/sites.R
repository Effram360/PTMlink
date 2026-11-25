# ptm_sites -- the reported sites of a protein with their canonical position and, where the release
# carries the structure sidecar, the modeled position and confidence of that residue.
#
# What other tools need from PTMlink is the mapping, not the picture: which residues of which
# accession this database serves, where they sit in the AlphaFold model of the canonical sequence,
# and how confident the model is there. A reader can hand that table to AlphaFold 3, a PTM predictor
# or a structure viewer without rebuilding the mapping. Coordinates are the alpha carbon of the
# residue in the AlphaFold Protein Structure Database model (AlphaFold 2, CC BY 4.0), centered as the
# sidecar ships them; a distance between two of them is a fact about the model, never about coupling.

#' The reported sites of one or more proteins, with their modeled position
#' @param x gene symbols, aliases or UniProt accessions
#' @param structure also read the gene's AlphaFold record, found as the website finds it (`structure/bundles.json`, or
#'   `structure/<GENE>.json.gz` in a data directory without bundles)
#' @return a data.frame, one row per reported site, glycosites first and each side in order along the protein: `query`, `gene`, `uniprot_acc`, `side`
#'   ("glyco" or "phospho"), `site`, `row_id`, `residue`, `position_canonical`, `residue_ok` (the residue at that
#'   position of the canonical sequence is the one the site claims), `n_cohorts_served`, `self_coupling_median`,
#'   `window_15` (the 15-residue window centered on the site), `frame_offset` (phosphosites: the CPTAC label's
#'   position minus the canonical one), `topology_class` (glycosites), `glyco_site_ambiguous`; and, with the
#'   structure sidecar, `model` (the AlphaFold entry), `model_coverage` (the modeled share of the chain),
#'   `modelled` (the position lies inside the modeled stretch), `model_residue`, `residue_matches_model`, `plddt`,
#'   `x`, `y`, `z`; and `structure_status`, which says why the structure columns are filled or NA ("model read", "no
#'   AlphaFold model in this data source", "not requested", or "could not be read: <reason>", with a message).
#' @examples
#' \dontrun{
#' # reads the website's data by default; optionally a local copy: ptm_set_data_dir("path/to/data")
#' s <- ptm_sites("EGFR")
#' s[s$side == "glyco", c("site", "position_canonical", "plddt", "residue_matches_model")]
#' }
#' @export
ptm_sites <- function(x, structure = TRUE) {
  x <- .ptm_chr(x, "x", .PTM_GENE_EXAMPLE)
  if (!isTRUE(structure) && !isFALSE(structure)) stop("structure is TRUE or FALSE", call. = FALSE)
  one <- function(q) {
    # one name asked: an ambiguous alias is an error; several asked: it is said and skipped
    r <- if (length(x) > 1L) tryCatch(.ptm_resolve(q), ptm_ambiguous = function(e) { message(q, " skipped: ", sub("^[^:]*: ", "", conditionMessage(e))); NULL }) else .ptm_resolve(q)
    if (is.null(r)) return(NULL)
    if (is.na(r$gene)) { message(q, if (grepl("^.+_[NSTY][0-9]+$", toupper(q))) " is a site key, not a gene; ptm_site() lists its partners" else " is not a recognized symbol, alias or UniProt accession"); return(NULL) }
    sh <- .ptm_or_absent(.ptm_shard(r$gene))
    if (is.null(sh)) { message(r$gene, " has no reported pair, so it has no reported site; ptm_status(\"", r$gene, "\") says what was measured"); return(NULL) }
    pick <- function(tab, nm, default = NA) if (is.data.frame(tab) && nm %in% names(tab)) tab[[nm]] else rep(default, if (is.data.frame(tab)) nrow(tab) else 0)
    g <- sh$glyco_sites; p <- sh$phospho_sites
    gf <- if (is.data.frame(g) && nrow(g)) data.frame(query = q, gene = r$gene, uniprot_acc = pick(g, "uniprot_acc"), side = "glyco", site = g$site, row_id = pick(g, "row_id"),
      residue = "N", position_canonical = as.integer(pick(g, "position")), residue_ok = as.logical(pick(g, "residue_ok")), n_cohorts_served = as.integer(pick(g, "n_cohorts_served")),
      self_coupling_median = as.numeric(pick(g, "self_coupling_median")), window_15 = pick(g, "window_15"), frame_offset = NA_integer_, topology_class = pick(g, "topology_class"),
      glyco_site_ambiguous = as.logical(pick(g, "glyco_site_ambiguous")), stringsAsFactors = FALSE) else NULL
    pf <- if (is.data.frame(p) && nrow(p)) data.frame(query = q, gene = r$gene, uniprot_acc = pick(p, "uniprot_acc"), side = "phospho", site = p$site, row_id = pick(p, "row_id"),
      residue = ifelse(is.na(pick(p, "residue_canonical")), pick(p, "residue_shipped"), pick(p, "residue_canonical")), position_canonical = as.integer(pick(p, "position_canonical")), residue_ok = as.logical(pick(p, "residue_ok")),
      n_cohorts_served = as.integer(pick(p, "n_cohorts_served")), self_coupling_median = as.numeric(pick(p, "self_coupling_median")), window_15 = pick(p, "window_15"),
      frame_offset = as.integer(pick(p, "position_gencode34")) - as.integer(pick(p, "position_canonical")), topology_class = NA_character_, glyco_site_ambiguous = NA, stringsAsFactors = FALSE) else NULL
    out <- rbind(gf, pf); if (is.null(out)) return(NULL)
    # along the protein: glycosites, then phosphosites, each by position (the files list S1026 before S991)
    out <- out[order(out$side != "glyco", out$position_canonical, out$site, method = "radix"), , drop = FALSE]
    # the model: found as the website finds it; a failure to read it is said, never hidden, and never fails the call
    st <- NULL; status <- "not requested"
    if (isTRUE(structure)) {
      st <- tryCatch(.ptm_structure(r$gene), error = function(e) e)
      if (inherits(st, "error")) {
        why <- conditionMessage(st); st <- NULL
        status <- paste("could not be read:", why)
        message("the structure model of ", r$gene, " could not be read, so its structure columns are NA: ", why)
      } else if (is.null(st)) {
        status <- "no AlphaFold model in this data source"
        message(r$gene, ": no AlphaFold model in this data source, so its structure columns are NA")
      } else status <- "model read"
    }
    if (is.null(st)) { out$model <- NA_character_; out$model_coverage <- NA_real_; out$modelled <- NA; out$model_residue <- NA_character_; out$residue_matches_model <- NA; out$plddt <- NA_real_; out$x <- NA_real_; out$y <- NA_real_; out$z <- NA_real_; out$structure_status <- status; return(out) }
    first <- as.integer(st$first); last <- as.integer(st$last); seq_len_ <- as.integer(st$seq_len); n <- as.integer(st$n)
    xyz <- if (is.matrix(st$xyz)) st$xyz else do.call(rbind, st$xyz); plddt <- as.numeric(st$plddt); seq <- as.character(st$seq)
    pos <- out$position_canonical; inside <- !is.na(pos) & pos >= first & pos <= last
    idx <- ifelse(inside, pos - first + 1L, NA_integer_)
    out$model <- if (!is.null(st$entity)) paste0(st$entity, " v", st$v) else NA_character_
    out$model_coverage <- if (length(seq_len_) && seq_len_ > 0) n / seq_len_ else NA_real_
    out$modelled <- inside
    out$model_residue <- ifelse(inside, substring(seq, idx, idx), NA_character_)
    out$residue_matches_model <- ifelse(inside, out$model_residue == out$residue, NA)
    out$plddt <- ifelse(inside, plddt[idx], NA_real_)
    out$x <- ifelse(inside, xyz[idx, 1], NA_real_); out$y <- ifelse(inside, xyz[idx, 2], NA_real_); out$z <- ifelse(inside, xyz[idx, 3], NA_real_)
    out$structure_status <- status
    out
  }
  parts <- Filter(Negate(is.null), lapply(x, one))
  if (!length(parts)) return(invisible(NULL))
  out <- do.call(rbind, parts); rownames(out) <- NULL
  class(out) <- c("ptm_sites", "data.frame")
  attr(out, "note") <- "coordinates are the alpha carbon in the AlphaFold Protein Structure Database model of the canonical sequence (AlphaFold 2, CC BY 4.0), centered; pLDDT is the model's own confidence; a distance between two sites is a fact about the model, never about coupling"
  out
}
