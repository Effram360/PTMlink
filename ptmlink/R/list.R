# ptm_list -- a list of proteins or sites against the reported pairs.
#
# The common case: a user has a list from a differential analysis and asks what the database says
# about it as a whole, and above all which sites on the other side of the pairs the list converges on.
# The answer has these parts, each defined here once and named the same in the result, the reference,
# the Python reader and the website's List page:
#
#   level          "protein" or "site". Auto uses protein level for gene-only lists. If any entry
#                  names a site, site entries stay literal and gene entries expand to their reported
#                  sites. Explicit protein level retains protein expansion for every entry.
#   member         a unit of the list on the tier: a protein at protein level, a site at site level.
#   internal pair  a reported pair with a member on both sides and two different proteins. Same-protein
#                  pairs are set apart: they belong to every list that contains the protein.
#   neighbourhood  a reported pair with a member on exactly one side; the site on the other side is a
#                  partner site and its protein a partner.
#   partner site   the ranked answer to "which sites does my list correlate with": a site outside the
#                  list reached by more members than its reach in the tier predicts. Of the N units on
#                  the members' side of the tier (proteins at protein level, sites at site level) the site
#                  pairs with K; the list holds n of those N and the site reaches x of the n. The
#                  hypergeometric P of x, BH-adjusted over every partner site tested in one family; the
#                  fold is (x/n) / (K/N). Beside the test: the strongest pair in the view.
#   convergent     the same test for a partner protein: of the N member-side units it pairs with K
#   partner        through any of its sites, and it reaches x of the list's n. One family over every
#                  partner protein; a separate family from the sites.
#   connectivity   the number of internal pairs against the number expected from the members' own
#                  numbers of reported pairs on each side (sum over member pairs on different proteins of
#                  k_glyco[i] * k_phospho[j] / M, the Chung-Lu expectation), with a permutation P from
#                  random lists of the same size drawn to match the members' degrees.
#
# The row tables obey the view (filters) like every other function. The tests count the SIGNIFICANT reported
# pairs (q_value < 0.05) between two different proteins, and nothing else: a pair that is reported only because a
# perturbation study measured it, with no significant co-variation in tumors, is not a pair the list "reaches".
# That universe is fixed: the view chooses which rows are shown, not what the tests count, so a test does not
# move when a display filter is toggled. Nothing here recomputes a score; the tests are counts over released
# rows, and a partner reached by one member is a pair-level finding that the neighborhood table, sorted on
# evidence, already carries.

#' Test a list of proteins or sites against the reported pairs
#'
#' Takes a list from a differential analysis and answers which sites and proteins on the other side of the pairs
#' the list co-varies with more often than expected, and whether its members co-vary with each other more often
#' than expected. The three tests count the significant reported pairs (q < 0.05) between two different
#' proteins; a pair reported only for its perturbation support is not counted. The row tables follow `filters`.
#'
#' For the v1.0.0 release a call takes about 30 seconds and needs about 7 GB of free memory; the first call in a
#' session also reads the table of pairs once.
#' @param x gene symbols, aliases, UniProt accessions or site keys (`GENE_N123`, `GENE_S123`); one entry per
#'   element, or several separated by newlines, commas, semicolons or whitespace (a pasted list). A
#'   data.frame is read as entries in its first column.
#' @param level `"auto"` (site level if any entry names a site, expanding gene entries to their reported sites; protein level for gene-only lists), `"protein"`
#'   (site keys contribute their protein) or `"site"` (every entry must be a site key)
#' @param filters the view applied to the row tables, see [ptm_default_filters()]
#' @param n_perm random lists drawn for the connectivity P; 0 skips the permutation and reports the
#'   analytic expectation only
#' @param seed the random seed for the permutation (`NULL` leaves the session's generator alone)
#' @return an object of class `ptm_list`, a list with: `level`; `members` (one row per entry: how it resolved,
#'   the site named if any, the status on each side, the number of reported pairs on each side,
#'   and how many internal and neighborhood rows it appears in); `internal` (pairs with a list member on both
#'   sides), `same_protein` and `neighbourhood` (pairs with a list member on one side): rows of reported pairs in
#'   the view, sorted on evidence, with `site_named_glyco` and `site_named_phospho`;
#'   `partner_sites` (one row per partner site: `site`, `gene`, `side` (the side
#'   it sits on), `members_reached` of `members_testable`, `tier_reach` of `tier_n` (how many units the site reaches
#'   among all significant pairs), `fold`, `p`, `q`, `n_pairs`, `members` (the members reached), `evidence_best`,
#'   `best_pair_id`; sorted on `p`); `convergent` (the same columns per partner protein, with `partner` and
#'   `partner_side`); `connectivity` (`members_in_tier`, `internal_pairs`, `expected_pairs`, `ratio`, `p_perm`,
#'   `n_perm`, `perm_median`); `universe` (the pairs the tests count: the rule and their number); `composition`
#'   (per table: rows, rows labeled with a tissue-composition program on both sides, the share); `withheld`;
#'   `provenance`; `note`.
#' @examples
#' \dontrun{
#' ptm_set_release_dir("path/to/release")   # the Parquet release, needed by ptm_list()
#' l <- ptm_list(c("EGFR", "ERBB2", "MET", "ITGB1", "CD44", "TFRC_N727"))
#' l                              # the summary
#' l$partner_sites[1:10, ]        # the sites the list converges on: fold, q, best pair
#' l$convergent[1:10, ]           # the same per partner protein
#' l$connectivity                 # is the list more connected than its degrees predict
#' ptm_list(c("MAPK1_T185", "MAPK1_Y187", "STAT3_Y705"))   # a phosphosite list: partner glycosites
#' ptm_plot_network(l)            # the internal pairs as a graph
#' }
#' @export
ptm_list <- function(x, level = c("auto", "protein", "site"), filters = ptm_default_filters(), n_perm = 200L, seed = 1L) {
  level <- match.arg(level); filters <- .ptm_filters(filters)
  if (!(is.numeric(n_perm) && length(n_perm) == 1L && is.finite(n_perm) && n_perm >= 0)) stop("n_perm is a whole number of random lists, 0 or more", call. = FALSE)
  if (!(is.numeric(seed) && length(seed) == 1L && is.finite(seed))) stop("seed is one whole number", call. = FALSE)
  ds <- .ptm_dataset(); schema <- ds$schema
  ent <- .ptm_entries(x)
  if (!nrow(ent)) stop("the list is empty: give gene symbols or site keys", call. = FALSE)
  is_site <- grepl("^.+_[NSTY][0-9]+$", toupper(ent$query))
  mixed <- level == "auto" && any(is_site) && !all(is_site)
  if (level == "auto") level <- if (any(is_site)) "site" else "protein"
  if (level == "site" && !all(is_site) && !mixed) stop("at site level every entry must be a site key (GENE_N123 or GENE_S123); ", paste(utils::head(ent$query[!is_site], 5), collapse = ", "), if (sum(!is_site) > 5) ", ..." else "", " are not. Pass level = \"protein\" to use their proteins", call. = FALSE)
  members <- do.call(rbind, lapply(seq_len(nrow(ent)), function(i) {
    q <- ent$query[i]; k <- if (is_site[i]) .ptm_site_key(q) else NULL
    r <- tryCatch(.ptm_resolve(if (is_site[i]) k$gene else q), error = function(e) list(gene = NA_character_, via = conditionMessage(e)))
    data.frame(query = q, kind = if (is_site[i]) "site" else "gene", gene = r$gene, resolved_via = if (is.na(r$gene) && r$via == "unknown") "not a recognized symbol" else r$via,
               site = if (is_site[i] && !is.na(r$gene)) paste0(r$gene, "_", sub("^.*_", "", k$key)) else NA_character_, side_named = if (is_site[i]) k$side else NA_character_, stringsAsFactors = FALSE)
  }))
  G <- unique(members$gene[!is.na(members$gene)]); S <- unique(members$site[!is.na(members$site) & !is.na(members$gene)])
  if (!length(G)) stop("no entry of the list is a recognized symbol, alias, UniProt accession or site key; ptm_status() says how names resolve", call. = FALSE)
  st <- ptm_status(G)
  for (nm in c("status_glyco", "status_phospho", "n_served_as_glyco", "n_served_as_phospho")) members[[nm]] <- if (nm %in% names(st)) st[[nm]][match(members$gene, st$gene)] else NA
  # n_reported: the reported pairs an entry brings, before any filter. A gene entry: its protein's pairs on both
  # sides. A site entry: that site's own pairs; a site in no reported pair is said, as ptm_site() says it.
  members$n_reported <- as.integer(ifelse(is.na(members$gene), NA, ifelse(is.na(members$n_served_as_glyco), 0L, members$n_served_as_glyco) + ifelse(is.na(members$n_served_as_phospho), 0L, members$n_served_as_phospho)))
  named <- which(members$kind == "site" & !is.na(members$site))
  if (length(named)) {
    ss <- unique(members$site[named])
    own <- as.data.frame(arrow::Scanner$create(ds, filter = .ptm_either_side("glyco_site", "phospho_site", ss, schema), projection = c("glyco_site", "phospho_site"))$ToTable(), stringsAsFactors = FALSE)
    cnt <- table(c(own$glyco_site[own$glyco_site %in% ss], own$phospho_site[own$phospho_site %in% ss]))
    members$n_reported[named] <- as.integer(ifelse(members$site[named] %in% names(cnt), cnt[members$site[named]], 0L))
    none <- unique(members$site[named][members$n_reported[named] == 0L])
    if (length(none)) message(paste(none, collapse = ", "), if (length(none) == 1L) " is" else " are", " in no reported pair of this release")
  }
  if (mixed) {
    expanded_genes <- unique(members$gene[members$kind == "gene" & !is.na(members$gene)])
    ef <- .ptm_either_side("glyco_gene", "phospho_gene", expanded_genes, schema)
    expansion <- as.data.frame(arrow::Scanner$create(ds, filter = ef, projection = c("glyco_site", "phospho_site", "glyco_gene", "phospho_gene"))$ToTable())
    S <- unique(c(S, expansion$glyco_site[expansion$glyco_gene %in% expanded_genes], expansion$phospho_site[expansion$phospho_gene %in% expanded_genes]))
  }
  # the unit a member is on the tier
  unit_of <- if (level == "site") members$site else members$gene
  units <- if (level == "site") S else G
  # the rows: every reported pair with a member on either side, in the view
  filt <- if (level == "site") .ptm_either_side("glyco_site", "phospho_site", S, schema) else .ptm_either_side("glyco_gene", "phospho_gene", G, schema)
  rows <- .ptm_scan(ds, filt, NULL, filters, NULL, list(list = ent$query, level = level, genes = G, sites = if (length(S)) S else NULL))
  prov <- attr(rows, "provenance"); withheld <- attr(rows, "withheld")
  rows$site_named_glyco <- rows$glyco_site %in% S; rows$site_named_phospho <- rows$phospho_site %in% S
  mg <- if (level == "site") rows$site_named_glyco else rows$glyco_gene %in% G
  mp <- if (level == "site") rows$site_named_phospho else rows$phospho_gene %in% G
  same <- rows$glyco_gene == rows$phospho_gene
  # each row table prints as ranked pairs do, with what it holds as its title
  strip <- function(r, what) { rownames(r) <- NULL; for (a in c("withheld", "n_matched", "n_view", "provenance")) attr(r, a) <- NULL; .ptm_as_pairs(r, what) }
  internal <- strip(rows[mg & mp & !same, , drop = FALSE], "pairs within the list"); same_protein <- strip(rows[same & (mg | mp), , drop = FALSE], "same-protein pairs of the list")
  neighbourhood <- strip(rows[xor(mg, mp) & !same, , drop = FALSE], "pairs with a partner outside the list")
  in_rows <- function(r, u) { if (level == "site") r$glyco_site == u | r$phospho_site == u else r$glyco_gene == u | r$phospho_gene == u }
  count_member <- function(tab, i) {
    if (mixed && members$kind[i] == "gene" && !is.na(members$gene[i])) return(sum(tab$glyco_gene == members$gene[i] | tab$phospho_gene == members$gene[i]))
    if (is.na(unit_of[i])) NA_integer_ else sum(in_rows(tab, unit_of[i]))
  }
  members$n_internal <- vapply(seq_len(nrow(members)), function(i) count_member(internal, i), integer(1))
  members$n_neighbourhood <- vapply(seq_len(nrow(members)), function(i) count_member(neighbourhood, i), integer(1))
  # the tests, on the tier as released
  tier <- .ptm_tier()
  tests <- .ptm_list_tests(tier, level, units, neighbourhood)
  connectivity <- .ptm_connectivity(tier, level, units, n_perm, seed)
  # the composition share of each table
  share <- function(r) { both <- if (nrow(r)) .ptm_glyco_labelled(r$glyco_programs, r$rep_glycans) & nzchar(ifelse(is.na(r$phospho_programs), "", r$phospho_programs)) else logical(0); c(rows = nrow(r), both_labelled = sum(both), share = if (nrow(r)) sum(both) / nrow(r) else NA_real_) }
  composition <- as.data.frame(rbind(internal = share(internal), same_protein = share(same_protein), neighbourhood = share(neighbourhood)))
  composition <- cbind(table = rownames(composition), composition, stringsAsFactors = FALSE); rownames(composition) <- NULL
  out <- list(level = level, mixed = mixed, sites = S, members = members, internal = internal, same_protein = same_protein, neighbourhood = neighbourhood,
              partner_sites = structure(tests$sites, class = c("ptm_partners", "data.frame")), convergent = structure(tests$proteins, class = c("ptm_partners", "data.frame")), connectivity = connectivity,
              universe = list(rule = .PTM_TEST_UNIVERSE, pairs = tier$M),
              composition = composition, withheld = withheld, provenance = prov,
              note = c(rows = "the row tables obey the view and are sorted on evidence, the only ranking; the score is symmetric and licenses no direction",
                       tests = "the partner-site, convergent-partner and connectivity tests count the significant reported pairs (q < 0.05) between two different proteins; a pair reported only for its perturbation support is not counted. They say a site or protein is reached by more members, or the list is more connected, than the members' own numbers of significant pairs predict, and nothing about mechanism; a partner reached by one member is a pair-level finding, already in the neighbourhood table",
                       composition = "a list from a tumor-versus-normal comparison is expected to be rich in tissue-composition labels; read the share before the pairs"))
  class(out) <- "ptm_list"
  out
}

# ---- entries --------------------------------------------------------------------------------------------
# the entries of a list: one per element, or several separated by newlines, commas, semicolons or whitespace
# (a pasted list). Every token is an entry; a data.frame is read as entries in its first column. A repeated
# entry is kept once, in the order it first appears.
.ptm_entries <- function(x) {
  if (is.data.frame(x)) x <- x[[1]]
  x <- as.character(x); x[is.na(x)] <- ""
  tok <- unlist(strsplit(x, "[\n\r,;[:space:]]+"), use.names = FALSE)
  tok <- tok[nzchar(tok)]
  data.frame(query = unique(tok), stringsAsFactors = FALSE)
}

# ---- the tier -------------------------------------------------------------------------------------------
# The pairs the tests count, reduced to what the tests need: every significant reported row (q_value < 0.05) with
# its two sites and two proteins as integer codes, same-protein rows dropped, and for every unit on each side its
# degree (pairs) and reach (distinct partner proteins; for a protein also distinct partner sites). Read once per
# release directory and kept in the package environment.
.PTM_TEST_UNIVERSE <- "significant reported pairs (q_value < 0.05) between two different proteins"
.ptm_tier <- function() {
  key <- paste(.ptm_release_dir(), .ptm_env$release_sha)
  if (!is.null(.ptm_env$tier) && identical(.ptm_env$tier$key, key)) return(.ptm_env$tier)
  ds <- .ptm_dataset()
  tab <- as.data.frame(arrow::Scanner$create(ds, projection = c("glyco_site", "phospho_site", "glyco_gene", "phospho_gene", "q_value"))$ToTable(), stringsAsFactors = FALSE)
  .ptm_env$tier <- .ptm_tier_from(tab, key, significant_only = TRUE)
  .ptm_env$tier
}
# the tier from a table of rows of reported pairs (any source: the Parquet, the build's TSV, the fixture) with the
# four site and gene columns and, for `significant_only`, q_value. The readers' tests pass significant_only = TRUE.
# The readers' tests and the website's List page file (tier_significant.json.gz, build/publish/web_tier.R) pass
# significant_only = TRUE; the default counts every reported row, which build/publish/reader_site_ids.R needs.
.ptm_tier_from <- function(tab, key = NA_character_, significant_only = FALSE) {
  if (isTRUE(significant_only)) {
    if (!"q_value" %in% names(tab)) stop("the table of pairs has no q_value column, so its significant pairs cannot be told apart", call. = FALSE)
    tab <- tab[is.finite(tab$q_value) & tab$q_value < 0.05, , drop = FALSE]
  }
  tab <- tab[tab$glyco_gene != tab$phospho_gene, , drop = FALSE]
  gs <- factor(tab$glyco_site); ps <- factor(tab$phospho_site); gg <- factor(tab$glyco_gene); pg <- factor(tab$phospho_gene)
  distinct <- function(a, b) { nb <- nlevels(b); code <- (as.numeric(a) - 1) * nb + as.numeric(b); u <- unique(code); list(a = as.integer((u - 1) %/% nb + 1), b = as.integer((u - 1) %% nb + 1)) }
  d_gg_pg <- distinct(gg, pg); d_ps_gg <- distinct(ps, gg); d_gs_pg <- distinct(gs, pg)
  tier <- list(key = key, M = nrow(tab), significant_only = isTRUE(significant_only),
               rows = list(gs = as.integer(gs), ps = as.integer(ps), gg = as.integer(gg), pg = as.integer(pg)),
               sites_glyco = levels(gs), sites_phospho = levels(ps), genes_glyco = levels(gg), genes_phospho = levels(pg),
               site_gene_glyco = levels(gg)[as.integer(gg)[match(seq_len(nlevels(gs)), as.integer(gs))]],
               site_gene_phospho = levels(pg)[as.integer(pg)[match(seq_len(nlevels(ps)), as.integer(ps))]],
               gene_glyco = data.frame(deg = tabulate(gg, nlevels(gg)), reach_genes = tabulate(d_gg_pg$a, nlevels(gg)), reach_sites = tabulate(d_ps_gg$b, nlevels(gg))),
               gene_phospho = data.frame(deg = tabulate(pg, nlevels(pg)), reach_genes = tabulate(d_gg_pg$b, nlevels(pg)), reach_sites = tabulate(d_gs_pg$b, nlevels(pg))),
               site_glyco = data.frame(deg = tabulate(gs, nlevels(gs)), reach_genes = tabulate(d_gs_pg$a, nlevels(gs))),
               site_phospho = data.frame(deg = tabulate(ps, nlevels(ps)), reach_genes = tabulate(d_ps_gg$a, nlevels(ps))))
  rownames(tier$gene_glyco) <- tier$genes_glyco; rownames(tier$gene_phospho) <- tier$genes_phospho
  rownames(tier$site_glyco) <- tier$sites_glyco; rownames(tier$site_phospho) <- tier$sites_phospho
  tier
}
# the tier as the website's List page reads it (its data/tier_significant.json.gz, written when the release is built,
# and by the fixture): the rule the pairs were counted under, the universe sizes, and per gene and per site the same
# counts under the same names
.ptm_tier_json <- function(tier, release = NA_character_, built_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")) {
  genes <- union(tier$genes_glyco, tier$genes_phospho)
  ig <- match(genes, tier$genes_glyco); ip <- match(genes, tier$genes_phospho)
  z <- function(v) { v[is.na(v)] <- 0L; as.integer(v) }
  gene_rows <- lapply(seq_along(genes), function(i) c(z(tier$gene_glyco$deg[ig[i]]), z(tier$gene_glyco$reach_genes[ig[i]]), z(tier$gene_glyco$reach_sites[ig[i]]),
                                                       z(tier$gene_phospho$deg[ip[i]]), z(tier$gene_phospho$reach_genes[ip[i]]), z(tier$gene_phospho$reach_sites[ip[i]])))
  sites <- c(tier$sites_glyco, tier$sites_phospho)
  site_rows <- c(lapply(seq_along(tier$sites_glyco), function(i) c(tier$site_glyco$deg[i], tier$site_glyco$reach_genes[i])),
                 lapply(seq_along(tier$sites_phospho), function(i) c(tier$site_phospho$deg[i], tier$site_phospho$reach_genes[i])))
  list(release = release, built_utc = built_utc, significant_only = isTRUE(tier$significant_only),
       rule = if (isTRUE(tier$significant_only)) .PTM_TEST_UNIVERSE else "reported pairs between two different proteins", M = tier$M,
       n = list(genes_glyco = length(tier$genes_glyco), genes_phospho = length(tier$genes_phospho), sites_glyco = length(tier$sites_glyco), sites_phospho = length(tier$sites_phospho)),
       columns = list(genes = c("deg_glyco", "reach_genes_glyco", "reach_sites_glyco", "deg_phospho", "reach_genes_phospho", "reach_sites_phospho"), sites = c("deg", "reach_genes")),
       genes = stats::setNames(gene_rows, genes), sites = stats::setNames(site_rows, sites))
}
# the codes of a set of units on one side of the tier (NA when absent), and a value per unit (0 when absent)
.ptm_codes <- function(tier, level, side, units) match(units, if (level == "site") (if (side == "glyco") tier$sites_glyco else tier$sites_phospho) else (if (side == "glyco") tier$genes_glyco else tier$genes_phospho))
.ptm_deg <- function(tier, level, side, units) { t <- tier[[paste0(if (level == "site") "site_" else "gene_", side)]]; k <- t$deg[.ptm_codes(tier, level, side, units)]; k[is.na(k)] <- 0L; k }

# ---- the convergence tests -----------------------------------------------------------------------------
# one family: for every tier row with a member unit on the member side and no member on the other, the partner
# unit's code `pu` and the member unit's rank `mr`; `reach` the partner units' distinct member-side units on the
# whole tier by code; `N` the member-side universe; `n` the testable members; `mnames` the member units by rank
.ptm_converge <- function(pu, mr, reach, pnames, N, n, mnames) {
  if (!length(pu) || !n) return(NULL)
  nu <- length(mnames); code <- (as.numeric(pu) - 1) * nu + mr; u <- unique(code)
  p_of <- as.integer((u - 1) %/% nu + 1); m_of <- as.integer((u - 1) %% nu + 1)
  x <- tabulate(p_of, nbins = length(reach)); np <- tabulate(pu, nbins = length(reach)); i <- which(x > 0)
  reached <- vapply(split(mnames[m_of], p_of)[as.character(i)], function(v) paste(sort(v), collapse = ";"), character(1))
  data.frame(partner = pnames[i], members_reached = x[i], members_testable = n, tier_reach = reach[i], tier_n = N, n_pairs = np[i], members = unname(reached), stringsAsFactors = FALSE)
}
# the two families (partner sites, partner proteins), each over both sides, with fold, P, BH q and the best
# pair in the view
.ptm_list_tests <- function(tier, level, units, neighbourhood) {
  r <- tier$rows
  one_side <- function(member_side) {
    partner_side <- if (member_side == "glyco") "phospho" else "glyco"
    ucol <- if (level == "site") (if (member_side == "glyco") "gs" else "ps") else (if (member_side == "glyco") "gg" else "pg")
    ocol <- if (level == "site") (if (member_side == "glyco") "ps" else "gs") else (if (member_side == "glyco") "pg" else "gg")
    universe <- if (level == "site") (if (member_side == "glyco") tier$sites_glyco else tier$sites_phospho) else (if (member_side == "glyco") tier$genes_glyco else tier$genes_phospho)
    codes <- .ptm_codes(tier, level, member_side, units); present <- units[!is.na(codes)]; codes <- codes[!is.na(codes)]
    n <- length(present); N <- length(universe)
    if (!n) return(NULL)
    mask_m <- logical(N); mask_m[codes] <- TRUE
    other_codes <- .ptm_codes(tier, level, partner_side, units); mask_o <- logical(length(if (level == "site") (if (partner_side == "glyco") tier$sites_glyco else tier$sites_phospho) else (if (partner_side == "glyco") tier$genes_glyco else tier$genes_phospho))); mask_o[other_codes[!is.na(other_codes)]] <- TRUE
    keep <- mask_m[r[[ucol]]] & !mask_o[r[[ocol]]]
    mr <- match(r[[ucol]][keep], codes)
    site_col <- if (partner_side == "glyco") "gs" else "ps"; gene_col <- if (partner_side == "glyco") "gg" else "pg"
    st <- tier[[paste0("site_", partner_side)]]; gt <- tier[[paste0("gene_", partner_side)]]
    reach_site <- if (level == "site") st$deg else st$reach_genes; reach_gene <- if (level == "site") gt$reach_sites else gt$reach_genes
    s <- .ptm_converge(r[[site_col]][keep], mr, reach_site, rownames(st), N, n, present)
    g <- .ptm_converge(r[[gene_col]][keep], mr, reach_gene, rownames(gt), N, n, present)
    if (!is.null(s)) { s$gene <- tier[[paste0("site_gene_", partner_side)]][match(s$partner, rownames(st))]; s$side <- partner_side }
    if (!is.null(g)) g$partner_side <- partner_side
    list(sites = s, proteins = g)
  }
  a <- one_side("glyco"); b <- one_side("phospho")
  finish <- function(tab, key) {
    empty <- data.frame(partner = character(0), members_reached = integer(0), members_testable = integer(0), tier_reach = integer(0), tier_n = integer(0), n_pairs = integer(0), members = character(0), stringsAsFactors = FALSE)
    if (is.null(tab) || !nrow(tab)) tab <- empty
    if (key == "site") { if (!"gene" %in% names(tab)) { tab$gene <- character(0); tab$side <- character(0) }; names(tab)[names(tab) == "partner"] <- "site"; side_col <- "side" }
    else { if (!"partner_side" %in% names(tab)) tab$partner_side <- character(0); side_col <- "partner_side" }
    tab$fold <- (tab$members_reached / tab$members_testable) / (tab$tier_reach / tab$tier_n)
    tab$p <- stats::phyper(tab$members_reached - 1, tab$tier_reach, tab$tier_n - tab$tier_reach, tab$members_testable, lower.tail = FALSE)
    tab$q <- stats::p.adjust(tab$p, method = "BH")
    v <- .ptm_partner_view(neighbourhood, level, units, key)
    i <- match(paste(tab[[key]], tab[[side_col]]), v$key)
    tab$evidence_best <- v$evidence_best[i]; tab$best_pair_id <- v$best_pair_id[i]
    cols <- if (key == "site") c("site", "gene", "side") else c("partner", "partner_side")
    tab <- tab[order(tab$p, -tab$fold, -ifelse(is.na(tab$evidence_best), -Inf, tab$evidence_best), tab[[key]], method = "radix"), c(cols, "members_reached", "members_testable", "tier_reach", "tier_n", "fold", "p", "q", "n_pairs", "members", "evidence_best", "best_pair_id")]
    rownames(tab) <- NULL; tab
  }
  list(sites = finish(rbind(a$sites, b$sites), "site"), proteins = finish(rbind(a$proteins, b$proteins), "partner"))
}
# what the view says about each partner: the strongest row (rows are sorted on evidence, so the first is the best)
.ptm_partner_view <- function(nb, level, units, key) {
  if (!nrow(nb)) return(data.frame(key = character(0), evidence_best = numeric(0), best_pair_id = character(0), stringsAsFactors = FALSE))
  on_glyco <- if (level == "site") nb$glyco_site %in% units else nb$glyco_gene %in% units
  col <- if (key == "site") "site" else "gene"
  partner <- ifelse(on_glyco, nb[[paste0("phospho_", col)]], nb[[paste0("glyco_", col)]]); side <- ifelse(on_glyco, "phospho", "glyco")
  pk <- paste(partner, side); first <- !duplicated(pk)
  data.frame(key = pk[first], evidence_best = nb$evidence[first], best_pair_id = nb$pair_id[first], stringsAsFactors = FALSE)
}

# ---- connectivity ---------------------------------------------------------------------------------------
# internal pairs observed on the tier, the Chung-Lu expectation from the members' degrees, and a permutation P
# from random lists drawn to match those degrees (every unit of the level binned by its total degree; each
# member replaced by a random unit from its bin)
.ptm_connectivity <- function(tier, level, units, n_perm = 200L, seed = 1L) {
  r <- tier$rows; gcol <- if (level == "site") "gs" else "gg"; pcol <- if (level == "site") "ps" else "pg"
  ng <- length(if (level == "site") tier$sites_glyco else tier$genes_glyco); np <- length(if (level == "site") tier$sites_phospho else tier$genes_phospho)
  count <- function(set) { cg <- .ptm_codes(tier, level, "glyco", set); cp <- .ptm_codes(tier, level, "phospho", set)
    mg <- logical(ng); mg[cg[!is.na(cg)]] <- TRUE; mp <- logical(np); mp[cp[!is.na(cp)]] <- TRUE; sum(mg[r[[gcol]]] & mp[r[[pcol]]]) }
  expected <- function(set) {
    a <- .ptm_deg(tier, level, "glyco", set); b <- .ptm_deg(tier, level, "phospho", set)
    gene <- if (level == "site") sub("_[NSTY][0-9]+$", "", set) else set
    ga <- tapply(a, gene, sum); gb <- tapply(b, gene, sum)
    if (tier$M) (sum(as.numeric(a)) * sum(as.numeric(b)) - sum(as.numeric(ga) * as.numeric(gb))) / tier$M else NA_real_
  }
  kg <- .ptm_deg(tier, level, "glyco", units); kp <- .ptm_deg(tier, level, "phospho", units); present <- units[kg + kp > 0]
  obs <- count(units); expct <- expected(units)
  out <- list(members_in_tier = length(present), internal_pairs = as.integer(obs), expected_pairs = expct, ratio = if (expct > 0) obs / expct else NA_real_,
              p_perm = NA_real_, n_perm = 0L, perm_median = NA_real_)
  if (n_perm > 0 && length(present) >= 2) {
    all_units <- if (level == "site") c(tier$sites_glyco, tier$sites_phospho) else union(tier$genes_glyco, tier$genes_phospho)
    tot <- .ptm_deg(tier, level, "glyco", all_units) + .ptm_deg(tier, level, "phospho", all_units)
    br <- unique(stats::quantile(tot, probs = seq(0, 1, by = 0.1))); bin <- if (length(br) < 2) rep(1L, length(tot)) else cut(tot, breaks = br, include.lowest = TRUE, labels = FALSE)
    by_bin <- split(all_units, bin); need <- table(bin[match(present, all_units)])
    if (!is.null(seed)) { old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL
      on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()), add = TRUE); set.seed(seed) }
    draw <- function() unlist(lapply(names(need), function(b) { pool <- by_bin[[b]]; if (length(pool) <= need[[b]]) pool else sample(pool, need[[b]]) }), use.names = FALSE)
    perm <- vapply(seq_len(n_perm), function(i) count(draw()), numeric(1))
    out$p_perm <- (1 + sum(perm >= obs)) / (1 + n_perm); out$n_perm <- as.integer(n_perm); out$perm_median <- stats::median(perm)
  }
  out
}
