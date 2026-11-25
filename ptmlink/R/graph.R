# The graph functions: reported pairs as an igraph object, a static plot, and node and edge tables for
# Cytoscape or any other viewer.
#
# One edge table underlies all three. An edge joins the protein whose glycosites are involved to the
# protein whose phosphosites are involved (gene level), or a glycosite to a phosphosite (site level).
# The graph is undirected because the statistic is symmetric. At the gene level an edge carries the
# strongest of its site pairs, as the website's network does; nothing is summed or recomputed.
#
# Input is whatever the package already returns: the edge table of ptm_network() (one or several,
# stacked), rows of reported pairs from ptm_gene(), ptm_site(), ptm_query() or ptm_list()$internal, or a ptm_list
# object (its internal pairs).

# the standard edge and node tables from any of the inputs above
.ptm_edge_table <- function(x, level = c("gene", "site")) {
  level <- match.arg(level)
  query_genes <- character(0)
  if (inherits(x, "ptm_list")) { query_genes <- unique(x$members$gene[!is.na(x$members$gene)]); x <- x$internal }
  if (!is.data.frame(x)) stop("x is the table of ptm_network(), a result of ptm_gene(), ptm_site() or ptm_query(), or a ptm_list()", call. = FALSE)
  ecols <- c("glyco_gene", "phospho_gene", "evidence", "rho", "k", "n_site_pairs", "anchored", "sign", "shadow_clean", "self_coupling", "glyco_site", "phospho_site", "pair_id")
  if (!nrow(x)) {
    e <- data.frame(from = character(0), to = character(0), glyco_gene = character(0), phospho_gene = character(0), evidence = numeric(0), rho = numeric(0), k = integer(0), n_site_pairs = integer(0), anchored = logical(0), sign = character(0),
                    shadow_clean = logical(0), self_coupling = numeric(0), glyco_site = character(0), phospho_site = character(0), pair_id = character(0), stringsAsFactors = FALSE)
    nodes <- data.frame(name = query_genes, kind = rep("protein", length(query_genes)), gene = query_genes, query = rep(TRUE, length(query_genes)), n_edges = rep(0L, length(query_genes)), stringsAsFactors = FALSE)
    return(list(edges = e, nodes = nodes, level = level))
  }
  if ("partner" %in% names(x) && "role" %in% names(x)) {
    if (level == "site") stop("site level needs rows of reported pairs (ptm_gene, ptm_site, ptm_query or ptm_list()$internal); a network table has one edge per protein pair", call. = FALSE)
    need <- c("gene", "partner", "role", "evidence_best", "rho_at_best", "k_at_best", "n_site_pairs", "anch_any", "sign_mixed", "shadow_clean_at_best", "self_coupling_at_best", "best_glyco_site", "best_phospho_site")
    miss <- setdiff(need, names(x)); if (length(miss)) stop("a network table needs the columns ", paste(miss, collapse = ", "), call. = FALSE)
    query_genes <- unique(c(query_genes, x$gene))
    e <- data.frame(glyco_gene = ifelse(x$role == "glyco", x$gene, x$partner), phospho_gene = ifelse(x$role == "glyco", x$partner, x$gene),
                    evidence = x$evidence_best, rho = x$rho_at_best, k = x$k_at_best, n_site_pairs = x$n_site_pairs, anchored = x$anch_any %in% TRUE,
                    sign = ifelse(x$sign_mixed %in% TRUE, "mixed", ifelse(x$rho_at_best < 0, "negative", "positive")),
                    shadow_clean = x$shadow_clean_at_best, self_coupling = x$self_coupling_at_best, glyco_site = x$best_glyco_site, phospho_site = x$best_phospho_site,
                    pair_id = NA_character_, stringsAsFactors = FALSE)
    e <- e[order(-e$evidence, method = "radix"), ]; e <- e[!duplicated(paste(e$glyco_gene, e$phospho_gene)), ]   # several stacked networks: one edge per protein pair, the strongest kept
  } else {
    need <- c("glyco_site", "phospho_site", "glyco_gene", "phospho_gene", "evidence"); miss <- setdiff(need, names(x))
    if (length(miss)) stop("rows of reported pairs need the columns ", paste(miss, collapse = ", "), call. = FALSE)
    col <- function(nm, default = NA) if (nm %in% names(x)) x[[nm]] else rep(default, nrow(x))
    rows <- data.frame(glyco_gene = x$glyco_gene, phospho_gene = x$phospho_gene, evidence = x$evidence, rho = as.numeric(col("rho_combined")), k = as.integer(col("k_cohorts")),
                       anchored = col("anch_any") %in% TRUE, shadow_clean = as.logical(col("shadow_clean")), self_coupling = as.numeric(col("self_coupling")),
                       glyco_site = x$glyco_site, phospho_site = x$phospho_site, pair_id = as.character(col("pair_id")), stringsAsFactors = FALSE)
    rows <- rows[order(-rows$evidence, method = "radix"), ]
    if (level == "site") { e <- rows; e$n_site_pairs <- 1L; e$sign <- ifelse(is.na(e$rho), NA_character_, ifelse(e$rho < 0, "negative", "positive")) }
    else {
      key <- paste(rows$glyco_gene, rows$phospho_gene); sp <- split(seq_len(nrow(rows)), key)
      e <- do.call(rbind, lapply(sp, function(i) { b <- i[1]; r <- rows[b, ]; pos <- sum(rows$rho[i] >= 0, na.rm = TRUE); neg <- sum(rows$rho[i] < 0, na.rm = TRUE)
        r$n_site_pairs <- length(i); r$anchored <- any(rows$anchored[i]); r$sign <- if (pos > 0 && neg > 0) "mixed" else if (neg > 0) "negative" else "positive"; r }))
      e <- e[order(-e$evidence, method = "radix"), ]
    }
  }
  e <- e[, ecols]
  if (level == "gene") {
    e <- cbind(from = e$glyco_gene, to = e$phospho_gene, e, stringsAsFactors = FALSE)
    genes <- unique(c(e$from, e$to, query_genes))
    nodes <- data.frame(name = genes, kind = "protein", gene = genes, query = genes %in% query_genes,
                        n_edges = vapply(genes, function(g) sum(e$from == g | e$to == g), integer(1)), stringsAsFactors = FALSE)
  } else {
    e <- cbind(from = e$glyco_site, to = e$phospho_site, e, stringsAsFactors = FALSE)
    sites <- unique(c(e$from, e$to)); gl <- sites %in% e$from
    nodes <- data.frame(name = sites, kind = ifelse(gl, "glycosite", "phosphosite"), gene = ifelse(gl, e$glyco_gene[match(sites, e$from)], e$phospho_gene[match(sites, e$to)]),
                        query = FALSE, n_edges = vapply(sites, function(s) sum(e$from == s | e$to == s), integer(1)), stringsAsFactors = FALSE)
    nodes$query <- nodes$gene %in% query_genes
  }
  rownames(e) <- NULL; rownames(nodes) <- NULL
  list(edges = e, nodes = nodes, level = level)
}

#' Reported pairs as an igraph object
#' @param x the edge table of [ptm_network()] (one or several stacked), rows of reported pairs from [ptm_gene()], [ptm_site()],
#'   [ptm_query()] or `ptm_list()$internal`, or a `ptm_list` (its internal pairs)
#' @param level "gene": one undirected edge per protein pair, from the protein whose glycosites are involved to the
#'   protein whose phosphosites are involved, carrying the strongest site pair; "site": one edge per reported pair
#'   between a glycosite node and a phosphosite node (rows of reported pairs only)
#' @return an undirected igraph graph. Vertex attributes: `kind` ("protein", "glycosite" or "phosphosite"), `gene`,
#'   `query` (a queried or listed protein), `n_edges`. Edge attributes: `glyco_gene`, `phospho_gene`, `evidence`,
#'   `rho`, `k`, `n_site_pairs`, `anchored` (perturbation evidence on any site pair), `sign` ("positive", "negative"
#'   or "mixed" across the site pairs), `shadow_clean`, `self_coupling`, `glyco_site`, `phospho_site` (the strongest
#'   site pair at the gene level), `pair_id` (site level). Graph attributes: `level`, `note`.
#' @examples
#' \dontrun{
#' g <- ptm_igraph(ptm_network("EGFR"))
#' igraph::E(g)$evidence
#' ptm_igraph(ptm_gene("EGFR", "glyco"), level = "site")
#' }
#' @export
ptm_igraph <- function(x, level = c("gene", "site")) {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("ptm_igraph() needs the igraph package: install.packages(\"igraph\")", call. = FALSE)
  t <- .ptm_edge_table(x, level)
  g <- igraph::graph_from_data_frame(t$edges, directed = FALSE, vertices = t$nodes)
  g <- igraph::set_graph_attr(g, "level", t$level)
  igraph::set_graph_attr(g, "note", "undirected: the statistic is symmetric and licenses no direction; an edge at the gene level carries its strongest site pair, not a consensus")
}

#' Plot reported pairs as a network
#'
#' A static drawing for a figure or a quick look. Edge width follows evidence; a dashed edge is a negative
#' correlation and a dotted one has site pairs of both signs; an orange edge carries perturbation evidence.
#' Queried or listed proteins are filled blue; at the site level glycosites are blue and phosphosites orange.
#' The layout is force-directed and carries no meaning; the website's network view is the designed drawing.
#' Labels are placed so that they overlap neither each other nor any node; a label moved away from its node is
#' joined to it by a thin line. The drawing has no title unless `main` is given.
#'
#' A drawing stays readable up to about a hundred links, so two limits apply, each said in a message when it
#' acts. A table from [ptm_network()] is drawn at the partners it marks `drawn` (its `top_n` strongest, 25 by
#' default; `ptm_network(gene, top_n = NULL)` marks all). Any input with more than `max_edges` links is drawn at
#' its `max_edges` strongest. [ptm_igraph()] and [ptm_write_network()] always keep every link.
#' @inheritParams ptm_igraph
#' @param layout an igraph layout function, or a two-column matrix of node positions (default [igraph::layout_with_fr])
#' @param label_cex label size
#' @param main an optional title; none by default
#' @param max_edges draw at most this many links, the strongest by evidence (default 100); `NULL` draws all
#' @param ... passed to [igraph::plot.igraph()]
#' @return the igraph graph, invisibly
#' @examples
#' \dontrun{
#' ptm_plot_network(ptm_network("EGFR", top_n = 15))     # the 15 strongest partners of EGFR
#' ptm_plot_network(ptm_gene("EGFR", "glyco"), level = "site", max_edges = 40)
#' ptm_plot_network(ptm_list(c("EGFR", "ERBB2", "MET", "ITGB1")))
#' }
#' @export
ptm_plot_network <- function(x, level = c("gene", "site"), layout = NULL, label_cex = 0.8, main = NULL, max_edges = 100L, ...) {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("ptm_plot_network() needs the igraph package: install.packages(\"igraph\")", call. = FALSE)
  if (!is.null(max_edges) && !(is.numeric(max_edges) && length(max_edges) == 1L && is.finite(max_edges) && max_edges >= 1 && max_edges == round(max_edges))) stop("max_edges is a whole number, 1 or more, or NULL to draw every link", call. = FALSE)
  if (is.data.frame(x) && all(c("partner", "role", "drawn") %in% names(x)) && nrow(x) && !any(x$drawn %in% TRUE))
    stop("the network marks no partner for drawing (top_n = 0); ptm_network(gene, top_n = 25) marks the strongest 25", call. = FALSE)
  # a network table is drawn at the partners it marks (top_n), never at every partner of the gene
  if (is.data.frame(x) && all(c("partner", "role", "drawn") %in% names(x)) && any(x$drawn %in% TRUE) && !all(x$drawn %in% TRUE)) {
    message("Drawing the ", sum(x$drawn %in% TRUE), " strongest of ", formatC(nrow(x), format = "d", big.mark = ","), " partner links (top_n); ptm_network(gene, top_n = NULL) marks all.")
    x <- x[x$drawn %in% TRUE, , drop = FALSE]
  }
  g <- if (inherits(x, "igraph")) x else ptm_igraph(x, level)
  if (igraph::ecount(g) == 0) stop("there is no link to draw: no pair passes the filters", call. = FALSE)
  if (!is.null(max_edges) && igraph::ecount(g) > max_edges) {
    ev <- igraph::edge_attr(g, "evidence"); ev[is.na(ev)] <- -Inf
    drop <- order(-ev)[-seq_len(max_edges)]
    message("Drawing the ", max_edges, " strongest of ", formatC(igraph::ecount(g), format = "d", big.mark = ","), " links; max_edges = NULL draws all.")
    g <- igraph::delete_edges(g, drop)
    lone <- igraph::degree(g) == 0 & !(igraph::vertex_attr(g, "query") %in% TRUE)
    if (any(lone)) g <- igraph::delete_vertices(g, which(lone))
  }
  lev <- igraph::graph_attr(g, "level"); if (is.null(lev)) lev <- "gene"
  ev <- igraph::edge_attr(g, "evidence"); rng <- range(ev, na.rm = TRUE)
  width <- if (diff(rng) > 0) 0.6 + 3.4 * (ev - rng[1]) / diff(rng) else 2
  anch <- igraph::edge_attr(g, "anchored") %in% TRUE; sgn <- igraph::edge_attr(g, "sign")
  ecol <- ifelse(anch, "#D55E00", "grey45"); elty <- ifelse(sgn %in% "negative", 2L, ifelse(sgn %in% "mixed", 3L, 1L))
  kind <- igraph::vertex_attr(g, "kind"); qry <- igraph::vertex_attr(g, "query") %in% TRUE
  vcol <- if (lev == "site") ifelse(kind == "glycosite", "#0072B2", "#E69F00") else ifelse(qry, "#0072B2", "grey82")
  vsize <- if (lev == "site") rep(5, igraph::vcount(g)) else ifelse(qry, 9, 6)
  # node positions scaled to [-1, 1], the units igraph draws in (a node's radius is its size / 200)
  if (is.null(layout)) layout <- igraph::layout_with_fr
  xy <- if (is.function(layout)) layout(g) else as.matrix(layout)
  if (!is.numeric(xy) || nrow(xy) != igraph::vcount(g) || ncol(xy) < 2L) stop("layout is a function of the graph, or a matrix with one row per node (", igraph::vcount(g), " here) and two columns", call. = FALSE)
  xy <- igraph::norm_coords(xy, -1, 1, -1, 1); xy[!is.finite(xy)] <- 0
  # the legends sit in a band below the nodes; the labels get a margin on every side
  leg_lines <- list(args = list(bty = "n", cex = 0.75, lwd = c(3, 1, 1, 1), lty = c(1, 1, 2, 3), col = c("grey45", "#D55E00", "grey45", "grey45"),
                                legend = c("width: evidence (-log10 P)", "perturbation evidence", "negative correlation", "site pairs of both signs")))
  leg_nodes <- list(bty = "n", cex = 0.75, pch = 21, pt.bg = if (lev == "site") c("#0072B2", "#E69F00") else c("#0072B2", "grey88"), pt.cex = 1.4,
                    legend = if (lev == "site") c("glycosite", "phosphosite") else c("queried or listed protein", "partner"))
  old <- graphics::par(mar = c(0.5, 0.5, if (is.null(main)) 0.5 else 3, 0.5)); on.exit(graphics::par(old), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(c(-1, 1), c(-1, 1), asp = 1)
  band <- do.call(graphics::legend, c(list("bottomleft", plot = FALSE), leg_lines$args))$rect$h
  graphics::plot.window(c(-1.3, 1.3), c(-1.2 - band, 1.2), asp = 1)
  if (!is.null(main)) graphics::title(main = main)
  igraph::plot.igraph(g, layout = xy, rescale = FALSE, add = TRUE, vertex.color = vcol, vertex.frame.color = "white", vertex.size = vsize,
                      vertex.label = NA, edge.width = width, edge.color = ecol, edge.lty = elty, ...)
  usr <- graphics::par("usr")
  lg1 <- do.call(graphics::legend, c(list("bottomleft"), leg_lines$args)); lg2 <- do.call(graphics::legend, c(list("bottomright"), leg_nodes))
  floor_y <- max(lg1$rect$top, lg2$rect$top)
  lab <- igraph::vertex_attr(g, "name")
  w <- graphics::strwidth(lab, cex = label_cex, family = "sans") * 1.06 + 0.012
  h <- graphics::strheight(lab, cex = label_cex, family = "sans") * 1.35 + 0.012
  r <- vsize / 200
  pos <- .ptm_place_labels(xy, w, h, r, c(usr[1], usr[2], floor_y, usr[4]))
  # a leader line from the node's edge to the nearest point of a label that was moved off its node
  qx <- pmin(pmax(xy[, 1], pos[, 1] - w / 2), pos[, 1] + w / 2); qy <- pmin(pmax(xy[, 2], pos[, 2] - h / 2), pos[, 2] + h / 2)
  gap <- sqrt((qx - xy[, 1])^2 + (qy - xy[, 2])^2); far <- gap > r + 0.04
  if (any(far)) {
    sx <- xy[far, 1] + (qx[far] - xy[far, 1]) * r[far] / gap[far]; sy <- xy[far, 2] + (qy[far] - xy[far, 2]) * r[far] / gap[far]
    graphics::segments(sx, sy, qx[far], qy[far], col = "grey60", lwd = 0.6)
  }
  graphics::rect(pos[, 1] - w / 2, pos[, 2] - h / 2, pos[, 1] + w / 2, pos[, 2] + h / 2, col = grDevices::adjustcolor("white", 0.8), border = NA)
  graphics::text(pos[, 1], pos[, 2], lab, cex = label_cex, family = "sans", col = "black")
  invisible(g)
}

# Label positions that overlap neither each other nor any node. Each label starts just outside its node, on the
# side facing away from the middle of the drawing, and is pushed apart from overlapping labels (along the axis
# that needs the smaller move) and away from overlapping nodes, until nothing overlaps or `rounds` pass.
# xy: node centres; w, h: label box sizes; r: node radii; region: x0, x1, y0, y1 the labels must stay within.
.ptm_place_labels <- function(xy, w, h, r, region, pad = 0.008, rounds = 500L) {
  n <- nrow(xy); x0 <- xy[, 1]; y0 <- xy[, 2]
  r <- rep_len(r, n)
  dx <- x0 - mean(range(x0)); dy <- y0 - mean(range(y0)); len <- sqrt(dx^2 + dy^2)
  ux <- ifelse(len > 1e-9, dx / len, 0); uy <- ifelse(len > 1e-9, dy / len, 1)
  off <- r + pad + 0.5 * (abs(ux) * w + abs(uy) * h)
  x <- x0 + ux * off; y <- y0 + uy * off
  clamp <- function(x, y) list(x = pmin(pmax(x, region[1] + w / 2), region[2] - w / 2), y = pmin(pmax(y, region[3] + h / 2), region[4] - h / 2))
  p <- clamp(x, y); x <- p$x; y <- p$y
  if (n < 2L || n > 600L) return(cbind(x, y))
  X0 <- matrix(x0, n, n, byrow = TRUE); Y0 <- matrix(y0, n, n, byrow = TRUE); R <- matrix(r, n, n, byrow = TRUE)
  tie <- ifelse(row(X0) < col(X0), -1, 1)
  for (it in seq_len(rounds)) {
    # label against label
    ddx <- outer(x, x, "-"); ddy <- outer(y, y, "-")
    ox <- outer(w, w, "+") / 2 + pad - abs(ddx); oy <- outer(h, h, "+") / 2 + pad - abs(ddy)
    hit <- ox > 0 & oy > 0; diag(hit) <- FALSE
    sx <- ifelse(ddx == 0, tie, sign(ddx)); sy <- ifelse(ddy == 0, tie, sign(ddy))
    mx <- rowSums(ifelse(hit & ox < oy, sx * ox / 2, 0)); my <- rowSums(ifelse(hit & ox >= oy, sy * oy / 2, 0))
    # label against node: the box's nearest point to each node centre must lie outside that node
    qx <- pmin(pmax(X0, x - w / 2), x + w / 2); qy <- pmin(pmax(Y0, y - h / 2), y + h / 2)
    ov <- R + pad - sqrt((qx - X0)^2 + (qy - Y0)^2); nhit <- ov > 0
    if (!any(hit) && !any(nhit)) break
    vx <- x - X0; vy <- y - Y0; vl <- pmax(sqrt(vx^2 + vy^2), 1e-9)
    mx <- mx + rowSums(ifelse(nhit, vx / vl * ov, 0)); my <- my + rowSums(ifelse(nhit, vy / vl * ov, 0))
    p <- clamp(x + mx * 0.9, y + my * 0.9); x <- p$x; y <- p$y
  }
  cbind(x, y)
}

#' Write reported pairs as node and edge tables, or GraphML, for Cytoscape and other viewers
#' @inheritParams ptm_igraph
#' @param path the file stem: `<path>_nodes.tsv` and `<path>_edges.tsv` are written for "tsv", `<path>.graphml` for "graphml"
#' @param format "tsv" (two tab-separated tables; Cytoscape's table import reads the edge table with `from` and `to` as
#'   source and target) or "graphml" (needs igraph)
#' @return the paths written, invisibly
#' @examples
#' \dontrun{
#' # writes egfr_network_nodes.tsv and egfr_network_edges.tsv
#' ptm_write_network(ptm_network("EGFR"), "egfr_network")
#' }
#' @export
ptm_write_network <- function(x, path, level = c("gene", "site"), format = c("tsv", "graphml")) {
  format <- match.arg(format)
  if (!(is.character(path) && length(path) == 1L && !is.na(path) && nzchar(path))) stop("path is the file name to write, without its extension, such as file.path(tempdir(), \"egfr\")", call. = FALSE)
  if (!dir.exists(dirname(path))) stop("no such directory: ", dirname(path), call. = FALSE)
  if (format == "graphml") {
    g <- ptm_igraph(x, level); f <- paste0(path, ".graphml"); igraph::write_graph(g, f, format = "graphml"); return(invisible(f))
  }
  t <- .ptm_edge_table(x, level)
  fn <- paste0(path, "_nodes.tsv"); fe <- paste0(path, "_edges.tsv")
  utils::write.table(t$nodes, fn, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  utils::write.table(t$edges, fe, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(c(nodes = fn, edges = fe))
}
