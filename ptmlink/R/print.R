# How results print at the console. Every function still returns the full object (a data.frame of every column, or a
# list); only print() is short: a header that says what was asked and what the view hid, the key columns of the first
# rows, and a footer that says how to see the rest. as.data.frame(x), View(x) and x$column reach everything. Each
# summary fits the console width (getOption("width")): a table that would wrap drops its least needed columns first.

.ptm_fmt_num <- function(x, digits = 2) ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "-")
.ptm_fmt_q <- function(q) ifelse(is.finite(q), ifelse(q < 0.001, formatC(q, format = "e", digits = 1), formatC(q, format = "f", digits = 3)), "-")
.ptm_fmt_int <- function(x) ifelse(is.finite(x), formatC(x, format = "d", big.mark = ","), "-")
.ptm_fmt_ci <- function(lo, hi) ifelse(is.finite(lo) & is.finite(hi), paste0("[", .ptm_fmt_num(lo), ", ", .ptm_fmt_num(hi), "]"), "-")
# The pooled 95 % interval, as the website and the article give it: +-1.96 standard errors (sum of weights)^-1/2 on
# the Fisher scale, transformed back.
.ptm_pooled_ci <- function(rho, sum_weight) {
  se <- ifelse(is.finite(sum_weight) & sum_weight > 0, 1 / sqrt(sum_weight), NA_real_)
  z <- atanh(pmin(pmax(rho, -0.999999), 0.999999))
  list(lo = tanh(z - 1.959964 * se), hi = tanh(z + 1.959964 * se))
}
.ptm_support <- function(x) {
  col <- function(nm) if (nm %in% names(x)) x[[nm]] %in% TRUE else rep(FALSE, nrow(x))
  lab <- cbind(ifelse(col("anch_dose"), "dose", ""), ifelse(col("anch_gene_editing"), "editing", ""), ifelse(col("anch_enzymatic"), "sialidase", ""))
  out <- apply(lab, 1, function(r) paste(r[nzchar(r)], collapse = "+"))
  out[!nzchar(out)] <- "-"
  out
}
.ptm_release_label <- function() {
  b <- tryCatch(ptm_release(), error = function(e) NULL)
  if (is.null(b) || is.null(b$release)) "PTMlink" else paste("PTMlink", b$release)
}
.ptm_width <- function() max(60L, as.integer(getOption("width", 80L)))
# A result cut down to columns the summary does not use (x[, c("a", "b")] keeps the class) prints as the plain
# data.frame it now is.
.ptm_print_n <- function(n) {
  if (!(is.numeric(n) && length(n) == 1L && !is.na(n) && n >= 1)) stop("n is the number of rows to show: a whole number, 1 or more", call. = FALSE)
  if (is.finite(n)) as.integer(n) else n
}
.ptm_print_plain <- function(x, need, n = Inf, hide = character(0)) {
  if (all(need %in% names(x))) return(FALSE)
  # a long table cut down to other columns shows its first rows, like the full summary, never thousands of lines
  d <- as.data.frame(x); shown <- min(n, nrow(d))
  hidden <- intersect(hide, names(d)); d <- d[, setdiff(names(d), hidden), drop = FALSE]
  print(d[seq_len(shown), , drop = FALSE], row.names = FALSE)
  if (nrow(d) > shown) cat(.ptm_fmt_int(nrow(d) - shown), " more rows: print(x, n = 50). Every row: as.data.frame(x).\n", sep = "")
  if (length(hidden)) cat("Also in the table: ", paste(hidden, collapse = ", "), ".\n", sep = "")
  TRUE
}
.ptm_say <- function(text, indent = 0L) cat(strwrap(text, width = .ptm_width() - 1L, indent = indent, exdent = indent + 2L), sep = "\n")
# Lines built from short units (a hint, a statistic), wrapped between units and never inside one; `sep` ends every
# unit but the last.
.ptm_pack <- function(units, sep = "", indent = 0L, exdent = indent + 2L) {
  units <- units[nzchar(units)]; n <- length(units)
  if (!n) return(invisible())
  if (n > 1) units[-n] <- paste0(units[-n], sep)
  w <- .ptm_width() - 1L; lines <- character(0); cur <- paste0(strrep(" ", indent), units[1])
  for (u in units[-1]) {
    if (nchar(cur) + 1L + nchar(u) > w) { lines <- c(lines, cur); cur <- paste0(strrep(" ", exdent), u) } else cur <- paste(cur, u)
  }
  cat(c(lines, cur), sep = "\n")
}
# The next steps under a summary: one line each, a short label and the code to type, the code in one column.
.ptm_hints <- function(h, lead = NULL) {
  h <- h[!is.na(h) & nzchar(h)]
  for (l in lead[nzchar(lead)]) cat(strwrap(l, width = .ptm_width() - 1L), sep = "\n")   # each note on its own line
  if (!length(h)) return(invisible())
  w <- max(nchar(names(h))) + 1L
  cat(paste0("  ", formatC(paste0(names(h), ":"), width = -w), "  ", h), sep = "\n")
}
.ptm_showing <- function(shown, total) if (total > shown) paste0("Showing ", .ptm_fmt_int(shown), " of ", .ptm_fmt_int(total), ".")
.ptm_cap <- function(w) paste0(toupper(substr(w, 1L, 1L)), substring(w, 2L))
.ptm_plural <- function(n, one, many = paste0(one, "s")) paste(.ptm_fmt_int(n), if (identical(as.numeric(n), 1)) one else many)
.ptm_withheld_line <- function(w, min_k = 2L) {
  if (is.null(w)) return(NULL)
  if (is.matrix(w)) w <- colSums(w)
  if (is.null(min_k) || length(min_k) != 1L || is.na(min_k)) min_k <- 2L
  words <- c(cis = "same-protein", k = if (min_k == 2L) "single-cohort" else paste0("in fewer than ", min_k, " cohorts"), composition = "tissue composition", anchored = "without perturbation support",
             significant = "not significant", heterogeneous = "heterogeneous", loco = "below the leave-one-cohort-out floor",
             adjusted = "adjusted interval includes zero", shadow = "parent proteins co-vary", self = "site follows its own protein")
  w <- w[w > 0 & names(w) %in% names(words)]
  if (!length(w)) return(NULL)
  u <- paste(.ptm_fmt_int(w), unname(words[names(w)])); u[1] <- paste("Hidden by filters:", u[1]); u[length(u)] <- paste0(u[length(u)], ".")
  u
}
# cells: named list of character columns; right: right-align each column; drop: column names to leave out, in order,
# while the table is wider than the console
.ptm_table <- function(cells, right, drop = character(0)) {
  width_of <- function(cl) pmax(nchar(names(cl)), vapply(cl, function(v) max(nchar(v), 0L), integer(1)))
  names(right) <- names(cells)
  for (d in drop) {
    if (2L + sum(width_of(cells)) + 2L * (length(cells) - 1L) <= .ptm_width() - 1L) break
    cells[[d]] <- NULL
  }
  right <- right[names(cells)]; w <- width_of(cells)
  pad <- function(v, i) if (right[i]) formatC(v, width = w[i]) else formatC(v, width = -w[i])
  head <- paste(vapply(seq_along(cells), function(i) pad(names(cells)[i], i), ""), collapse = "  ")
  body <- vapply(seq_len(length(cells[[1]])), function(r) paste(vapply(seq_along(cells), function(i) pad(cells[[i]][r], i), ""), collapse = "  "), "")
  cat(sub(" +$", "", paste0("  ", c(head, body))), sep = "\n")
}

#' Printing PTMlink results
#'
#' Results print as a short summary that fits the console: what was asked, what the default view hid, the key
#' columns of the first rows, and how to see the rest. The objects themselves are unchanged: every column stays in
#' the data.frame (`as.data.frame(x)`, `View(x)`, `x$evidence`), and `ptm_columns()` defines each column.
#' @param x a result of [ptm_gene()], [ptm_site()], [ptm_query()], [ptm_pair()], [ptm_explain()], [ptm_status()],
#'   [ptm_network()], [ptm_list()] or one of its tables, [ptm_release()], [ptm_cohorts()], [ptm_sites()],
#'   [ptm_columns()], [ptm_glycans()], [ptm_default_filters()], [ptm_cite()], [ptm_reference()] or [ptm_tools()]
#' @param n rows to show; for several genes or sites asked at once, the top rows of each are shown
#' @param ... ignored
#' @return `x`, invisibly
#' @name ptm_print
NULL

#' @rdname ptm_print
#' @export
print.ptm_pairs <- function(x, n = 10L, ...) {
  n <- .ptm_print_n(n)
  # the two site identifiers are the two halves of pair_id: a cut-down table shows pair_id alone
  if (.ptm_print_plain(x, c("glyco_site", "phospho_site", "rho_combined", "evidence", "q_value", "k_cohorts"), n,
                       hide = if ("pair_id" %in% names(x)) c("glyco_row_id", "phospho_row_id") else character(0))) return(invisible(x))
  what <- attr(x, "query"); n_query <- attr(x, "n_query"); n_view <- attr(x, "n_view"); by <- attr(x, "by")
  qcol <- intersect(c("query_gene", "query_site"), names(x))[1]
  groups <- if (!is.na(qcol)) unique(x[[qcol]]) else character(0)
  stacked <- length(groups) > 1L
  # the header states only what is true of the rows in hand: a subset says "3 of 8,579", a re-sorted table does not
  # claim the ranking
  desc <- function(e) length(e) < 2L || !is.unsorted(rev(e[is.finite(e)]))
  order_note <- if (nrow(x) < 2L) "" else if (desc(x$evidence)) ", sorted by evidence" else if (stacked && all(vapply(groups, function(g) desc(x$evidence[x[[qcol]] == g]), logical(1))))
    paste0(", sorted by evidence within each ", if (is.null(by)) "query" else by) else ""
  count <- if (!is.null(n_query) && n_query > nrow(x)) paste0(.ptm_fmt_int(nrow(x)), " of ", .ptm_plural(n_query, "pair"))
           else if (!is.null(n_view) && n_view > nrow(x)) paste0("the top ", .ptm_fmt_int(nrow(x)), " of ", .ptm_plural(n_view, "pair"))
           else .ptm_plural(nrow(x), "pair")
  .ptm_pack(c(paste0(if (!is.null(what)) paste0(.ptm_cap(what), ": "), count, order_note), paste0("(", .ptm_release_label(), ")")))   # never broken inside the release label
  line <- .ptm_withheld_line(attr(x, "withheld"), attr(x, "filters")$min_k)
  if (!is.null(line)) .ptm_pack(line, ",")
  if (!nrow(x)) {
    hidden <- attr(x, "withheld"); matched <- attr(x, "n_matched")
    cat(if (!is.null(n_view) && n_view > 0) "No rows returned: the limit was 0; a higher limit returns them.\n"
        else if (!is.null(matched) && matched == 0) "No reported pair meets the conditions.\n"
        else if (!is.null(hidden) && sum(hidden) > 0) "No pair passes the filters. To relax them: ?ptm_default_filters\n"
        else "No reported pair.\n")
    return(invisible(x))
  }
  # several genes or sites asked at once: the top rows of each, so no query is pushed off the screen by the first
  per <- if (stacked) max(1L, n %/% length(groups)) else n
  show <- if (stacked) utils::head(unlist(lapply(groups, function(g) utils::head(which(x[[qcol]] == g), per))), max(n, length(groups))) else seq_len(min(n, nrow(x)))
  s <- x[show, , drop = FALSE]
  # a field whose columns were not asked for (ptm_query(columns = ...)) is left out, never shown as "-"
  has_ci <- "sum_weight" %in% names(s); has_support <- any(c("anch_dose", "anch_gene_editing", "anch_enzymatic") %in% names(s))
  ci <- if (has_ci) .ptm_pooled_ci(s$rho_combined, s$sum_weight) else NULL
  cells <- list(glyco_site = s$glyco_site, phospho_site = s$phospho_site, rho = .ptm_fmt_num(s$rho_combined),
                `95% CI` = if (has_ci) .ptm_fmt_ci(ci$lo, ci$hi), evidence = .ptm_fmt_num(s$evidence, 1), q = .ptm_fmt_q(s$q_value),
                cohorts = .ptm_fmt_int(s$k_cohorts), perturbation = if (has_support) sub("^-$", "", .ptm_support(s)))
  right <- c(glyco_site = FALSE, phospho_site = FALSE, rho = TRUE, `95% CI` = TRUE, evidence = TRUE, q = TRUE, cohorts = TRUE, perturbation = FALSE)
  keep <- !vapply(cells, is.null, logical(1)); cells <- cells[keep]; right <- unname(right[keep])
  if (!is.na(qcol)) { cells <- c(list(query = s[[qcol]]), cells); right <- c(FALSE, right) }
  cat("\n")
  .ptm_table(cells, right, drop = c("95% CI", "q"))
  cat("\n")
  .ptm_hints(c(`more rows` = if (nrow(x) > nrow(s)) "print(x, n = 50)", `all columns` = "View(x)",
               `one pair` = paste0("ptm_pair(\"", s$glyco_site[1], "\", \"", s$phospho_site[1], "\")"),
               `how it was read` = if (!is.null(attr(x, "provenance"))) "attr(x, \"provenance\")", filters = "?ptm_default_filters"),
             lead = if (stacked && nrow(x) > nrow(s)) paste0("Showing the top ", per, " of each ", if (is.null(by)) "query" else by, ": ", .ptm_fmt_int(nrow(s)), " of ", .ptm_fmt_int(nrow(x)), ".")
                    else .ptm_showing(nrow(s), nrow(x)))
  invisible(x)
}

#' @rdname ptm_print
#' @export
print.ptm_status <- function(x, ...) {
  if (.ptm_print_plain(x, c("query", "gene", "resolved_via", "status_glyco", "status_phospho"))) return(invisible(x))
  for (i in seq_len(nrow(x))) {
    r <- x[i, , drop = FALSE]
    if (is.na(r$gene)) { cat(r$query, ": not a recognized symbol, alias or UniProt accession\n", sep = ""); next }
    asked <- !identical(toupper(r$query), toupper(r$gene))
    how <- sub(" of .*$", "", r$resolved_via)
    cat(r$gene, if (asked) paste0(" (asked as ", r$query, if (!is.na(how) && how != "current symbol") paste0(", ", how), ")"), "\n", sep = "")
    side <- function(lab, status, ncoh, nrep) {
      extra <- if ("n_cohorts_glyco_detected" %in% names(r)) paste0("; detected in ", .ptm_fmt_int(ncoh), if (identical(as.numeric(ncoh), 1)) " cohort, " else " cohorts, ",
                                                                   .ptm_fmt_int(nrep), if (identical(as.numeric(nrep), 1)) " reported pair" else " reported pairs") else ""
      cat("  ", formatC(lab, width = -14), status, extra, "\n", sep = "")
    }
    side("glycosites:", r$status_glyco, r$n_cohorts_glyco_detected, r$n_served_as_glyco)
    side("phosphosites:", r$status_phospho, r$n_cohorts_phospho_detected, r$n_served_as_phospho)
  }
  invisible(x)
}

#' @rdname ptm_print
#' @export
print.ptm_pair <- function(x, ...) {
  if (is.null(x$row)) {
    if (!is.null(x$missing) && nrow(x$missing)) { u <- paste(x$missing$glyco_site, x$missing$phospho_site, sep = " x "); u[1] <- paste("Not reported:", u[1]); u[length(u)] <- paste0(u[length(u)], "."); .ptm_pack(u, ","); .ptm_say(paste0(if (nrow(x$missing) == 1L) "It is " else "They are ", .PTM_ABSENT_PAIR, ".")) }
    else .ptm_say(paste0("Not a reported pair: ", x$note))
    return(invisible(x))
  }
  if (nrow(x$row) > 1) {
    rows <- .ptm_as_pairs(x$row, "the pairs asked")   # printed in the order asked; the header claims no ranking unless they are in it
    print(rows, n = nrow(rows))
    if (!is.null(x$missing) && nrow(x$missing)) { u <- paste(x$missing$glyco_site, x$missing$phospho_site, sep = " x "); u[1] <- paste("Not reported:", u[1]); .ptm_pack(u, ",") }
    return(invisible(x))
  }
  r <- x$row
  ci <- .ptm_pooled_ci(r$rho_combined, if ("sum_weight" %in% names(r)) r$sum_weight else NA)
  sup <- .ptm_support(r)
  words <- c(dose = "tunicamycin dose response", editing = "glycogene editing", sialidase = "sialidase")
  has <- function(nm) nm %in% names(r)
  kv <- c(`pair id` = if (has("pair_id") && !is.na(r$pair_id)) r$pair_id,
          `pooled rho` = paste0(.ptm_fmt_num(r$rho_combined), ", 95% CI ", .ptm_fmt_ci(ci$lo, ci$hi)),
          evidence = paste0(.ptm_fmt_num(r$evidence, 1), ", q = ", .ptm_fmt_q(r$q_value)),
          cohorts = paste0(r$k_cohorts, ", ", .ptm_fmt_int(x$n_total), " tumors", if (is.finite(r$I2)) paste0("; I2 = ", .ptm_fmt_num(r$I2))),
          perturbation = paste0(if (sup == "-") "none" else paste(words[strsplit(sup, "+", fixed = TRUE)[[1]]], collapse = ", "),
                                if (has("anch_pair_observed") && isTRUE(r$anch_pair_observed)) " (the pair itself was measured)"),
          `parent proteins` = if (has("rho_prot")) paste("rho", .ptm_fmt_num(r$rho_prot)),
          `own-protein coupling` = if (has("self_coupling")) .ptm_fmt_num(r$self_coupling),
          `adjusted 95% CI` = if (has("rho_adjusted_lo")) paste(.ptm_fmt_ci(r$rho_adjusted_lo, r$rho_adjusted_hi), "(an annotation; never used to rank)"))
  cat(r$glyco_site, " x ", r$phospho_site, "\n", sep = "")
  .ptm_definitions(names(kv), NULL, unname(kv), cap = 22L)
  co <- x$cohorts
  if (!is.null(co) && nrow(co)) {
    cat("\n")
    .ptm_table(list(cohort = co$label, rho = .ptm_fmt_num(co$rho), `95% CI` = .ptm_fmt_ci(co$lo, co$hi), tumors = .ptm_fmt_int(co$n)),
               right = c(FALSE, TRUE, TRUE, TRUE))
  }
  cat("\n")
  .ptm_hints(c(cautions = paste0("ptm_explain(\"", r$glyco_site, "\", \"", r$phospho_site, "\")"), glycoforms = "x$glycoforms",
               `cohort table` = "x$cohorts", `all columns` = "x$row"))
  invisible(x)
}

# a rule's name in plain words, for the "Does not apply" line
.ptm_rule_words <- function(rule) { w <- gsub("_", " ", rule); w[rule == "protein_shadow"] <- "parent-protein co-variation"; w }
# One short line per reading rule, for the console; the full sentence stays in x$statement. `extra` is the part of
# the statement that carries this pair's values (the statement after the rule's fixed sentence).
.ptm_rule_short <- function(rule, extra, statement) {
  lead <- c(direction = "Co-variation only: neither site is shown to regulate the other.",
            significance = "Significant after Benjamini-Hochberg correction:",
            perturbation = "Perturbation support, which admits a pair but never raises its rank.",
            same_protein = "Both sites are on one protein; read the pair as shared abundance first.",
            protein_shadow = "The parent proteins co-vary, so the sites may simply follow them.",
            self_coupling = "One site follows its own protein's abundance.",
            tissue_composition = "Both sites track tissue composition, which alone can make them co-vary.",
            heterogeneity = "The cohorts disagree beyond sampling noise.",
            one_cohort = "Measured in one cohort only.",
            adjusted = "The program-adjusted 95% interval excludes zero (an annotation, never a ranking).",
            glycoform_pooling = "Different glycoforms are pooled across cohorts; see ptm_pair()$glycoforms.",
            sequon_ambiguity = "The peptide carried more than one sequon; the site label is an attribution.")
  if (!rule %in% names(lead)) return(if (nzchar(extra)) extra else statement)
  trimws(paste(lead[[rule]], extra))
}

#' @rdname ptm_print
#' @export
print.ptm_explain <- function(x, ...) {
  if (.ptm_print_plain(x, c("glyco_site", "phospho_site", "rule", "applies", "statement"))) return(invisible(x))
  base <- vapply(.ptm_rules(), function(r) r$statement, ""); names(base) <- vapply(.ptm_rules(), function(r) r$id, "")
  extra_of <- function(e) vapply(seq_len(nrow(e)), function(i) {
    b <- base[e$rule[i]]
    if (!is.na(b) && startsWith(e$statement[i], b)) trimws(substring(e$statement[i], nchar(b) + 1L)) else ""
  }, "")
  for (key in unique(paste(x$glyco_site, x$phospho_site, sep = " x "))) {
    e <- x[paste(x$glyco_site, x$phospho_site, sep = " x ") == key, , drop = FALSE]
    e$extra <- extra_of(e)
    cat(key, "\n", sep = "")
    blank <- is.na(e$applies) & startsWith(e$statement, "Not assessed")   # a tested rule whose value this pair lacks
    yes <- e[e$applies %in% TRUE, , drop = FALSE]; info <- e[is.na(e$applies) & !blank, , drop = FALSE]; no <- e[e$applies %in% FALSE, , drop = FALSE]
    bullet <- function(d) for (i in seq_len(nrow(d))) cat(strwrap(.ptm_rule_short(d$rule[i], d$extra[i], d$statement[i]), width = .ptm_width() - 1L, initial = "  - ", prefix = "    "), sep = "\n")
    if (nrow(yes)) { cat("Applies:\n"); bullet(yes) }
    if (nrow(info)) { cat("Context:\n"); bullet(info) }
    if (nrow(no)) { u <- .ptm_rule_words(no$rule); u[1] <- paste("Does not apply:", u[1]); u[length(u)] <- paste0(u[length(u)], "."); .ptm_pack(u, ",") }
    if (any(blank)) { u <- .ptm_rule_words(e$rule[blank]); u[1] <- paste("Not assessed, no value:", u[1]); u[length(u)] <- paste0(u[length(u)], "."); .ptm_pack(u, ",") }
  }
  cat("\n"); .ptm_hints(c(`full wording` = "x$statement"))
  invisible(x)
}

#' @rdname ptm_print
#' @export
print.ptm_network <- function(x, n = 10L, ...) {
  n <- .ptm_print_n(n)
  if (.ptm_print_plain(x, c("partner", "n_site_pairs", "best_glyco_site", "best_phospho_site", "rho_at_best", "evidence_best"), n)) return(invisible(x))
  # the header says what is true of the rows in hand: stacked tables name every gene, a subset counts its own pairs
  gene <- if ("gene" %in% names(x) && nrow(x)) unique(x$gene) else attr(x, "gene")
  cat("Network of ", paste(gene, collapse = ", "), ": ", .ptm_plural(nrow(x), "partner link"), " from ", .ptm_plural(if (nrow(x)) sum(x$n_site_pairs) else 0L, "site pair"),
      " (", .ptm_release_label(), ")\n", sep = "")
  cat("One row per partner, with its strongest site pair; links have no direction.\n")
  line <- .ptm_withheld_line(attr(x, "withheld"), attr(x, "filters")$min_k)
  if (!is.null(line)) .ptm_pack(line, ",")
  if (!nrow(x)) return(invisible(x))
  s <- x[seq_len(min(n, nrow(x))), , drop = FALSE]
  # the query protein's site and the partner's, each in its own column; with one query gene its name heads the column
  # and the gene prefix is dropped from both sites
  mine <- if ("role" %in% names(s)) s$role == "glyco" else rep(TRUE, nrow(s))
  own <- ifelse(mine, s$best_glyco_site, s$best_phospho_site); other <- ifelse(mine, s$best_phospho_site, s$best_glyco_site)
  one <- length(gene) == 1L
  if (one) { own <- sub("^.*_", "", own); other <- sub("^.*_", "", other) }
  cells <- list(partner = s$partner, `site pairs` = .ptm_fmt_int(s$n_site_pairs), own, other,
                rho = .ptm_fmt_num(s$rho_at_best), evidence = .ptm_fmt_num(s$evidence_best, 1), perturbation = ifelse(s$anch_any %in% TRUE, "yes", ""))
  names(cells)[3:4] <- c(if (one) paste(gene, "site") else "site", "partner site")
  cat("\n")
  .ptm_table(cells, right = c(FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE), drop = c("perturbation", "site pairs"))
  cat("\n")
  .ptm_hints(c(`more rows` = if (nrow(x) > nrow(s)) "print(x, n = 50)", `all columns` = "View(x)", plot = "ptm_plot_network(x)"),
             lead = .ptm_showing(nrow(s), nrow(x)))
  invisible(x)
}

#' @rdname ptm_print
#' @export
print.ptm_release <- function(x, ...) {
  cat("PTMlink release ", x$release, if (!is.null(x$released) && !is.na(x$released)) paste0(", released ", x$released),
      if (!is.null(x$built_utc)) paste0(" (data built ", substr(x$built_utc, 1, 10), ")"), "\n", sep = "")
  has <- function(nm) !is.null(x[[nm]]) && length(x[[nm]]) == 1L && is.finite(x[[nm]])
  inf <- x$inference
  kv <- c(`reported pairs` = if (has("n_served_rows")) paste0(.ptm_fmt_int(x$n_served_rows), if (has("n_genes")) paste0(" on ", .ptm_plural(x$n_genes, "gene")),
                                                              if (has("n_served_anchored")) paste0("; ", .ptm_fmt_int(x$n_served_anchored), " with perturbation support")),
          cohorts = if (has("n_cohorts")) paste0(.ptm_plural(x$n_cohorts, "tumor cohort"), if (has("n_samples")) paste0(", ", .ptm_plural(x$n_samples, "tumor")),
                                                 if (has("n_normal_cohorts") && x$n_normal_cohorts > 0) paste0("; ", .ptm_plural(x$n_normal_cohorts, "normal adjacent tissue cohort"), ", never pooled")),
          statistics = if (!is.null(inf)) if (identical(inf$inference_revision, "student_normal_bh")) "Student tests pooled by weighted Stouffer" else "Fisher-normal pooling",
          `q values` = if (!is.null(inf)) "Benjamini-Hochberg over all scored site pairs")
  if (length(kv)) .ptm_definitions(names(kv), NULL, unname(kv))
  cat("\n"); .ptm_hints(c(`all fields` = "str(x)", cohorts = "ptm_cohorts()", cite = "ptm_cite()"))
  invisible(x)
}

# ---- reference tables: the cohorts, sites, columns, glycans, filters, citation, reference and tool definitions ----

# a two- or three-column listing whose last column wraps under itself, for definitions that do not fit one line
# A key longer than `cap` sits alone on its line, so one long name does not squeeze every definition.
.ptm_definitions <- function(keys, mids = NULL, text, cap = 20L) {
  wk <- min(max(nchar(keys)), cap); wm <- if (is.null(mids)) 0L else max(nchar(mids))
  lead <- 2L + wk + 2L + if (wm) wm + 2L else 0L
  avail <- max(20L, .ptm_width() - 1L - lead)
  for (i in seq_along(keys)) {
    body <- strwrap(text[i], width = avail); if (!length(body)) body <- ""
    long <- nchar(keys[i]) > wk
    if (long) cat("  ", keys[i], "\n", sep = "")
    first <- paste0("  ", formatC(if (long) "" else keys[i], width = -wk), "  ", if (wm) paste0(formatC(mids[i], width = -wm), "  "))
    cat(sub(" +$", "", c(paste0(first, body[1]), if (length(body) > 1) paste0(strrep(" ", lead), body[-1]))), sep = "\n")
  }
}

#' @rdname ptm_print
#' @export
print.ptm_cohorts <- function(x, ...) {
  if (.ptm_print_plain(x, c("cohort", "label", "role"))) return(invisible(x))
  tum <- x$role %in% "tumour"
  cat(.ptm_cap(.ptm_plural(nrow(x), "cohort")), " (", .ptm_release_label(), ")\n", sep = "")
  .ptm_say(paste0(sum(tum), " tumor cohorts are pooled", if (any(!tum)) paste0("; ", sum(!tum), " normal adjacent tissue cohorts are shown beside them and never pooled"), "."))
  cat("\n")
  col <- function(nm) if (nm %in% names(x)) x[[nm]] else rep(NA, nrow(x))
  tissue <- sub(",\\s*normal adjacent tissue.*$", "", col("tumour_type"))   # the role column already says normal tissue
  .ptm_table(list(cohort = x$label, `tumor type` = ifelse(is.na(tissue), "-", tissue),
                  samples = .ptm_fmt_int(col("n_samples_glyco_phospho")), role = ifelse(tum, "tumor", "normal tissue"),
                  covariates = ifelse(is.na(col("covariates")), "-", col("covariates"))),
             right = c(FALSE, FALSE, TRUE, FALSE, FALSE), drop = c("covariates", "tumor type"))
  cat("\n")
  .ptm_hints(c(`sources, enrichment` = "View(x)"), lead = "Samples: those with both a glycoproteome and a phosphoproteome.")
  invisible(x)
}

#' @rdname ptm_print
#' @export
print.ptm_sites <- function(x, n = 30L, ...) {
  n <- .ptm_print_n(n)
  if (.ptm_print_plain(x, c("gene", "side", "site"), n)) return(invisible(x))
  genes <- unique(x$gene)
  acc <- if ("uniprot_acc" %in% names(x)) x$uniprot_acc[match(genes, x$gene)] else rep(NA, length(genes))
  who <- paste(ifelse(is.na(acc), genes, paste0(genes, " (", acc, ")")), collapse = ", ")
  ng <- sum(x$side == "glyco"); np <- sum(x$side == "phospho")
  .ptm_say(paste0(who, ": ", .ptm_plural(ng, "glycosite"), ", ", .ptm_plural(np, "phosphosite"), " reported (", .ptm_release_label(), ")"))
  s <- x[seq_len(min(n, nrow(x))), , drop = FALSE]
  col <- function(nm) if (nm %in% names(s)) s[[nm]] else rep(NA, nrow(s))
  cells <- list(site = s$site, side = s$side, cohorts = .ptm_fmt_int(col("n_cohorts_served")),
                `own-protein coupling` = .ptm_fmt_num(as.numeric(col("self_coupling_median"))),
                pLDDT = ifelse(is.finite(as.numeric(col("plddt"))), formatC(as.numeric(col("plddt")), format = "f", digits = 0), "-"),
                model = ifelse(col("residue_matches_model") %in% TRUE, "matches", ifelse(col("residue_matches_model") %in% FALSE, "differs", "-")))
  cat("\n")
  .ptm_table(cells, right = c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE), drop = c("model", "own-protein coupling"))
  cat("\n")
  .ptm_hints(c(`more rows` = if (nrow(x) > nrow(s)) "print(x, n = 100)", `coordinates, all columns` = "View(x)"),
             lead = c("pLDDT: the AlphaFold model's confidence at the site.", .ptm_showing(nrow(s), nrow(x))))
  invisible(x)
}

#' @rdname ptm_print
#' @export
print.ptm_columns <- function(x, ...) {
  if (.ptm_print_plain(x, c("column", "type", "meaning"))) return(invisible(x))
  only_full <- if ("in_shards" %in% names(x)) !(x$in_shards %in% TRUE) else rep(FALSE, nrow(x))
  cat(nrow(x), " columns of a reported pair\n\n", sep = "")
  .ptm_definitions(paste0(x$column, ifelse(only_full, "*", "")), ifelse(is.na(x$type), "-", x$type), x$meaning)
  cat("\n")
  .ptm_hints(c(`one column` = "x[x$column == \"evidence\", ]"), lead = if (any(only_full)) "* only in the full release table, read by ptm_query() and ptm_list().")
  invisible(x)
}

#' @rdname ptm_print
#' @export
print.ptm_glycans <- function(x, n = 10L, ...) {
  n <- .ptm_print_n(n)
  if (.ptm_print_plain(x, c("composition", "glycan_class"), n)) return(invisible(x))
  cat(.ptm_fmt_int(nrow(x)), " glycan compositions, the most often measured first (", .ptm_release_label(), ")\n\n", sep = "")
  s <- x[seq_len(min(n, nrow(x))), , drop = FALSE]
  col <- function(nm) if (nm %in% names(s)) s[[nm]] else rep(NA, nrow(s))
  .ptm_table(list(composition = s$composition, class = s$glycan_class, notation = ifelse(is.na(col("community_notation")), "-", col("community_notation")),
                  `mass (Da)` = .ptm_fmt_num(as.numeric(col("monoisotopic_mass")), 1), cohorts = .ptm_fmt_int(col("n_cohorts_observed"))),
             right = c(FALSE, FALSE, FALSE, TRUE, TRUE), drop = c("mass (Da)", "notation"))
  cat("\n")
  .ptm_hints(c(`more rows` = if (nrow(x) > nrow(s)) "print(x, n = 50)", `all columns` = "View(x)"),
             lead = c("Composition: N HexNAc, H hexose, F fucose, S NeuAc, G NeuGc.", .ptm_showing(nrow(s), nrow(x))))
  invisible(x)
}

#' @rdname ptm_print
#' @export
print.ptm_filters <- function(x, ...) {
  def <- tryCatch(.ptm_reference_raw()$filters, error = function(e) list())
  meaning <- setNames(vapply(def, function(f) f$meaning, ""), vapply(def, function(f) f$name, ""))
  val <- vapply(x, function(v) if (length(v) != 1L || is.na(v)) "none" else as.character(v), "")
  cat("PTMlink view filters\n\n")
  .ptm_definitions(names(x), val, ifelse(names(x) %in% names(meaning), meaning[names(x)], ""))
  cat("\n")
  cat("To change one, pass a modified copy:\n  f <- modifyList(ptm_default_filters(), list(min_k = 1))\n  ptm_gene(\"EGFR\", \"glyco\", filters = f)\n")
  invisible(x)
}

#' @rdname ptm_print
#' @export
print.ptm_cite <- function(x, ...) {
  cat("To cite PTMlink:\n\n")
  cat(strwrap(x$text, width = .ptm_width() - 1L, indent = 2L, exdent = 2L), sep = "\n")
  cat("\n")
  .ptm_hints(c(BibTeX = "toBibtex(citation(\"ptmlink\"))"), lead = if (!is.null(x$release) && !is.na(x$release)) paste0("Name the data release you used: ", x$release, "."))
  invisible(x)
}

#' @rdname ptm_print
#' @export
print.ptm_reference <- function(x, ...) {
  n <- function(v) length(v)
  cat("PTMlink reference for agents and tool loaders (ptmlink ", if (!is.null(x$version)) x$version else "", ")\n", sep = "")
  cat("  ", n(x$rules), " reading rules, ", n(x$functions), " functions, ", n(x$columns), " columns, ", n(x$filters), " filters\n\n", sep = "")
  rules <- x$rules
  .ptm_definitions(vapply(rules, function(r) r$id, ""), NULL, vapply(rules, function(r) r$statement, ""))
  cat("\n")
  .ptm_hints(c(parts = "x$rules, x$functions, x$columns, x$filters", JSON = "ptm_reference(as = \"json\")"))
  invisible(x)
}

#' @rdname ptm_print
#' @export
print.ptm_tools <- function(x, ...) {
  nm <- vapply(x, function(t) if (!is.null(t$name)) t$name else if (!is.null(t[["function"]]$name)) t[["function"]]$name else "", "")
  cat(length(x), " tool definitions, ", attr(x, "format"), " format\n\n", sep = "")
  .ptm_pack(paste0(nm, c(rep(",", length(nm) - 1L), "")), indent = 2L, exdent = 2L)
  cat("\n")
  .ptm_hints(c(`one tool` = "x[[1]]", `assistant instructions` = "attr(x, \"instructions\")", JSON = "ptm_tools(as = \"json\")"))
  invisible(x)
}

# ---- a list of proteins or sites ----------------------------------------------------------------------------

#' @rdname ptm_print
#' @export
print.ptm_list <- function(x, n = 5L, ...) {
  n <- .ptm_print_n(n)
  m <- x$members; cn <- x$connectivity; site <- identical(x$level, "site")
  unit <- if (site) "site" else "protein"
  n_units <- if (site) length(x$sites) else length(unique(m$gene[!is.na(m$gene)]))
  cat("List of ", .ptm_plural(n_units, unit), " against the reported pairs (", .ptm_release_label(), ")\n", sep = "")
  .ptm_say(paste0("Pairs in view: ", .ptm_fmt_int(nrow(x$internal)), " within the list, ", .ptm_fmt_int(nrow(x$neighbourhood)), " with a partner outside it",
                  if (nrow(x$same_protein)) paste0(", ", .ptm_fmt_int(nrow(x$same_protein)), " on one protein"), "."))
  line <- .ptm_withheld_line(x$withheld, x$provenance$filters$min_k)
  if (!is.null(line)) .ptm_pack(line, ",")
  if (isTRUE(x$mixed)) .ptm_say(paste0("Read at site level: ", .ptm_plural(sum(m$kind == "gene" & !is.na(m$gene)), "gene entry", "gene entries"), " expanded to reported sites."))
  lost <- m$query[is.na(m$gene)]
  if (length(lost)) { u <- lost; u[1] <- paste("Not recognized:", u[1]); u[length(u)] <- paste0(u[length(u)], "."); .ptm_pack(u, ",") }

  # the entries: how each resolved and how many pairs it brings
  cat("\nMembers\n")
  s <- utils::head(m, 10L)
  col <- function(nm) if (nm %in% names(s)) s[[nm]] else rep(NA, nrow(s))
  # counted by site, an entry shows its own reported pairs; counted by protein, the protein's pairs on each side
  if (site) .ptm_table(list(entry = s$query, protein = ifelse(is.na(s$gene), "-", s$gene), `reported pairs` = .ptm_fmt_int(col("n_reported")),
                            `within list` = .ptm_fmt_int(col("n_internal")), `outside list` = .ptm_fmt_int(col("n_neighbourhood"))), right = c(FALSE, FALSE, TRUE, TRUE, TRUE))
  else .ptm_table(list(entry = s$query, protein = ifelse(is.na(s$gene), "-", s$gene), `glyco pairs` = .ptm_fmt_int(col("n_served_as_glyco")),
                       `phospho pairs` = .ptm_fmt_int(col("n_served_as_phospho")), `within list` = .ptm_fmt_int(col("n_internal")), `outside list` = .ptm_fmt_int(col("n_neighbourhood"))),
                  right = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE), drop = c("phospho pairs", "glyco pairs"))
  if (nrow(m) > nrow(s)) cat("  ", .ptm_fmt_int(nrow(m) - nrow(s)), " more entries\n", sep = "")

  # do the members pair with each other more often than their own pair counts predict
  num <- function(v, d = 1) if (length(v) != 1L || !is.finite(v)) "-" else formatC(v, format = "f", digits = d, big.mark = ",")
  cat("\nConnectivity among the members\n")
  if (is.finite(cn$expected_pairs) && cn$expected_pairs == 0 && cn$internal_pairs == 0)
    .ptm_say(paste0("No pair within the list is possible: ", if (length(unique(m$gene[!is.na(m$gene)])) < 2L) "it holds one protein, and pairs within a protein are not counted."
                    else "its members with significant pairs are all glycosites, or all phosphosites."), indent = 2L)
  else .ptm_say(paste0(.ptm_plural(cn$internal_pairs, "significant pair"), " within the list; ", num(cn$expected_pairs), " expected from the members' own pair counts",
                  if (is.finite(cn$ratio)) paste0(" (", num(cn$ratio, 2), " times as many)"), ".",
                  if (is.finite(cn$p_perm)) paste0(" Permutation P = ", .ptm_fmt_q(cn$p_perm), " over ", .ptm_fmt_int(cn$n_perm), " random lists.")), indent = 2L)

  # the partners outside the list that more members reach than expected
  top <- function(tab, title, key, head) {
    cat("\n", title, "\n", sep = "")
    if (!nrow(tab)) { cat("  none\n"); return(invisible()) }
    d <- as.data.frame(tab); shared <- d$members_reached >= 2L
    # a partner reached by one member is a single pair, already in the row tables: the summary shows the shared ones,
    # and says so when there is none (a fold computed on one member is not a finding)
    if (!any(shared)) {
      .ptm_say(paste0("None is reached by two or more members; ", .ptm_fmt_int(nrow(d)), " reached by one member each (single pairs, in x$", if (key == "site") "partner_sites" else "convergent", ")."), indent = 2L)
      return(invisible())
    }
    d <- d[shared, , drop = FALSE]
    k <- utils::head(d, n)
    cells <- list(k[[key]], side = k[[if (key == "site") "side" else "partner_side"]], members = paste(.ptm_fmt_int(k$members_reached), "of", .ptm_fmt_int(k$members_testable)),
                  fold = .ptm_fmt_num(k$fold, 1), q = .ptm_fmt_q(k$q), `best evidence` = .ptm_fmt_num(k$evidence_best, 1))
    names(cells)[1] <- head
    .ptm_table(cells, right = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE), drop = c("best evidence"))
    sig <- sum(tab$q < 0.05, na.rm = TRUE)
    .ptm_say(paste0(if (nrow(d) > nrow(k)) paste0("Top ", nrow(k), " by P of ") else "", .ptm_fmt_int(nrow(d)), if (nrow(d) < nrow(tab)) paste0(" reached by two or more members (", .ptm_fmt_int(nrow(tab)), " reached in all)") else " reached",
                    "; ", .ptm_fmt_int(sig), " at q < 0.05."), indent = 2L)
  }
  top(x$partner_sites, "Partner sites reached by more members than expected", "site", "site")
  top(x$convergent, "Partner proteins reached by more members than expected", "partner", "protein")

  cat("\n")
  .ptm_say(paste0("The three tests count ", if (!is.null(x$universe$pairs)) paste0(.ptm_fmt_int(x$universe$pairs), " "), "significant pairs (q < 0.05) between two different proteins."))
  cs <- x$composition; out_share <- cs$share[cs$table == "neighbourhood"]
  if (length(out_share) && is.finite(out_share) && out_share >= 0.05)
    .ptm_say(paste0(formatC(100 * out_share, format = "f", digits = 0), "% of the pairs with a partner outside the list carry a tissue-composition label on both sites."))
  cat("\n")
  .ptm_hints(c(members = "x$members", `pairs within` = "x$internal", `pairs outside` = "x$neighbourhood", `partner sites` = "x$partner_sites",
               `partner proteins` = "x$convergent", plot = if (nrow(x$internal)) "ptm_plot_network(x)"))
  invisible(x)
}

#' @rdname ptm_print
#' @export
print.ptm_partners <- function(x, n = 10L, ...) {
  n <- .ptm_print_n(n)
  key <- if ("site" %in% names(x)) "site" else "partner"
  if (.ptm_print_plain(x, c(key, "members_reached", "members_testable", "tier_reach", "tier_n", "fold", "p", "q"), n)) return(invisible(x))
  what <- if (key == "site") "partner site" else "partner protein"
  cat(.ptm_cap(.ptm_plural(nrow(x), what)), if (nrow(x) > 1L) ", sorted by P", "; ", .ptm_fmt_int(sum(x$q < 0.05, na.rm = TRUE)), " at q < 0.05\n", sep = "")
  if (!nrow(x)) return(invisible(x))
  s <- utils::head(as.data.frame(x), n)
  cells <- list(s[[key]], side = s[[if (key == "site") "side" else "partner_side"]], members = paste(.ptm_fmt_int(s$members_reached), "of", .ptm_fmt_int(s$members_testable)),
                `all pairs` = paste(.ptm_fmt_int(s$tier_reach), "of", .ptm_fmt_int(s$tier_n)), fold = .ptm_fmt_num(s$fold, 1), P = .ptm_fmt_q(s$p), q = .ptm_fmt_q(s$q),
                `best evidence` = .ptm_fmt_num(s$evidence_best, 1))
  names(cells)[1] <- if (key == "site") "site" else "protein"
  cat("\n")
  .ptm_table(cells, right = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), drop = c("best evidence", "P", "all pairs"))
  cat("\n")
  .ptm_hints(c(`more rows` = if (nrow(x) > nrow(s)) "print(x, n = 50)", `members reached` = "x$members", `all columns` = "View(x)"),
             lead = c("members: list members that pair with it. all pairs: its reach among every significant pair.", .ptm_showing(nrow(s), nrow(x))))
  invisible(x)
}
