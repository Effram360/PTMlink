# Source this file after library(ptmlink), then call plot_pair_cohorts(ptm_pair(...)).
# Uses the estimates and intervals returned by the reader; does not refit the data. No title unless `main` is given.
plot_pair_cohorts <- function(pair, main = NULL) {
  co <- pair$cohorts
  if (is.null(pair$row) || is.null(co) || !nrow(co)) stop("This pair has no reported cohort estimates.")
  stopifnot(all(c("cohort", "rho", "n", "lo", "hi") %in% names(co)))
  co <- co[order(co$rho), , drop = FALSE]
  old <- par(mar = c(4, 11, if (is.null(main)) 1 else 3, 2) + 0.1)
  on.exit(par(old))
  y <- seq_len(nrow(co))
  plot(co$rho, y, xlim = c(-1, 1), ylim = c(0.5, nrow(co) + 0.5), yaxt = "n", ylab = "",
       xlab = "Partial Spearman correlation (95% interval)", main = main,
       pch = 19, col = "#0072B2", cex = 1.2, bty = "n")
  abline(v = 0, lty = 2, col = "grey70")
  segments(co$lo, y, co$hi, y, col = "#0072B2", lwd = 2)
  # the article's display names (GC, ccRCC, Rare kidney cancers); `cohort` keeps the analysis key
  lab <- if ("label" %in% names(co)) ifelse(is.na(co$label), co$cohort, co$label) else co$cohort
  axis(2, at = y, labels = paste0(lab, "  (n = ", co$n, ")"), las = 1, tick = FALSE)
  invisible(co)
}
