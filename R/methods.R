#' @title Methods for federated NoLimits fits
#' @description Standard extractors for the object returned by
#'   [ds.nolimitsFit()]. `coef()`, `vcov()` and `confint()` report the natural
#'   (untransformed) parameter scale by default; pass
#'   `scale = "transformed"` for the scale the optimiser worked on.
#' @param x,object an object of class `"ds_nolimits_fit"`.
#' @param scale `"natural"` (the default) or `"transformed"`.
#' @param parm a character vector of parameter names, or NULL for all of them.
#' @param level the confidence level.
#' @param digits the number of significant digits to print.
#' @param ... ignored.
#' @return `coef()` a named numeric vector; `vcov()` a matrix or NULL;
#'   `confint()` a two-column matrix; `summary()` an object of class
#'   `"summary.ds_nolimits_fit"`; `print()` its argument, invisibly.
#' @name ds_nolimits_fit-methods
NULL

#' @rdname ds_nolimits_fit-methods
#' @export
print.ds_nolimits_fit <- function(x, digits = 4L, ...) {
  cat("Federated NoLimits fit\n")
  cat("  estimator: ", x$estimator,
      if (identical(x$estimator, "ghq")) paste0(" (level ", x$ghq.level, ")") else "",
      if (x$estimator %in% c("mcem", "saem") && !is.null(x$outer.iters))
        paste0(" (", x$outer.iters, " outer iterations",
               if (identical(x$estimator, "saem") && !is.null(x$coordinator))
                 paste0(", coordinator ", x$coordinator) else "", ")") else "",
      "\n", sep = "")
  cat("  sites:     ", nrow(x$sites), " (",
      sum(x$sites$n.subjects), " subjects, ", sum(x$sites$n.obs), " observations)\n", sep = "")
  cat("  loglik:    ", format(x$loglik, digits = digits), "\n", sep = "")
  cat("  rounds:    ", x$rounds, "\n", sep = "")
  if (!is.null(x$model.hash)) {
    cat("  model:     ", substr(x$model.hash, 1L, 12L), " (agreed by all sites)\n",
        sep = "")
  }
  cat("  converged: ", x$converged,
      if (!isTRUE(x$converged) && !is.na(x$message)) paste0(" (", x$message, ")") else "",
      "\n", sep = "")
  if (!is.null(x$dp)) {
    cat("\nDifferential privacy (per server, add/remove one subject):\n")
    print(x$dp, row.names = FALSE, digits = digits)
    cat("  the estimates are noised; divergence from an exact fit is expected\n")
  }
  cat("\nCoefficients (natural scale):\n")
  print(round(x$theta.natural, digits))
  invisible(x)
}

#' @rdname ds_nolimits_fit-methods
#' @export
coef.ds_nolimits_fit <- function(object, scale = c("natural", "transformed"), ...) {
  switch(match.arg(scale),
         natural = object$theta.natural,
         transformed = object$theta.transformed)
}

#' @rdname ds_nolimits_fit-methods
#' @export
vcov.ds_nolimits_fit <- function(object, scale = c("natural", "transformed"), ...) {
  switch(match.arg(scale),
         natural = object$vcov.natural,
         transformed = object$vcov.transformed)
}

#' @rdname ds_nolimits_fit-methods
#' @export
confint.ds_nolimits_fit <- function(object, parm = NULL, level = 0.95,
                                    scale = c("natural", "transformed"), ...) {
  scale <- match.arg(scale)
  est <- coef(object, scale = scale)
  V <- vcov(object, scale = scale)
  if (is.null(V)) {
    cli::cli_abort("no covariance matrix is available; refit with {.code se = TRUE}")
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    cli::cli_abort("{.arg level} must be a single number strictly between 0 and 1")
  }
  if (is.null(parm)) parm <- names(est)
  bad <- setdiff(parm, names(est))
  if (length(bad)) cli::cli_abort("unknown parameter{?s}: {.val {bad}}")

  z <- stats::qnorm(1 - (1 - level) / 2)
  sd <- .safe_sqrt(diag(V)[parm])
  out <- cbind(est[parm] - z * sd, est[parm] + z * sd)
  colnames(out) <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                                 trim = TRUE), " %")
  rownames(out) <- parm
  out
}

#' @rdname ds_nolimits_fit-methods
#' @export
summary.ds_nolimits_fit <- function(object, scale = c("natural", "transformed"), ...) {
  scale <- match.arg(scale)
  est <- coef(object, scale = scale)
  V <- vcov(object, scale = scale)
  sd <- if (is.null(V)) rep(NA_real_, length(est)) else .safe_sqrt(diag(V))
  z <- est / sd
  tab <- data.frame(
    Estimate = est,
    `Std. Error` = sd,
    `z value` = z,
    `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)),
    row.names = names(est),
    check.names = FALSE
  )
  out <- list(fit = object, scale = scale, coefficients = tab)
  class(out) <- "summary.ds_nolimits_fit"
  out
}

#' @rdname ds_nolimits_fit-methods
#' @export
print.summary.ds_nolimits_fit <- function(x, digits = 4L, ...) {
  print(x$fit, digits = digits)
  cat("\nCoefficients (", x$scale, " scale):\n", sep = "")
  stats::printCoefmat(as.matrix(x$coefficients), digits = digits, has.Pvalue = TRUE,
                      P.values = TRUE, na.print = "NA")
  cat("\nPer-site contributions:\n")
  print(x$fit$sites, row.names = FALSE)
  invisible(x)
}
