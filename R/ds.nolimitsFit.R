## Minimum dsNoLimits version the wire contract below is pinned to.
.DSNOLIMITS_MIN_VERSION <- "0.0.5.9000"

#' @title Fit a nonlinear mixed-effects model across DataSHIELD servers
#' @description Fits a NoLimits nonlinear mixed-effects model to longitudinal
#'   data held on several servers, without any individual-level data leaving
#'   them. Each round the client sends one parameter vector and each site
#'   returns its marginal log-likelihood and gradient; the client sums them and
#'   drives an L-BFGS-B optimisation. For the Laplace, FOCEI, Gauss-Hermite and
#'   naive-pooled estimators the marginal log-likelihood is a sum of independent
#'   per-subject terms, so the federated result is numerically the pooled one.
#'   The same holds for maximum likelihood and maximum a posteriori on a model
#'   without random effects.
#' @details Round budget: one assign, one prep-info aggregate, one aggregate per
#'   optimiser evaluation (typically 10 to 30), plus `2 * p` aggregates when
#'   `se = TRUE`.
#'
#'   Standard errors are Wald: the transformed-scale Hessian is a central
#'   difference of the federated gradient at the optimum, with step
#'   `eps^(1/3) * max(|theta_j|, 1)`, and the natural-scale covariance follows by
#'   the delta method from the same probes. If the observed information is not
#'   positive definite a warning is issued and the affected standard errors are
#'   `NA`; the fit is never failed for it.
#'
#'   Per-site log-likelihood contributions are returned in `$sites`. This is a
#'   deliberate release of the same class as a study-level deviance in an SLMA
#'   meta-analysis. With `estimator = "map"` they are the likelihood-only parts,
#'   so they sum to `$loglik` minus the log-prior rather than to `$loglik`: the
#'   log-prior is data-free, identical at every site, and is added exactly once
#'   (a plain sum would count it once per site). The client checks that the
#'   sites' prior blocks are bit-identical on every round.
#'
#'   With `dp = TRUE` the fit runs against the servers' differentially private
#'   gradient method instead. Every parameter of the privacy mechanism - the
#'   clip, the noise multiplier, the epsilon budget, delta and the cap on the
#'   number of rounds - is a SERVER option set by the data owner; no argument of
#'   this function reaches any of them. The optimiser becomes a fixed schedule of
#'   `dp.rounds` Adam steps: no line search, no convergence test, because the
#'   servers release no objective value and every quantity that could gate one is
#'   noisy. `se` is forced to `FALSE`, `converged` is `NA`, `loglik` is `NA`, and
#'   the fitted object gains a `dp` element reporting each server's epsilon,
#'   delta, release count and remaining budget. Divergence from the exact fit is
#'   expected; that is what the noise is for.
#'
#'   Server functions called: `nolimitsPrepDS`, `nolimitsPrepStringDS`,
#'   `nolimitsPrepInfoDS`, `nolimitsObjGradDS`, and with `dp = TRUE`
#'   `nolimitsObjGradDpDS` and `nolimitsNaturalDS`.
#' @param model a character string. With `model.mode = "registered"` (the
#'   default) the name of a model pre-registered by every data owner, as listed
#'   by [ds.nolimitsModels()]. With `model.mode = "string"` the Julia model
#'   source itself, which is transmitted base64url-encoded and executed on each
#'   node; sites must have enabled the permissive privacy control level.
#' @param data a character string, the name of a server-side data frame.
#' @param id.col a character string, the subject identifier column of `data`.
#' @param time.col a character string, the time column of `data`.
#' @param evid.col `NULL` (the default) for data without dose events, or a
#'   character string naming the event-identifier column of a NONMEM-style
#'   event table (`EVID` 0 = observation, 1 = dose, 2 = reset; dose columns
#'   `AMT`, `RATE` and `CMT`). Needed for ODE models with doses or infusions
#'   given as events.
#' @param estimator one of `"laplace"`, `"focei"`, `"ghq"`, `"pooled"`,
#'   `"mle"`, `"map"`, `"mcem"` or `"saem"`. `"mle"` and `"map"` are for models
#'   with NO random effects and fail on any other model. `"mcem"` and `"saem"` are
#'   sampling-based EM for a random-effects model: a nested loop with a stateful
#'   local E-step (the per-subject posterior draws never leave a site) and a
#'   federated M-step. `"saem"` additionally runs a stateful closed-form M-step on
#'   a deterministic coordinator site (the lowest site name), summing
#'   per-subject-additive sufficient statistics across sites. Both force
#'   `se = FALSE`, report `loglik = NA` and `converged = NA`, and honour the
#'   `em.*` arguments. With `dp = TRUE` each M-step becomes DP-Adam over the
#'   servers' per-subject-clipped, noised gradients (and, for `"saem"`, the stats
#'   round releases noised sufficient statistics while the coordinator's
#'   closed-form update is post-processing): the most privacy-expensive path,
#'   spending `outer.iters * 2 * dp.rounds` releases for `"mcem"` and
#'   `outer.iters * (1 + numerical_parts * dp.rounds)` for `"saem"`.
#' @param ghq.level the Gauss-Hermite quadrature level, 1 to 9; used only by
#'   `estimator = "ghq"`.
#' @param model.mode `"registered"` or `"string"`, see `model`.
#' @param se `TRUE` (the default) to compute Wald standard errors. Forced to
#'   `FALSE` when `dp = TRUE`.
#' @param dp `FALSE` (the default) for the exact fit; `TRUE` to fit against the
#'   servers' differentially private gradient releases. Requires the data owner
#'   to have permitted `nolimitsObjGradDpDS` and configured a budget ledger.
#'   `estimator` must then be one of `"laplace"`, `"focei"`, `"ghq"` or
#'   `"mle"`.
#' @param dp.rounds the number of Adam steps, and therefore of privacy releases,
#'   when `dp = TRUE`. Each server refuses a round number above its own
#'   `dsNoLimits.dp.maxT` option.
#' @param dp.lr the Adam learning rate when `dp = TRUE`. One unit of the
#'   preconditioned scale is roughly one unit of relative change in a parameter,
#'   so the default moves each parameter by about 5 per cent per round.
#' @param em.outer.iters MCEM/SAEM only: the number of outer EM iterations, a
#'   fixed budget with no convergence test. Defaults to 15 for MCEM and 20 for
#'   SAEM when not set.
#' @param em.sample.schedule MCEM/SAEM only: the number of posterior draws per
#'   subject per E-step.
#' @param em.mstep.maxfun MCEM/SAEM only: the L-BFGS-B iteration cap for each
#'   inner numerical M-step part (an approximate M-step is enough).
#' @param em.seed MCEM/SAEM only: the base E-step RNG seed. Each site adds its stable
#'   integer id (the rank of its sorted name, 0-based), so the draws are a pure
#'   function of the seed, the site and the subject, and the federated fit
#'   reproduces the pooled one up to Monte-Carlo noise.
#' @param max.rounds the L-BFGS-B iteration cap.
#' @param newobj a character string, the server-side symbol the prepared model
#'   is bound to. Re-using it lets a second fit with another estimator skip the
#'   preparation round.
#' @param checks `TRUE` (the default) to verify the server package version on
#'   every server before starting. Costs two requests per server.
#' @param datasources a list of [DSI::DSConnection-class] objects obtained after
#'   login. If not specified the default set of connections is used: see
#'   [DSI::datashield.connections_default()].
#' @return an object of class `"ds_nolimits_fit"`: a list with elements
#'   `theta.natural`, `theta.transformed`, `loglik`, `converged`, `message`,
#'   `rounds`, `estimator`, `ghq.level`, `model.hash`, `se`, `vcov.natural`,
#'   `vcov.transformed`, `sites`, `versions`, `timings` and, when `dp = TRUE`,
#'   `dp`. `model.hash` is the canonical model source hash every site agreed on.
#'   See
#'   [print.ds_nolimits_fit()], [coef.ds_nolimits_fit()],
#'   [vcov.ds_nolimits_fit()] and [confint.ds_nolimits_fit()].
#' @author Manuel Huth
#' @examples
#' \dontrun{
#' fit <- ds.nolimitsFit("warfarin_pk", data = "D", id.col = "ID", time.col = "t",
#'                       datasources = conns)
#' summary(fit)
#' }
#' @export
ds.nolimitsFit <- function(model = NULL, data = NULL,
                           id.col = "ID", time.col = "t", evid.col = NULL,
                           estimator = c("laplace", "focei", "ghq", "pooled",
                                         "mle", "map", "mcem", "saem"),
                           ghq.level = 1L,
                           model.mode = c("registered", "string"),
                           se = TRUE, max.rounds = 100L,
                           dp = FALSE, dp.rounds = 50L, dp.lr = 0.05,
                           em.outer.iters = 15L, em.sample.schedule = 100L,
                           em.mstep.maxfun = 12L, em.seed = 20260824L,
                           newobj = "nolimits.prep",
                           checks = TRUE, datasources = NULL) {

  ## 1a. client-side validation, all of it before any DSI call
  estimator <- match.arg(estimator)
  model.mode <- match.arg(model.mode)
  dp <- .check_flag(dp, "dp")
  is.mcem <- identical(estimator, "mcem")
  is.saem <- identical(estimator, "saem")
  is.em <- is.mcem || is.saem
  # SAEM's closed-form M-step converges faster, so it defaults to a longer outer
  # budget than MCEM. Branch only when the analyst did not set it.
  if (is.em && missing(em.outer.iters)) em.outer.iters <- if (is.saem) 20L else 15L
  if (is.em) {
    em.outer.iters <- .check_count(em.outer.iters, "em.outer.iters", 1L, 10000L)
    em.sample.schedule <- .check_count(em.sample.schedule, "em.sample.schedule", 1L, 100000L)
    em.mstep.maxfun <- .check_count(em.mstep.maxfun, "em.mstep.maxfun", 1L, 100000L)
    em.seed <- .check_count(em.seed, "em.seed", 0L, .Machine$integer.max)
  }
  if (dp) {
    if (estimator %in% c("pooled", "map")) {
      cli::cli_abort(c(
        "x" = "{.code estimator = \"{estimator}\"} has no per-subject term and no bounded sensitivity.",
        "i" = "differentially private fits use {.val laplace}, {.val focei}, {.val ghq} or {.val mle}"
      ))
    }
    dp.rounds <- .check_count(dp.rounds, "dp.rounds", 1L, 100000L)
    if (!is.numeric(dp.lr) || length(dp.lr) != 1L || !is.finite(dp.lr) ||
        dp.lr <= 0) {
      cli::cli_abort("{.arg dp.lr} must be a single positive number")
    }
  }
  if (!is.character(model) || length(model) != 1L || is.na(model) || !nzchar(model)) {
    cli::cli_abort("{.arg model} must be a single non-empty character string")
  }
  if (model.mode == "registered" && !grepl("^[A-Za-z][A-Za-z0-9._-]*$", model)) {
    cli::cli_abort(c(
      "x" = "{.arg model} must be a registered model name, not model source.",
      "i" = "use {.code model.mode = \"string\"} to send Julia source."
    ))
  }
  .check_name(data, "data")
  .check_name(id.col, "id.col")
  .check_name(time.col, "time.col")
  if (!is.null(evid.col)) .check_name(evid.col, "evid.col")
  .check_name(newobj, "newobj")
  ghq.level <- .check_count(ghq.level, "ghq.level", 1L, 9L)
  max.rounds <- .check_count(max.rounds, "max.rounds", 1L, 100000L)
  se <- .check_flag(se, "se")
  checks <- .check_flag(checks, "checks")
  if (dp && se) {
    cli::cli_inform(c(
      "i" = "{.code se = TRUE} is ignored under {.code dp = TRUE}.",
      "*" = "Wald standard errors need {.val 2p} extra gradient releases and a Hessian of a noisy gradient; neither is meaningful here."
    ))
    se <- FALSE
  }
  if (is.em && se) {
    cli::cli_inform(c(
      "i" = "{.code se = TRUE} is ignored for {.code estimator = \"{estimator}\"}.",
      "*" = "The sampling-based EM estimators are stochastic: they track a Monte-Carlo Q at varying draws, not a comparable marginal likelihood, so Wald standard errors are not defined here."
    ))
    se <- FALSE
  }

  ## 1b. connections
  datasources <- .setDatasources(datasources)
  sites <- names(datasources)
  ## Stable per-site id: rank of the sorted site name, 0-based. Makes the E-step
  ## draws a pure function of (em.seed, site.id, subject), so federated == pooled
  ## up to Monte-Carlo noise.
  site.ids <- stats::setNames(match(sites, sort(sites)) - 1L, sites)
  ## SAEM's stateful closed-form M-step runs on ONE deterministic coordinator: the
  ## lowest site name (site id 0), the same rule as the MAP prior carrier.
  coordinator <- sort(sites)[1L]

  ## The federated loop is thousands of short round trips; the default polling
  ## backoff dominates the wall clock if left alone.
  old.opts <- options(datashield.polling.sleep.0 = 0.02,
                      datashield.polling.sleep.1 = 0.1,
                      datashield.progress = FALSE)
  on.exit(options(old.opts), add = TRUE)

  t.start <- proc.time()[["elapsed"]]

  ## 2. server package version guard
  if (checks) .require_server_pkg(datasources, "dsNoLimits", .DSNOLIMITS_MIN_VERSION)

  ## 3. prepare the model on every site
  cally <- if (model.mode == "registered") {
    call("nolimitsPrepDS", data, model, id.col, time.col)
  } else {
    call("nolimitsPrepStringDS", data, .ds_encode(model), id.col, time.col)
  }
  # Only added when used, so calls without events are unchanged.
  if (!is.null(evid.col)) cally$evid.col <- evid.col
  # cally is always a call object, so use the explicit expression form: plain
  # datashield.assign() dispatches on the value's type and would misroute a call
  # bound to a variable. No behaviour change today; removes a latent misdispatch.
  DSI::datashield.assign.expr(datasources, newobj, cally)

  ## 4. read back everything except the Julia handle
  infos <- .agg(datasources, call("nolimitsPrepInfoDS", newobj), "model preparation")
  .check_complete(infos, sites, "model preparation")

  ## 5. every site must be running the same model: model-derived names/theta0/
  ## scale agree, AND every site resolved the byte-identical model source.
  ref <- agree(infos)
  model.hash <- .agree_model(infos)
  theta0 <- as.numeric(ref$theta0)
  s <- as.numeric(ref$scale)
  par.names <- as.character(ref$names)
  p <- length(theta0)
  if (p < 1L || length(s) != p || length(par.names) != p) {
    cli::cli_abort("the servers returned an inconsistent parameter description")
  }
  t.prep <- proc.time()[["elapsed"]]

  ## 6. optimise on the preconditioned scale z, theta = theta0 + s * z
  dp.report <- NULL
  if (dp && !is.em) {
    ## Fixed-schedule Adam on the private gradients. The servers clip and noise
    ## in the preconditioned coordinates, so their gradient is already dL/dz and
    ## must NOT be multiplied by s again.
    spec <- list(newobj = newobj, estimator = estimator, ghq.level = ghq.level,
                 p = p, sites = sites,
                 aggregate_fun = function(cl) {
                   .agg(datasources, cl, "differentially private gradient")
                 })
    run <- .dp_optimize(p, theta0, s, dp.lr, dp.rounds, spec)
    theta.t <- theta0 + s * run$z
    nat <- .agg(datasources, call("nolimitsNaturalDS", newobj,
                                  .ds_num_encode(theta.t)),
                "natural-scale parameters")
    .check_complete(nat, sites, "natural-scale parameters")
    theta.nat <- as.numeric(nat[[1L]]$natural)
    if (length(theta.nat) != p) {
      cli::cli_abort("the servers returned an unusable natural-scale parameter vector")
    }
    dp.report <- .dp_report(run$last$per.site)
    n.rounds <- run$rounds + 1L
    loglik <- NA_real_
    per.site <- NULL
    converged <- NA
    fit.message <- paste0("fixed-schedule Adam, ", dp.rounds,
                          " differentially private releases; no convergence ",
                          "test is possible on noisy gradients")
  } else if (is.mcem) {
    ## Nested federated MCEM: a local, stateful E-step per site each outer
    ## iteration, then federated L-BFGS-B over the q1 then q2 coordinates at the
    ## fixed draws. The q1/q2 partition is model-derived, so it must agree across
    ## sites; the E-step is the one place a per-site argument (the site id) forces
    ## a short loop over connections.
    em <- .agree_em(infos, par.names)
    em.spec <- list(newobj = newobj, outer.iters = em.outer.iters,
                    sample.schedule = em.sample.schedule, maxfun = em.mstep.maxfun,
                    seed = em.seed, dp.rounds = dp.rounds, dp.lr = dp.lr,
                    agg_fun = function(cl) .agg(datasources, cl, "Q-function"),
                    dp_agg = function(cl) .agg(datasources, cl,
                                               "differentially private gradient"))
    run <- if (dp) {
      .em_run_dp(p, theta0, s, em$q1i, em$q2i, datasources, sites, site.ids, em.spec)
    } else {
      .em_run(p, theta0, s, em$q1i, em$q2i, datasources, sites, site.ids, em.spec)
    }
    theta.t <- run$theta.t
    nat <- .agg(datasources, call("nolimitsNaturalDS", newobj,
                                  .ds_num_encode(theta.t)),
                "natural-scale parameters")
    .check_complete(nat, sites, "natural-scale parameters")
    theta.nat <- as.numeric(nat[[1L]]$natural)
    if (length(theta.nat) != p) {
      cli::cli_abort("the servers returned an unusable natural-scale parameter vector")
    }
    n.rounds <- run$rounds
    loglik <- NA_real_
    per.site <- NULL
    converged <- NA
    if (dp) {
      dp.report <- .dp_report(run$last$per.site)
      fit.message <- paste0("MCEM, ", em.outer.iters, " outer iterations, ",
                            run$releases, " differentially private releases; ",
                            "no convergence test (fixed budget, noisy gradients)")
    } else {
      fit.message <- paste0("MCEM, ", em.outer.iters,
                            " outer iterations; no convergence test (fixed budget)")
    }
  } else if (is.saem) {
    ## Nested federated SAEM (hybrid M-step). Each outer iteration: a local
    ## E-step per site, a stats round summed coordinate-wise, a stateful
    ## closed-form M-step on the coordinator (lowest site name), then federated
    ## L-BFGS-B over the numerical params of each part. The closed-form/numerical
    ## split, like the q1/q2 partition, is model-derived and agreed across sites.
    em <- .agree_em(infos, par.names, estimator = "saem")
    em.spec <- list(newobj = newobj, outer.iters = em.outer.iters,
                    sample.schedule = em.sample.schedule, maxfun = em.mstep.maxfun,
                    seed = em.seed, dp.rounds = dp.rounds, dp.lr = dp.lr,
                    agg_fun = function(cl) .agg(datasources, cl, "Q-function"),
                    dp_agg = function(cl) .agg(datasources, cl,
                                               "differentially private gradient"))
    run <- if (dp) {
      .saem_run_dp(p, theta0, s, par.names, em, datasources, sites, site.ids,
                   coordinator, em.spec)
    } else {
      .saem_run(p, theta0, s, par.names, em, datasources, sites, site.ids,
                coordinator, em.spec)
    }
    theta.t <- run$theta.t
    nat <- .agg(datasources, call("nolimitsNaturalDS", newobj,
                                  .ds_num_encode(theta.t)),
                "natural-scale parameters")
    .check_complete(nat, sites, "natural-scale parameters")
    theta.nat <- as.numeric(nat[[1L]]$natural)
    if (length(theta.nat) != p) {
      cli::cli_abort("the servers returned an unusable natural-scale parameter vector")
    }
    n.rounds <- run$rounds
    loglik <- NA_real_
    per.site <- NULL
    converged <- NA
    if (dp) {
      dp.report <- .dp_report(run$last$per.site)
      fit.message <- paste0("SAEM, ", em.outer.iters, " outer iterations ",
                            "(coordinator ", coordinator, "), ", run$releases,
                            " differentially private releases; no convergence ",
                            "test (fixed budget, noisy gradients)")
    } else {
      fit.message <- paste0("SAEM, ", em.outer.iters,
                            " outer iterations (coordinator ", coordinator,
                            "); no convergence test (fixed budget)")
    }
  } else {
    fed <- .fed_closure(theta0, s, newobj, estimator, ghq.level, sites,
                        function(cl) .agg(datasources, cl, "objective and gradient"))
    # factr / pgtol pinned to the SciPy defaults the reference federated runs used.
    opt <- stats::optim(
      rep(0, p),
      fn = function(z) fed$eval(z)$fn,
      gr = function(z) fed$eval(z)$gr,
      method = "L-BFGS-B",
      control = list(maxit = max.rounds, factr = 1e7, pgtol = 1e-5)
    )

    ## Final evaluation at the optimum; the memo makes this free whenever optim
    ## ended there, which it usually does. Done before the SE probes so the memo
    ## is still warm.
    final <- fed$eval(opt$par)
    if (!isTRUE(final$finite)) {
      cli::cli_abort("the objective is not finite at the reported optimum")
    }
    theta.t <- theta0 + s * opt$par
    theta.nat <- final$natural
    loglik <- final$value
    per.site <- final$per.site
    converged <- identical(opt$convergence, 0L)
    fit.message <- if (is.null(opt$message)) NA_character_ else opt$message
  }
  names(theta.t) <- par.names
  names(theta.nat) <- par.names
  t.optim <- proc.time()[["elapsed"]]

  ## 7. Wald standard errors from central differences of the federated gradient
  vcov.t <- NULL
  vcov.nat <- NULL
  se.nat <- NULL
  if (se) {
    h <- .Machine$double.eps^(1 / 3) * pmax(abs(theta.t), 1)
    G <- matrix(NA_real_, p, p)
    J <- matrix(NA_real_, p, p)
    ok <- TRUE
    for (j in seq_len(p)) {
      tp <- theta.t; tp[j] <- tp[j] + h[j]
      tm <- theta.t; tm[j] <- tm[j] - h[j]
      up <- fed$probe(tp)
      dn <- fed$probe(tm)
      if (!up$finite || !dn$finite) {
        ok <- FALSE
        next
      }
      G[, j] <- (up$grad - dn$grad) / (2 * h[j])
      J[, j] <- (up$natural - dn$natural) / (2 * h[j])
    }
    if (!ok) {
      warning("the objective is not finite at some standard-error probe; ",
              "standard errors are not available.", call. = FALSE)
    } else {
      H <- (G + t(G)) / 2
      if (inherits(try(chol(-H), silent = TRUE), "try-error")) {
        warning("the observed information is not positive definite; ",
                "standard errors may be unreliable or NA.", call. = FALSE)
      }
      vcov.t <- try(solve(-H), silent = TRUE)
      if (inherits(vcov.t, "try-error")) {
        warning("the observed information is singular; standard errors are NA.",
                call. = FALSE)
        vcov.t <- NULL
      } else {
        dimnames(vcov.t) <- list(par.names, par.names)
        vcov.nat <- J %*% vcov.t %*% t(J)
        dimnames(vcov.nat) <- list(par.names, par.names)
        se.nat <- stats::setNames(.safe_sqrt(diag(vcov.nat)), par.names)
      }
    }
    if (is.null(se.nat)) se.nat <- stats::setNames(rep(NA_real_, p), par.names)
  }
  t.se <- proc.time()[["elapsed"]]

  ## 8. assemble
  out <- list(
    theta.natural = theta.nat,
    theta.transformed = theta.t,
    loglik = loglik,
    converged = converged,
    message = fit.message,
    rounds = if (dp || is.em) n.rounds else fed$rounds(),
    estimator = estimator,
    ghq.level = ghq.level,
    # EM-only (MCEM/SAEM): the fixed outer-iteration count and the EM config.
    outer.iters = if (is.em) em.outer.iters else NULL,
    # SAEM-only: the deterministic coordinator site (transparency).
    coordinator = if (is.saem) coordinator else NULL,
    em = if (is.em) list(outer.iters = em.outer.iters,
                         sample.schedule = em.sample.schedule,
                         mstep.maxfun = em.mstep.maxfun, seed = em.seed) else NULL,
    # The agreed canonical model hash: the analyst's receipt that every site ran
    # the byte-identical model.
    model.hash = model.hash,
    se = se.nat,
    vcov.natural = vcov.nat,
    vcov.transformed = vcov.t,
    sites = data.frame(
      site = sites,
      n.subjects = vapply(infos[sites], function(x) as.integer(.field1(x, "n.subjects", NA_integer_)), integer(1)),
      n.obs = vapply(infos[sites], function(x) as.integer(.field1(x, "n.obs", NA_integer_)), integer(1)),
      # No per-site log-likelihood under DP or EM: DP releases no objective value,
      # and MCEM/SAEM track a Monte-Carlo Q at varying draws, not a loglik.
      loglik.contribution = if (dp || is.em) rep(NA_real_, length(sites)) else
        vapply(per.site[sites], function(x) as.numeric(.field1(x, "value", NA_real_)), numeric(1)),
      row.names = NULL,
      stringsAsFactors = FALSE
    ),
    dp = dp.report,
    versions = lapply(infos, `[[`, "versions"),
    timings = c(prep = t.prep - t.start, optim = t.optim - t.prep,
                se = t.se - t.optim, total = t.se - t.start)
  )
  class(out) <- "ds_nolimits_fit"
  out
}
