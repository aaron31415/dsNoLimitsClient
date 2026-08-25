# The model fixture and the pooled oracle for tiers 4 and 5.
#
# ONE model is reused everywhere: every distinct @Model pays Julia code
# generation again. Closed-form one-compartment IV with a log-normal random
# effect on clearance - no ODE solver runs at all. This string is the copy of
# the constant in the dsNoLimits test suite; the two must stay identical,
# because a federated fit and its pooled oracle have to be the same model.
NLDS_TEST_MODEL <- '
@fixedEffects begin
    cl    = RealNumber(1.0)
    v     = RealNumber(10.0)
    omega = RealNumber(0.3, scale=:log)
    sigma = RealNumber(0.5, scale=:log)
end

@covariates begin
    t    = Covariate()
    Dose = ConstantCovariate(constant_on=:ID)
end

@randomEffects begin
    eta_cl = RandomEffect(LogNormal(0.0, omega); column=:ID)
end

@formulas begin
    cli = cl * eta_cl
    conc ~ Normal(Dose / v * exp(-(cli / v) * t), sigma)
end
'

# The second fixture: NO random effects, a weakly informative prior on every
# parameter. `mle` and `map` require a model without random effects. Fitting it
# to the RE-simulated fixture data is a deliberately misspecified fit; the oracle
# is the SAME fit pooled, so the equivalence statement is unaffected. p = 3, so
# the 6-subject sites need nfilter.glm >= 0.5 (the harness uses 0.8).
NLDS_TEST_MODEL_NORE <- '
@fixedEffects begin
    cl    = RealNumber(1.0,  prior=LogNormal(log(1.0), 0.5))
    v     = RealNumber(10.0, prior=LogNormal(log(10.0), 0.5))
    sigma = RealNumber(0.5, scale=:log, prior=LogNormal(log(0.5), 0.5))
end

@covariates begin
    t    = Covariate()
    Dose = ConstantCovariate(constant_on=:ID)
end

@formulas begin
    conc ~ Normal(Dose / v * exp(-(cl / v) * t), sigma)
end
'

# 6 subjects x 6 observations each, synthetic, vendored from the reference
# federated application's example data.
nlc_site <- function(which) {
  utils::read.csv(test_path("fixtures", paste0("site_", which, ".csv")),
                  stringsAsFactors = FALSE)
}

nlc_sites <- function(which = c("a", "b", "c")) {
  stats::setNames(lapply(which, nlc_site), which)
}

# Seeded generator, ~30 subjects. Large enough that the CANONICAL
# nfilter.glm = 0.33 saturation gate passes on every site of an unequal split,
# which the 6-subject CSV fixtures are not.
nlc_sim_data <- function(n = 30L, seed = 20260823L) {
  withr::with_seed(seed, {
    times <- c(0.5, 1, 2, 4, 8, 12)
    eta <- exp(stats::rnorm(n, 0, 0.3))
    do.call(rbind, lapply(seq_len(n), function(i) {
      cli <- 1 * eta[i]
      data.frame(ID = sprintf("S%03d", i), t = times, Dose = 100,
                 conc = 100 / 10 * exp(-(cli / 10) * times) +
                   stats::rnorm(length(times), 0, 0.5),
                 stringsAsFactors = FALSE)
    }))
  })
}

# Unequal split by subject: 5 / 10 / 15.
nlc_sim_split <- function(d = nlc_sim_data(), sizes = c(5L, 10L, 15L)) {
  ids <- unique(d$ID)
  stopifnot(sum(sizes) == length(ids))
  g <- rep(seq_along(sizes), times = sizes)
  stats::setNames(lapply(seq_along(sizes),
                         function(k) d[d$ID %in% ids[g == k], , drop = FALSE]),
                  paste0("u", seq_along(sizes)))
}

# TRUE only if NoLimits actually boots: Julia present AND NoLimits.jl installed
# and loadable. Cached so booting is attempted once. On CI, where no Julia /
# NoLimits environment is provisioned, this is FALSE and the Julia tiers skip
# instead of erroring.
.nlc_ready <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- isTRUE(tryCatch(
        requireNamespace("NoLimitsR", quietly = TRUE) &&
          requireNamespace("JuliaConnectoR", quietly = TRUE) &&
          isTRUE(JuliaConnectoR::juliaSetupOk()) &&
          !is.null(NoLimitsR::nolimits()),
        error = function(e) FALSE))
    }
    cached
  }
})

nlc_skip_no_julia <- function() {
  skip_on_cran()
  skip_if_not_installed("dsNoLimits")
  skip_if_not(.nlc_ready(), "Julia / NoLimits not available")
}

## ---------------------------------------------------------------- oracle ----
## The pooled reference: NoLimitsR's own fit_model on the rbind-ed data, which
## is an independent implementation of the estimator. Only the objective probe
## reuses the server package's Julia helper - the parameter transform it wraps
## is shared by construction, so a second copy would add no coverage, and the
## quantity under test is the federated SUM against a single-process
## evaluation, not the transform.

.nlc_cache <- new.env(parent = emptyenv())

nlc_dm <- function(d, model = NLDS_TEST_MODEL) {
  nlc_skip_no_julia()
  nl <- NoLimitsR::nolimits()
  nl$DataModel(NoLimitsR::nl_model(model), NoLimitsR::nl_data(d),
               primary_id = NoLimitsR::jl_sym("ID"),
               time_col = NoLimitsR::jl_sym("t"))
}

# [value; grad(p)] of the pooled objective at a transformed-scale theta. For
# "map" that is the likelihood plus the log-prior ONCE, which is what a
# single-process fit_model(dm, MAP()) optimises - the whole point of the
# federated split being tested.
nlc_objgrad <- function(dm, theta, estimator = "laplace", ghq.level = 1L) {
  nlc_skip_no_julia()
  dsNoLimits:::.nlds_helpers()
  flat <- as.numeric(dsNoLimits:::.nlds_fun("nlds_objgrad")(
    dm, as.numeric(theta),
    dsNoLimits:::.nlds_fun("nlds_method")(estimator, as.integer(ghq.level))
  ))
  out <- list(value = flat[1L], grad = flat[seq(2L, (length(flat) + 1L) / 2L)])
  if (identical(estimator, "map")) {
    pr <- as.numeric(dsNoLimits:::.nlds_fun("nlds_prior")(dm, as.numeric(theta)))
    out$value <- out$value + pr[1L]
    out$grad <- out$grad + pr[-1L]
  }
  out
}

# Pooled fit, memoised on `key` because a fit costs tens of seconds.
nlc_pooled <- function(d, estimator = "laplace", ghq.level = 1L, key = NULL,
                       model = NLDS_TEST_MODEL) {
  nlc_skip_no_julia()
  key <- paste(key %||% deparse(substitute(d)), estimator, ghq.level)
  if (!is.null(.nlc_cache[[key]])) return(.nlc_cache[[key]])
  nl <- NoLimitsR::nolimits()
  dm <- nlc_dm(d, model)
  est <- switch(estimator,
                laplace = nl$Laplace(), focei = nl$FOCEI(),
                ghq = nl$GHQuadrature(level = as.integer(ghq.level)),
                pooled = nl$Pooled(),
                mle = nl$MLE(), map = nl$MAP())
  fit <- nl$fit_model(dm, est)
  NoLimitsR::nl_eval("nlc_vec(x) = Vector{Float64}(x)")
  vec <- JuliaConnectoR::juliaEval("NoLimitsSession.nlc_vec")
  theta.t <- as.numeric(vec(nl$get_params(fit, scale = NoLimitsR::jl_sym("transformed"))))
  out <- list(dm = dm, theta.t = theta.t,
              value = nlc_objgrad(dm, theta.t, estimator, ghq.level)$value)
  .nlc_cache[[key]] <- out
  out
}

# Pooled MCEM reference: fit_model(dm, MCEM(...)) in one process, seeded so the
# oracle is deterministic. MCEM is stochastic and the per-site RNG partition
# differs from the pooled one, so federated matches this only up to Monte-Carlo
# noise (the machine-precision claim is the Layer A additivity test). Returns the
# natural-scale parameters, named.
nlc_mcem_pooled <- function(d, sample.schedule = 100L, outer.iters = 15L,
                            seed = 20260824L, key = NULL) {
  nlc_skip_no_julia()
  key <- paste("mcem", key %||% deparse(substitute(d)), sample.schedule,
               outer.iters, seed)
  if (!is.null(.nlc_cache[[key]])) return(.nlc_cache[[key]])
  nl <- NoLimitsR::nolimits()
  dm <- nlc_dm(d, NLDS_TEST_MODEL)
  JuliaConnectoR::juliaEval(sprintf("import Random; Random.seed!(%d)",
                                    as.integer(seed)))
  method <- nl$MCEM(sample_schedule = as.integer(sample.schedule),
                    maxiters = as.integer(outer.iters))
  fit <- nl$fit_model(dm, method)
  NoLimitsR::nl_eval("nlc_vec(x) = Vector{Float64}(x)")
  vec <- JuliaConnectoR::juliaEval("NoLimitsSession.nlc_vec")
  theta.t <- as.numeric(vec(nl$get_params(fit, scale = NoLimitsR::jl_sym("transformed"))))
  dsNoLimits:::.nlds_helpers()
  nat <- as.numeric(dsNoLimits:::.nlds_fun("nlds_natural_vec")(dm, theta.t))
  names(nat) <- as.character(dsNoLimits:::.nlds_fun("nlds_names")(dm))
  out <- list(dm = dm, theta.t = theta.t, natural = nat)
  .nlc_cache[[key]] <- out
  out
}

# Pooled SAEM reference: fit_model(dm, SAEM(...)) in one process, seeded so the
# oracle is deterministic. The demo SAEM method (no SA burn-in, plain-maximization
# numerical M-step, convergence window past the budget) matches the federated
# driver's method exactly. Returns the natural-scale parameters, named. Like MCEM,
# SAEM is stochastic and the per-site RNG partition differs from the pooled one, so
# federated matches this only up to Monte-Carlo noise.
nlc_saem_pooled <- function(d, sample.schedule = 100L, outer.iters = 20L,
                            seed = 20260824L, key = NULL) {
  nlc_skip_no_julia()
  key <- paste("saem", key %||% deparse(substitute(d)), sample.schedule,
               outer.iters, seed)
  if (!is.null(.nlc_cache[[key]])) return(.nlc_cache[[key]])
  nl <- NoLimitsR::nolimits()
  dm <- nlc_dm(d, NLDS_TEST_MODEL)
  JuliaConnectoR::juliaEval(sprintf("import Random; Random.seed!(%d)",
                                    as.integer(seed)))
  dsNoLimits:::.nlds_helpers()
  method <- dsNoLimits:::.nlds_fun("nlds_saem_method")(as.integer(outer.iters))
  # Match the federated E-step sample schedule (MCEM sampler drives SAEM's E-step).
  fit <- nl$fit_model(dm, method)
  NoLimitsR::nl_eval("nlc_vec(x) = Vector{Float64}(x)")
  vec <- JuliaConnectoR::juliaEval("NoLimitsSession.nlc_vec")
  theta.t <- as.numeric(vec(nl$get_params(fit, scale = NoLimitsR::jl_sym("transformed"))))
  nat <- as.numeric(dsNoLimits:::.nlds_fun("nlds_natural_vec")(dm, theta.t))
  names(nat) <- as.character(dsNoLimits:::.nlds_fun("nlds_names")(dm))
  out <- list(dm = dm, theta.t = theta.t, natural = nat)
  .nlc_cache[[key]] <- out
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# Relative difference, safe at zero.
nlc_reldiff <- function(a, b) abs(a - b) / max(abs(a), abs(b), 1)

# The observed aggregate-round count of the 3-site Laplace fit with se = TRUE:
# optimiser evaluations plus 2p standard-error probes. Pinned as a regression
# tripwire, not as a specification.
NLC_EXPECTED_ROUNDS <- 19L

# Report a measured acceptance number into the test log, so the gate numbers do
# not have to be re-derived by rerunning the suite.
nlc_report <- function(...) cat("\n[equiv]", ..., "\n")
