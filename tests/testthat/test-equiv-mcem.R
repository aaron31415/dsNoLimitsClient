# Tier 5 for MCEM, the nested federated EM estimator.
#
# Layer A (the real federated-learning proof, machine precision): the M-step
# Q-function at FIXED draws is a per-subject sum, so summing the per-subject Q
# contributions reproduces the population Q for both parts, value and gradient.
# This validates the primitive summation without the federated loop.
#
# Layer B (Monte-Carlo): a 3-site DSLite federated MCEM fit reproduces the pooled
# fit_model(dm, MCEM()) at the same seed and outer budget. Loose by design: the
# per-site RNG partition differs from the pooled one, so the parameters agree only
# up to Monte-Carlo noise (the exactness is Layer A).

test_that("MCEM Q is per-subject additive to machine precision (Layer A)", {
  nlc_skip_no_julia()
  d <- do.call(rbind, unname(nlc_sites()))
  dm <- nlc_dm(d, NLDS_TEST_MODEL)
  dsNoLimits:::.nlds_helpers()
  sf <- function(n) dsNoLimits:::.nlds_fun(n)

  theta0 <- as.numeric(sf("nlds_theta0")(dm))
  q1 <- as.character(sf("nlds_mcem_q1")(dm))
  q2 <- as.character(sf("nlds_mcem_q2")(dm))
  expect_setequal(c(q1, q2), as.character(sf("nlds_names")(dm)))

  # One E-step (state = nothing) at theta0 gives the fixed draws; each batch is
  # one subject, so summing per-idx equals summing whole sites.
  em <- sf("nlds_mcem_estep")(dm, theta0, 100, as.integer(20260824L), NULL)
  nb <- as.integer(sf("nlds_mcem_nbatches")(em))
  expect_equal(nb, length(unique(d$ID)))

  worst <- 0
  for (pr in list(list(p = "q1", fn = q1), list(p = "q2", fn = q2))) {
    if (!length(pr$fn)) next
    pooled <- as.numeric(sf("nlds_mcem_q")(dm, theta0, em, pr$p, pr$fn))
    fed <- rep(0, length(pooled))
    for (i in seq_len(nb)) {
      fed <- fed + as.numeric(sf("nlds_mcem_q_idx")(dm, theta0, em,
                                                    as.integer(i), pr$p, pr$fn))
    }
    v.rel <- abs(fed[1] - pooled[1]) / max(abs(pooled[1]), 1e-300)
    g.rel <- sqrt(sum((fed[-1] - pooled[-1])^2)) /
      max(sqrt(sum(pooled[-1]^2)), 1e-300)
    nlc_report(sprintf("mcem additivity %s: value_rel=%.3e grad_rel=%.3e (len %d)",
                       pr$p, v.rel, g.rel, length(pr$fn)))
    worst <- max(worst, v.rel, g.rel)
    expect_lt(v.rel, 1e-10)
    expect_lt(g.rel, 1e-10)
  }
  nlc_report(sprintf("mcem Layer A worst additivity rel = %.3e", worst))
})

test_that("federated MCEM reproduces the pooled fit (Layer B)", {
  nlc_skip_no_julia()

  sites <- nlc_sites()
  pooled <- do.call(rbind, unname(sites))
  conns <- ds_lite_conns(sites, aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options())
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  fed <- ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                        estimator = "mcem", checks = FALSE, datasources = conns,
                        em.outer.iters = 15L, em.sample.schedule = 100L,
                        em.seed = 20260824L)

  # The fit object reflects a stochastic, fixed-budget estimator.
  expect_identical(fed$estimator, "mcem")
  expect_identical(fed$outer.iters, 15L)
  expect_true(is.na(fed$loglik))
  expect_true(is.na(fed$converged))
  expect_null(fed$se)
  expect_true(all(is.na(fed$sites$loglik.contribution)))
  expect_equal(names(fed$theta.natural), c("cl", "v", "omega", "sigma"))

  ref <- nlc_mcem_pooled(pooled, 100L, 15L, 20260824L, key = "csv3")
  fn <- names(ref$natural)
  rel <- abs(fed$theta.natural[fn] - ref$natural[fn]) / pmax(abs(ref$natural[fn]), 1e-8)
  nlc_report(sprintf("mcem Layer B: worst param rel=%.3e (%s), rounds=%d",
                     max(rel), fn[which.max(rel)], fed$rounds))
  # Monte-Carlo tolerance, matching the reference federated application.
  expect_lt(max(rel), 5e-2)
})
