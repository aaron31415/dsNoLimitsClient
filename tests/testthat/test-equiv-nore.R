# Tier 5 for the no-random-effects estimators, `mle` and `map`.
#
# `mle` is additive like the other four. `map` is not: its log-prior is one
# data-free term that a plain sum over S sites would count S times. The sites
# therefore return the likelihood-only part plus a shared prior block, and the
# client adds the prior once. The ADDITIVITY test below is what catches a
# double-counted prior: fitting {a} and fitting {a, b} must agree with the
# corresponding pooled fits, and a doubled prior would break the two-site case
# while leaving the single-site case correct.

nlc_nore_fit <- function(conns, est, ...) {
  ds.nolimitsFit(NLDS_TEST_MODEL_NORE, data = "D", model.mode = "string",
                 estimator = est, se = FALSE, checks = FALSE,
                 datasources = conns, ...)
}

test_that("federated == pooled for mle and map on a 3-site split", {
  nlc_skip_no_julia()

  sites <- nlc_sites()
  pooled <- do.call(rbind, unname(sites))
  conns <- ds_lite_conns(sites, aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options())
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  for (est in c("mle", "map")) {
    fed <- nlc_nore_fit(conns, est)
    ref <- nlc_pooled(pooled, est, 1L, key = "csv3nore", model = NLDS_TEST_MODEL_NORE)

    # 1. the transport is the pooled objective: the sum of the three sites'
    #    likelihood parts plus the prior once IS the pooled objective there.
    at.fed <- nlc_objgrad(ref$dm, fed$theta.transformed, est)$value
    expect_equal(fed$loglik, at.fed, tolerance = 1e-9)

    # 2. per-site contributions. For map they are likelihood-only, so they sum
    #    to the objective MINUS the log-prior, and that gap is exactly one
    #    prior evaluation - not two, not three.
    lik.only <- nlc_objgrad(ref$dm, fed$theta.transformed, "mle")$value
    expect_equal(sum(fed$sites$loglik.contribution), lik.only, tolerance = 1e-9)
    expect_equal(sum(fed$sites$n.subjects), 18L)
    expect_equal(sum(fed$sites$n.obs), nrow(pooled))

    # 3. objective agreement at the two optima, and parameters.
    obj.rel <- nlc_reldiff(at.fed, ref$value)
    dtheta <- abs(fed$theta.transformed - ref$theta.t)
    nlc_report(sprintf("%-8s obj.rel=%.3e  max|dtheta|=%.3e  rounds=%d",
                       est, obj.rel, max(dtheta), fed$rounds))
    expect_lt(obj.rel, 1e-6)
    expect_lt(max(dtheta), 1e-3)

    expect_true(fed$converged)
    expect_equal(names(fed$theta.transformed), c("cl", "v", "sigma"))
    expect_lt(fed$rounds, 60L)
  }
})

test_that("map really is a different objective from mle", {
  nlc_skip_no_julia()
  d <- do.call(rbind, unname(nlc_sites()))
  a <- nlc_pooled(d, "mle", 1L, key = "csv3nore", model = NLDS_TEST_MODEL_NORE)
  b <- nlc_pooled(d, "map", 1L, key = "csv3nore", model = NLDS_TEST_MODEL_NORE)
  expect_false(isTRUE(all.equal(a$theta.t, b$theta.t)))
  # The prior is worth a real amount of objective at the MLE optimum.
  pr <- nlc_objgrad(a$dm, a$theta.t, "map")$value -
    nlc_objgrad(a$dm, a$theta.t, "mle")$value
  nlc_report(sprintf("map prior at the mle optimum = %.6f", pr))
  expect_true(is.finite(pr) && abs(pr) > 1e-6)
})

test_that("single-site federation == the plain pooled fit for mle and map", {
  nlc_skip_no_julia()

  one <- nlc_sites("b")
  conns <- ds_lite_conns(one, aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options())
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  for (est in c("mle", "map")) {
    fed <- nlc_nore_fit(conns, est)
    ref <- nlc_pooled(one$b, est, 1L, key = "site-b-nore",
                      model = NLDS_TEST_MODEL_NORE)
    at.fed <- nlc_objgrad(ref$dm, fed$theta.transformed, est)$value
    obj.rel <- nlc_reldiff(at.fed, ref$value)
    nlc_report(sprintf("%-8s single-site obj.rel=%.3e  max|dtheta|=%.3e",
                       est, obj.rel, max(abs(fed$theta.transformed - ref$theta.t))))
    expect_equal(fed$loglik, at.fed, tolerance = 1e-9)
    expect_lt(obj.rel, 1e-6)
    expect_lt(max(abs(fed$theta.transformed - ref$theta.t)), 1e-3)
  }
})

test_that("additivity: {a} and {a,b} both match their pooled fits", {
  nlc_skip_no_julia()

  sites <- nlc_sites(c("a", "b"))
  for (est in c("mle", "map")) {
    for (k in 1:2) {
      part <- sites[seq_len(k)]
      conns <- ds_lite_conns(part, aggregate = DS_AGG, assign = DS_ASN,
                             options = ds_test_options())
      fed <- nlc_nore_fit(conns, est)
      DSI::datashield.logout(conns)

      d <- do.call(rbind, unname(part))
      ref <- nlc_pooled(d, est, 1L, key = paste0("nore-add-", k),
                        model = NLDS_TEST_MODEL_NORE)
      at.fed <- nlc_objgrad(ref$dm, fed$theta.transformed, est)$value
      obj.rel <- nlc_reldiff(at.fed, ref$value)
      nlc_report(sprintf("%-8s additivity k=%d obj.rel=%.3e  max|dtheta|=%.3e",
                         est, k, obj.rel,
                         max(abs(fed$theta.transformed - ref$theta.t))))
      # A prior counted once per site agrees with the pooled fit at k = 1 and
      # disagrees at k = 2. Both assertions have to hold.
      expect_equal(fed$loglik, at.fed, tolerance = 1e-9)
      expect_lt(obj.rel, 1e-6)
      expect_lt(max(abs(fed$theta.transformed - ref$theta.t)), 1e-3)
    }
  }
})

test_that("a random-effects model with mle or map fails legibly", {
  nlc_skip_no_julia()

  conns <- ds_lite_conns(nlc_sites("a"), aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options())
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  for (est in c("mle", "map")) {
    e <- tryCatch(
      suppressWarnings(
        ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                       estimator = est, se = FALSE, checks = FALSE,
                       datasources = conns)),
      error = function(e) e)
    expect_s3_class(e, "error")
    # DSI aborts the aggregate itself, so the reason travels in
    # datashield.errors() rather than in the condition message. The failure is
    # reported generically: the Julia precondition text is not relayed, by
    # design, so that data values can never leak through an error string.
    msg <- paste(unlist(DSI::datashield.errors()), collapse = " ")
    expect_match(msg, "FAILED: the objective could not be evaluated on this server")
  }
})

test_that("map standard errors match a local Hessian of the pooled objective", {
  nlc_skip_no_julia()

  sites <- nlc_sites()
  pooled <- do.call(rbind, unname(sites))
  conns <- ds_lite_conns(sites, aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options())
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  fed <- ds.nolimitsFit(NLDS_TEST_MODEL_NORE, data = "D", model.mode = "string",
                        estimator = "map", se = TRUE, checks = FALSE,
                        datasources = conns)
  expect_false(is.null(fed$vcov.transformed))

  # The same 2p-probe formula, on the pooled objective, in one process.
  dm <- nlc_pooled(pooled, "map", 1L, key = "csv3nore",
                   model = NLDS_TEST_MODEL_NORE)$dm
  th <- as.numeric(fed$theta.transformed)
  p <- length(th)
  h <- .Machine$double.eps^(1 / 3) * pmax(abs(th), 1)
  G <- matrix(NA_real_, p, p)
  for (j in seq_len(p)) {
    tp <- th; tp[j] <- tp[j] + h[j]
    tm <- th; tm[j] <- tm[j] - h[j]
    G[, j] <- (nlc_objgrad(dm, tp, "map")$grad -
                 nlc_objgrad(dm, tm, "map")$grad) / (2 * h[j])
  }
  ref <- solve(-((G + t(G)) / 2))
  rel <- max(abs(fed$vcov.transformed - ref) / pmax(abs(ref), 1e-12))
  nlc_report(sprintf("map se worst rel=%.3e  se.t=%s", rel,
                     paste(sprintf("%.5f", sqrt(diag(fed$vcov.transformed))),
                           collapse = " ")))
  expect_lt(rel, 1e-4)
  # 2p extra rounds on top of the optimiser evaluations.
  nlc_report(sprintf("map rounds with se = %d", fed$rounds))
})
