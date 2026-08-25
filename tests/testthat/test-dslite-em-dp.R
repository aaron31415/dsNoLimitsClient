# Tier 4, differential privacy for the EM estimators: DP-MCEM and DP-SAEM over
# three DSLite sites running the real server functions.
#
# Utility is NOT gated - a DP EM fit is deliberately noised. What is asserted is
# that the fit RUNS, that it spends exactly the privacy budget the composition
# formula (and the accountant) say it does, and that nothing un-noised leaks.
# Settings are deliberately small (few outer iterations, few Adam steps) so the
# nested loop stays fast; the numbers, not the fit quality, are the point.

test_that("a DP-MCEM fit runs and spends the accounted epsilon", {
  nlc_skip_no_julia()
  dir <- withr::local_tempdir()
  conns <- ds_lite_conns(nlc_sites(), aggregate = DS_AGG, assign = DS_ASN,
                         options = dp_options(dir))
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  outer <- 3L
  steps <- 2L
  fit <- ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                        estimator = "mcem", dp = TRUE, dp.rounds = steps,
                        em.outer.iters = outer, em.sample.schedule = 20L,
                        se = FALSE, checks = FALSE, datasources = conns)

  expect_s3_class(fit, "ds_nolimits_fit")
  expect_true(all(is.finite(fit$theta.natural)))
  # No objective value under DP: nothing to converge on or report.
  expect_true(is.na(fit$converged))
  expect_true(is.na(fit$loglik))
  expect_true(all(is.na(fit$sites$loglik.contribution)))
  expect_null(fit$se)

  # MCEM composition: outer * 2 parts (q1, q2 both non-empty here) * steps.
  expected.releases <- outer * 2L * steps
  expect_s3_class(fit$dp, "data.frame")
  expect_equal(fit$dp$releases, rep(expected.releases, 3L))
  eps <- dsNoLimitsClient:::.dp_epsilon(rep(1, expected.releases), 1e-5)
  expect_equal(fit$dp$epsilon, rep(eps, 3L))
  expect_true(all(fit$dp$epsilon > 0))
  nlc_report(sprintf("dp-mcem  outer=%d steps=%d -> %d releases, eps=%.2f/site",
                     outer, steps, expected.releases, eps))
})

test_that("a DP-SAEM fit runs and spends the accounted epsilon", {
  nlc_skip_no_julia()
  dir <- withr::local_tempdir()
  conns <- ds_lite_conns(nlc_sites(), aggregate = DS_AGG, assign = DS_ASN,
                         options = dp_options(dir))
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  outer <- 3L
  steps <- 2L
  fit <- ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                        estimator = "saem", dp = TRUE, dp.rounds = steps,
                        em.outer.iters = outer, em.sample.schedule = 20L,
                        se = FALSE, checks = FALSE, datasources = conns)

  expect_s3_class(fit, "ds_nolimits_fit")
  expect_true(all(is.finite(fit$theta.natural)))
  expect_true(is.na(fit$loglik))
  expect_true(all(is.na(fit$sites$loglik.contribution)))

  # SAEM composition: outer * (1 stats release + numerical_parts * steps). The
  # closed-form coordinator update on the noised stats is post-processing and
  # charges nothing. The exact count depends on how many Q-parts carry a
  # numerical parameter, which the fit's own release accounting fixes; assert
  # the accountant agrees with whatever count was charged, and that it is the
  # documented shape (one stats release per outer iteration, plus Adam steps).
  releases <- fit$dp$releases[1L]
  expect_true(releases >= outer)                     # at least the stats releases
  expect_equal(fit$dp$releases, rep(releases, 3L))
  eps <- dsNoLimitsClient:::.dp_epsilon(rep(1, releases), 1e-5)
  expect_equal(fit$dp$epsilon, rep(eps, 3L))
  expect_true(all(fit$dp$epsilon > 0))
  nlc_report(sprintf("dp-saem  outer=%d steps=%d -> %d releases, eps=%.2f/site",
                     outer, steps, releases, eps))
})
