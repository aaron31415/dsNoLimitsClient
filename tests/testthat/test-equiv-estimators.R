# Tier 5 (a): the acceptance gate. A 3-site federated fit must be the pooled
# fit, for every estimator whose marginal log-likelihood is a sum of
# independent per-subject terms.
#
# Two things are asserted, and they catch different bugs:
#   * the federated SUM at a point equals the pooled objective at that point
#     (exact to summation order) - this is the transport being additive;
#   * the two optima agree in objective value and in parameters - this is the
#     optimisation driven by that transport landing in the same place.

test_that("federated == pooled for all four estimators on a 3-site split", {
  nlc_skip_no_julia()

  sites <- nlc_sites()
  pooled <- do.call(rbind, unname(sites))
  conns <- ds_lite_conns(sites, aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options())
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  # Which coordinates of theta are comparable between two different optimisers.
  # The naive-pooled objective has an EXACT flat ridge in (cl, omega) for this
  # model: measured on the pooled objective, replacing the reference optimum's
  # cl AND omega by the federated ones leaves the value unchanged in the last
  # bit, while replacing either alone costs 3.9e-3. GHQ at level 1 coincides
  # with the naive-pooled objective here (both optima evaluate to
  # -137.086476936145) and inherits the same ridge. So for those two the
  # objective is the assertion and only v and sigma are comparable iterates;
  # for Laplace and FOCEI the whole vector is.
  identified <- list(laplace = 1:4, focei = 1:4, ghq = c(2L, 4L), pooled = c(2L, 4L))

  for (est in c("laplace", "focei", "ghq", "pooled")) {
    fed <- ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                          estimator = est, ghq.level = 1L, se = FALSE,
                          checks = FALSE, datasources = conns)
    ref <- nlc_pooled(pooled, est, 1L, key = "csv3")

    # 1. additivity of the transport: the sum of the three sites' values at the
    #    federated optimum IS the pooled objective there.
    at.fed <- nlc_objgrad(ref$dm, fed$theta.transformed, est)$value
    expect_equal(fed$loglik, at.fed, tolerance = 1e-9)

    # 2. per-site contributions sum to the total.
    expect_equal(sum(fed$sites$loglik.contribution), fed$loglik, tolerance = 1e-12)
    expect_equal(sum(fed$sites$n.subjects), 18L)
    expect_equal(sum(fed$sites$n.obs), nrow(pooled))

    # 3. objective agreement at the two optima (relative), and parameters.
    obj.rel <- nlc_reldiff(at.fed, ref$value)
    dtheta <- abs(fed$theta.transformed - ref$theta.t)
    par.max <- max(dtheta[identified[[est]]])
    nlc_report(sprintf("%-8s obj.rel=%.3e  max|dtheta|=%.3e (all %.3e)  rounds=%d",
                       est, obj.rel, par.max, max(dtheta), fed$rounds))
    expect_lt(obj.rel, 1e-6)
    expect_lt(par.max, 1e-3)

    expect_true(fed$converged)
    expect_equal(names(fed$theta.transformed), c("cl", "v", "omega", "sigma"))
  }
})
