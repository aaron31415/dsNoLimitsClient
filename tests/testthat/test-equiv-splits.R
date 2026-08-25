# Tier 5 (b), (c), (d): the split-shape tests. Laplace only - a pooling bug is
# a property of the transport, not of the estimator, and the estimator matrix
# is covered in test-equiv-estimators.R.

test_that("a one-site federation is the plain pooled fit on that site", {
  nlc_skip_no_julia()

  a <- nlc_site("a")
  conns <- ds_lite_conns(list(a = a), aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options())
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  fed <- ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                        estimator = "laplace", se = FALSE, checks = FALSE,
                        datasources = conns)
  ref <- nlc_pooled(a, "laplace", key = "site.a")

  # With one site there is nothing to pool, so any pooling formula that is not
  # an identity fails right here.
  expect_equal(fed$loglik, nlc_objgrad(ref$dm, fed$theta.transformed)$value,
               tolerance = 1e-12)
  obj.rel <- nlc_reldiff(nlc_objgrad(ref$dm, fed$theta.transformed)$value, ref$value)
  par.max <- max(abs(fed$theta.transformed - ref$theta.t))
  nlc_report(sprintf("1-site   obj.rel=%.3e  max|dtheta|=%.3e", obj.rel, par.max))
  expect_lt(obj.rel, 1e-6)
  expect_lt(par.max, 1e-3)
  expect_equal(nrow(fed$sites), 1L)
})

test_that("additivity: {a, b} federated equals the pooled fit on rbind(a, b)", {
  nlc_skip_no_julia()

  a <- nlc_site("a")
  b <- nlc_site("b")
  conns <- ds_lite_conns(list(a = a, b = b), aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options())
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  fed <- ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                        estimator = "laplace", se = FALSE, checks = FALSE,
                        datasources = conns)
  ref <- nlc_pooled(rbind(a, b), "laplace", key = "site.ab")

  # The only test that catches a dropped between-site term.
  at.fed <- nlc_objgrad(ref$dm, fed$theta.transformed)$value
  expect_equal(fed$loglik, at.fed, tolerance = 1e-9)
  obj.rel <- nlc_reldiff(at.fed, ref$value)
  par.max <- max(abs(fed$theta.transformed - ref$theta.t))
  nlc_report(sprintf("2-site   obj.rel=%.3e  max|dtheta|=%.3e", obj.rel, par.max))
  expect_lt(obj.rel, 1e-6)
  expect_lt(par.max, 1e-3)

  # And each site's contribution is its own fit's value at the shared optimum.
  a.only <- nlc_objgrad(nlc_pooled(a, "laplace", key = "site.a")$dm,
                        fed$theta.transformed)$value
  expect_equal(fed$sites$loglik.contribution[fed$sites$site == "a"], a.only,
               tolerance = 1e-9)
})

test_that("an unequal 5/10/15 split passes at the canonical nfilter.glm", {
  nlc_skip_no_julia()

  d <- nlc_sim_data()
  sp <- nlc_sim_split(d)
  # The smallest site is 5 subjects, so the saturation gate needs 4 <= glm * 5:
  # 0.8 is exactly at the threshold and the site is legal but tiny.
  conns <- ds_lite_conns(sp, aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options(nfilter.glm = 0.8))
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  fed <- ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                        estimator = "laplace", se = FALSE, checks = FALSE,
                        datasources = conns)
  ref <- nlc_pooled(d, "laplace", key = "sim30")

  at.fed <- nlc_objgrad(ref$dm, fed$theta.transformed)$value
  expect_equal(fed$loglik, at.fed, tolerance = 1e-9)
  obj.rel <- nlc_reldiff(at.fed, ref$value)
  par.max <- max(abs(fed$theta.transformed - ref$theta.t))
  nlc_report(sprintf("unequal  obj.rel=%.3e  max|dtheta|=%.3e", obj.rel, par.max))
  expect_lt(obj.rel, 1e-6)
  expect_lt(par.max, 1e-3)
  expect_equal(fed$sites$n.subjects, c(5L, 10L, 15L))
})

test_that("a generated 3-way split preps at the canonical nfilter.glm = 0.33", {
  nlc_skip_no_julia()

  # 45 subjects, 15 per site: the canonical default admits 0.33 * 15 = 4.95
  # parameters and the model has 4. The 6-subject CSV fixtures cannot clear
  # this gate at any split, which is why the generator exists.
  sp <- nlc_sim_split(nlc_sim_data(n = 45L), sizes = c(15L, 15L, 15L))
  conns <- ds_lite_conns(sp, aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options(nfilter.glm = 0.33))
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  DSI::datashield.assign(
    conns, "nolimits.prep",
    call("nolimitsPrepStringDS", "D",
         dsNoLimitsClient:::.ds_encode(NLDS_TEST_MODEL), "ID", "t"))
  expect_setequal(DSI::datashield.symbols(conns)$u1, c("D", "nolimits.prep"))

  info <- DSI::datashield.aggregate(conns, call("nolimitsPrepInfoDS", "nolimits.prep"))
  expect_equal(info$u1$n.subjects, 15L)
  expect_equal(info$u3$n.subjects, 15L)
})
