# Tier 5 for SAEM, the nested federated EM estimator with a hybrid M-step.
#
# Layer A (the federated-learning proof, machine precision): the DE-NORMALIZED
# per-subject-additive sufficient statistics sum coordinate-wise to the population
# statistics. This is the proof that the client's element-wise sum of the per-site
# flat payloads IS the pooled sufficient statistics the coordinator re-normalizes.
#
# Layer B (Monte-Carlo): a 3-site DSLite federated SAEM fit reproduces the pooled
# fit_model(dm, SAEM()) at the same seed and outer budget. Loose by design (the
# per-site RNG partition differs from the pooled one); the exactness is Layer A.
# The 3-site split exercises the coordinator path (>= 2 sites, so a non-coordinator
# site exists).

test_that("SAEM sufficient statistics are per-subject additive (Layer A)", {
  nlc_skip_no_julia()
  d <- do.call(rbind, unname(nlc_sites()))
  dm <- nlc_dm(d, NLDS_TEST_MODEL)
  dsNoLimits:::.nlds_helpers()
  sf <- function(n) dsNoLimits:::.nlds_fun(n)

  theta0 <- as.numeric(sf("nlds_theta0")(dm))
  cf <- as.character(sf("nlds_saem_cf")(dm))
  num <- as.character(sf("nlds_saem_num")(dm))
  expect_setequal(c(cf, num), as.character(sf("nlds_names")(dm)))
  expect_true(length(intersect(cf, num)) == 0L)
  expect_true("omega" %in% cf)

  # ONE E-step (state = nothing) at theta0 gives the shared fixed draws; each batch
  # is one subject, so summing per-idx equals summing whole sites.
  em <- sf("nlds_mcem_estep")(dm, theta0, 100, as.integer(20260824L), NULL)
  nb <- as.integer(sf("nlds_mcem_nbatches")(em))
  expect_equal(nb, length(unique(d$ID)))

  pooled <- as.numeric(sf("nlds_saem_stats_flat")(dm, theta0, em))
  fed <- rep(0, length(pooled))
  for (i in seq_len(nb)) {
    fed <- fed + as.numeric(sf("nlds_saem_stats_flat_idx")(dm, theta0, em,
                                                           as.integer(i)))
  }
  rel <- max(abs(fed - pooled) / pmax(abs(pooled), 1e-300))
  nlc_report(sprintf("saem Layer A stats additivity rel = %.3e (flat len %d, %d subjects)",
                     rel, length(pooled), nb))
  expect_lt(rel, 1e-10)
})

test_that("federated SAEM reproduces the pooled fit (Layer B)", {
  nlc_skip_no_julia()

  sites <- nlc_sites()
  pooled <- do.call(rbind, unname(sites))
  conns <- ds_lite_conns(sites, aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options())
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  fed <- ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                        estimator = "saem", checks = FALSE, datasources = conns,
                        em.outer.iters = 20L, em.sample.schedule = 100L,
                        em.seed = 20260824L)

  # A stochastic, fixed-budget estimator that exercised the coordinator path.
  expect_identical(fed$estimator, "saem")
  expect_identical(fed$outer.iters, 20L)
  expect_identical(fed$coordinator, sort(names(conns))[1L])
  expect_true(fed$coordinator %in% names(conns))
  expect_gt(length(names(conns)), 1L)
  expect_true(is.na(fed$loglik))
  expect_true(is.na(fed$converged))
  expect_null(fed$se)
  expect_true(all(is.na(fed$sites$loglik.contribution)))
  expect_equal(names(fed$theta.natural), c("cl", "v", "omega", "sigma"))

  ref <- nlc_saem_pooled(pooled, 100L, 20L, 20260824L, key = "csv3")
  fn <- names(ref$natural)
  rel <- abs(fed$theta.natural[fn] - ref$natural[fn]) / pmax(abs(ref$natural[fn]), 1e-8)
  nlc_report(sprintf("saem Layer B: worst param rel=%.3e (%s), rounds=%d",
                     max(rel), fn[which.max(rel)], fed$rounds))
  # Monte-Carlo tolerance, matching the reference federated application.
  expect_lt(max(rel), 5e-2)
})
