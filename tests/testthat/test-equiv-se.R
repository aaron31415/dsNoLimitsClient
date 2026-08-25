# Tier 5 (e) and (f): the round budget and the Wald standard errors, plus one
# full fit through model.mode = "registered".

test_that("Wald SEs match a local FD Hessian of the pooled objective", {
  nlc_skip_no_julia()

  sites <- nlc_sites()
  pooled <- do.call(rbind, unname(sites))
  conns <- ds_lite_conns(sites, aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options())
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  fed <- ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                        estimator = "laplace", se = TRUE, checks = FALSE,
                        datasources = conns)
  ref <- nlc_pooled(pooled, "laplace", key = "csv3")

  ## (e) round count. The budget is 1 assign + 1 prep info + optimiser rounds
  ## + 2p SE probes; only the aggregate rounds are counted by fit$rounds.
  nlc_report(sprintf("laplace+se rounds=%d", fed$rounds))
  expect_lt(fed$rounds, 60L)
  # Pinned to the observed count as a regression tripwire: a change here is the
  # earliest symptom of a broken memo (optim calls fn and gr separately at the
  # same point), a changed start value or a changed optimiser control list.
  expect_equal(fed$rounds, NLC_EXPECTED_ROUNDS)

  ## (f) the same central-difference formula, computed locally on the pooled
  ## objective at the SAME point the federated run used.
  theta <- fed$theta.transformed
  p <- length(theta)
  h <- .Machine$double.eps^(1 / 3) * pmax(abs(theta), 1)
  G <- matrix(NA_real_, p, p)
  for (j in seq_len(p)) {
    tp <- theta; tp[j] <- tp[j] + h[j]
    tm <- theta; tm[j] <- tm[j] - h[j]
    G[, j] <- (nlc_objgrad(ref$dm, tp)$grad - nlc_objgrad(ref$dm, tm)$grad) / (2 * h[j])
  }
  H <- (G + t(G)) / 2
  se.local <- sqrt(diag(solve(-H)))

  se.fed <- sqrt(diag(fed$vcov.transformed))
  rel <- max(abs(se.fed - se.local) / se.local)
  nlc_report(sprintf("se.rel=%.3e  se=%s", rel,
                     paste(sprintf("%.5f", se.fed), collapse = ",")))
  expect_lt(rel, 1e-4)

  expect_true(all(is.finite(fed$se)))
  expect_equal(names(fed$se), c("cl", "v", "omega", "sigma"))
  expect_equal(dim(fed$vcov.natural), c(4L, 4L))
  # confint() is the analyst-facing consumer of the same matrix.
  ci <- stats::confint(fed)
  expect_true(all(ci[, 1L] < coef(fed) & coef(fed) < ci[, 2L]))
})

test_that("a registered-mode fit equals the string-mode fit", {
  nlc_skip_no_julia()

  sites <- nlc_sites()
  opts <- ds_model_dir_options(ds_test_options())
  conns <- ds_lite_conns(sites, aggregate = DS_AGG, assign = DS_ASN, options = opts)
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  # checks = FALSE: DSLite reports a package name for a method only when the
  # config was built by defaultDSConfiguration(), and this harness registers
  # methods explicitly. The checks = TRUE path runs in test-dslite-manifest.R,
  # which uses the real manifest.
  fed <- ds.nolimitsFit("onecomp_iv", data = "D", model.mode = "registered",
                        estimator = "laplace", se = FALSE, checks = FALSE,
                        datasources = conns)
  ref <- nlc_pooled(do.call(rbind, unname(sites)), "laplace", key = "csv3")

  at.fed <- nlc_objgrad(ref$dm, fed$theta.transformed)$value
  expect_equal(fed$loglik, at.fed, tolerance = 1e-9)
  expect_lt(nlc_reldiff(at.fed, ref$value), 1e-6)
  expect_lt(max(abs(fed$theta.transformed - ref$theta.t)), 1e-3)

  # The fit carries the agreed model hash, equal to the canonical hash of the
  # registered source, and it prints in short form.
  expect_match(fed$model.hash, "^[0-9a-f]{64}$")
  expect_identical(fed$model.hash, dsNoLimits:::.nlds_model_hash(NLDS_TEST_MODEL))
  expect_output(print(fed), substr(fed$model.hash, 1L, 12L), fixed = TRUE)
})
