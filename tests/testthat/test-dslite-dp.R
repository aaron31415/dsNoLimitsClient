# Tier 4, differential privacy: the whole DP path over three DSLite sites
# running the real server functions.
#
# Utility is NOT gated here. A DP fit is deliberately noised, so its distance
# from the exact fit is a property of the mechanism, not a defect; what is
# asserted is that the fit RUNS, that it spends exactly the privacy budget the
# accountant says it does, and that the budget survives a logout.

test_that("a differentially private fit runs and spends the accounted epsilon", {
  nlc_skip_no_julia()
  dir <- withr::local_tempdir()
  conns <- ds_lite_conns(nlc_sites(), aggregate = DS_AGG, assign = DS_ASN,
                         options = dp_options(dir))
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  rounds <- 10L
  fit <- ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                        estimator = "laplace", dp = TRUE, dp.rounds = rounds,
                        se = FALSE, checks = FALSE, datasources = conns)

  expect_s3_class(fit, "ds_nolimits_fit")
  expect_length(fit$theta.natural, 4L)
  expect_true(all(is.finite(fit$theta.natural)))
  expect_named(fit$theta.natural, c("cl", "v", "omega", "sigma"))
  # No objective value is released under DP, so there is nothing to converge on
  # and nothing to report as a log-likelihood.
  expect_true(is.na(fit$converged))
  expect_true(is.na(fit$loglik))
  expect_true(all(is.na(fit$sites$loglik.contribution)))
  expect_null(fit$se)
  expect_null(fit$vcov.natural)
  # dp.rounds releases plus the one data-free natural-scale round.
  expect_equal(fit$rounds, rounds + 1L)

  expect_s3_class(fit$dp, "data.frame")
  expect_equal(fit$dp$site, c("a", "b", "c"))
  expect_equal(fit$dp$releases, rep(rounds, 3L))
  expect_equal(fit$dp$delta, rep(1e-5, 3L))
  # THE accounting assertion: every site charged exactly `rounds` releases at
  # the owner's noise multiplier, and the client's accountant agrees with the
  # server's to the last bit.
  expected <- dsNoLimitsClient:::.dp_epsilon(rep(1, rounds), 1e-5)
  expect_equal(fit$dp$epsilon, rep(expected, 3L))
  expect_equal(fit$dp$remaining.budget, rep(1e6 - expected, 3L))
  nlc_report(sprintf("dp smoke  sigma=1 T=%d delta=1e-5 -> eps=%.4f per site",
                     rounds, expected))

  expect_output(print(fit), "Differential privacy")

  # The final natural-scale vector came from nolimitsNaturalDS. It must be the
  # SAME quantity the exact path returns inside its own `natural` element, and
  # it must read no data - reusing this session's prep symbol proves both.
  enc <- dsNoLimitsClient:::.ds_num_encode(unname(fit$theta.transformed))
  nat <- DSI::datashield.aggregate(conns, call("nolimitsNaturalDS",
                                               "nolimits.prep", enc))
  ref <- DSI::datashield.aggregate(
    conns, call("nolimitsObjGradDS", "nolimits.prep", enc, "laplace", 1, 1))
  expect_named(nat$a, "natural")
  expect_identical(nat$a$natural, ref$a$natural)
  expect_equal(unname(fit$theta.natural), nat$a$natural)
  # Data-free: every site returns the same vector for the same argument.
  expect_identical(nat$b$natural, nat$a$natural)
  expect_identical(nat$c$natural, nat$a$natural)
})

test_that("logging out and back in does not reset the privacy budget", {
  nlc_skip_no_julia()
  dir <- withr::local_tempdir()
  tables <- nlc_sites("a")

  first <- ds_lite_conns(tables, aggregate = DS_AGG, assign = DS_ASN,
                         options = dp_options(dir))
  fit1 <- ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                         estimator = "laplace", dp = TRUE, dp.rounds = 3L,
                         se = FALSE, checks = FALSE, datasources = first)
  DSI::datashield.logout(first)
  expect_equal(fit1$dp$releases, 3L)

  # A brand-new DSLite server, a brand-new session, the same ledger directory
  # and the same data: the budget picks up where it left off, because it lives
  # on disk and is keyed by the data itself.
  second <- ds_lite_conns(tables, aggregate = DS_AGG, assign = DS_ASN,
                          options = dp_options(dir))
  on.exit(DSI::datashield.logout(second), add = TRUE)
  fit2 <- ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                         estimator = "laplace", dp = TRUE, dp.rounds = 3L,
                         se = FALSE, checks = FALSE, datasources = second)
  expect_equal(fit2$dp$releases, 6L)
  expect_gt(fit2$dp$epsilon, fit1$dp$epsilon)
  expect_lt(fit2$dp$remaining.budget, fit1$dp$remaining.budget)
  expect_equal(fit2$dp$epsilon,
               dsNoLimitsClient:::.dp_epsilon(rep(1, 6), 1e-5))
  # One ledger file for this analysis, holding all six releases.
  expect_length(list.files(dir, pattern = "^dsnolimits-.*csv$"), 1L)
})

test_that("an exhausted budget stops the fit rather than degrading it", {
  nlc_skip_no_julia()
  dir <- withr::local_tempdir()
  # Strictly between the cost of one release and the cost of two, so that the
  # second release is refused whatever rounding the option round trip applies:
  # DSLite pastes the value into a string and parses it back. Under the
  # post-release gate the second release's post-epsilon exceeds the budget.
  budget <- mean(c(dsNoLimitsClient:::.dp_epsilon(1, 1e-5),
                   dsNoLimitsClient:::.dp_epsilon(rep(1, 2), 1e-5)))
  conns <- ds_lite_conns(nlc_sites("a"), aggregate = DS_AGG, assign = DS_ASN,
                         options = dp_options(dir, budget = budget))
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  # A refused release is a server ERROR, and DSI aborts the call on one: the
  # fit stops rather than quietly continuing on a subset of the servers.
  expect_error(
    ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                   estimator = "laplace", dp = TRUE, dp.rounds = 5L,
                   se = FALSE, checks = FALSE, datasources = conns),
    "DataSHIELD errors")
  expect_match(paste(unlist(DSI::datashield.errors()), collapse = " "),
               "budget for this data set is exhausted")
  # Exactly the one release the budget covered was charged, and no more.
  ledger <- list.files(dir, pattern = "^dsnolimits-.*csv$", full.names = TRUE)
  expect_length(ledger, 1L)
  expect_equal(nrow(utils::read.csv(ledger)), 1L)
})

test_that("a server with no ledger directory refuses the private fit", {
  nlc_skip_no_julia()
  # Explicitly EMPTY, not merely absent: DSLite applies server options to the
  # process-wide option list, so a value another test set would still be live.
  conns <- ds_lite_conns(nlc_sites("a"), aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options(
                           dsNoLimits.dp.ledgerDir = .q(""),
                           default.dsNoLimits.dp.ledgerDir = .q("")))
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  expect_error(
    ds.nolimitsFit(NLDS_TEST_MODEL, data = "D", model.mode = "string",
                   estimator = "laplace", dp = TRUE, dp.rounds = 2L,
                   se = FALSE, checks = FALSE, datasources = conns),
    "DataSHIELD errors")
  expect_match(paste(unlist(DSI::datashield.errors()), collapse = " "),
               "not provisioned for differentially private")
})
