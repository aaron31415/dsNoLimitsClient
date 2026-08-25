## Every check here must fire BEFORE a byte leaves the machine. The fake
## connections pass .setDatasources() but cannot serve a request, so a test that
## sees the expected validation message has proven no request was attempted.

conns <- fake_datasources(2L)

test_that(".setDatasources rejects anything that is not a connection list", {
  expect_error(dsNoLimitsClient:::.setDatasources(list()), "DSConnection")
  expect_error(dsNoLimitsClient:::.setDatasources(list(1, 2)), "DSConnection")
  expect_error(dsNoLimitsClient:::.setDatasources("site1"), "DSConnection")
  expect_error(dsNoLimitsClient:::.setDatasources(unname(conns)), "DSConnection")
  expect_identical(dsNoLimitsClient:::.setDatasources(conns), conns)
})

test_that("estimator and model.mode use match.arg semantics", {
  expect_error(
    ds.nolimitsFit("m", data = "D", estimator = "nosuch", datasources = conns),
    "should be one of"
  )
  expect_error(
    ds.nolimitsFit("m", data = "D", model.mode = "inline", datasources = conns),
    "should be one of"
  )
  # A partial match is accepted by match.arg and must get past validation.
  expect_error(
    ds.nolimitsFit(NULL, data = "D", estimator = "foc", datasources = conns),
    "model"
  )
  # The no-random-effects estimators are in the set, and "m" is ambiguous
  # between them, which is what match.arg must say.
  for (est in c("mle", "map")) {
    expect_error(
      ds.nolimitsFit(NULL, data = "D", estimator = est, datasources = conns),
      "model"
    )
  }
  expect_error(
    ds.nolimitsFit("m", data = "D", estimator = "m", datasources = conns),
    "should be one of"
  )
})

test_that("model is validated before any call", {
  expect_error(ds.nolimitsFit(NULL, data = "D", datasources = conns), "model")
  expect_error(ds.nolimitsFit(c("a", "b"), data = "D", datasources = conns), "model")
  expect_error(ds.nolimitsFit("", data = "D", datasources = conns), "model")
  expect_error(ds.nolimitsFit(NA_character_, data = "D", datasources = conns), "model")
  # Julia source in registered mode is a client-side error, not a server one.
  expect_error(
    ds.nolimitsFit("conc ~ Normal(mu, sigma)", data = "D", datasources = conns),
    "registered model name"
  )
})

test_that("symbol-shaped arguments are validated before any call", {
  expect_error(ds.nolimitsFit("m", data = NULL, datasources = conns), "data")
  expect_error(ds.nolimitsFit("m", data = "D$x", datasources = conns), "data")
  expect_error(ds.nolimitsFit("m", data = "2D", datasources = conns), "data")
  expect_error(ds.nolimitsFit("m", data = "D", id.col = "id col", datasources = conns), "id.col")
  expect_error(ds.nolimitsFit("m", data = "D", time.col = NA, datasources = conns), "time.col")
  expect_error(ds.nolimitsFit("m", data = "D", newobj = "", datasources = conns), "newobj")
})

test_that("numeric and flag arguments are validated before any call", {
  expect_error(ds.nolimitsFit("m", data = "D", ghq.level = 0, datasources = conns), "ghq.level")
  expect_error(ds.nolimitsFit("m", data = "D", ghq.level = 10, datasources = conns), "ghq.level")
  expect_error(ds.nolimitsFit("m", data = "D", ghq.level = 1.5, datasources = conns), "ghq.level")
  expect_error(ds.nolimitsFit("m", data = "D", ghq.level = "1", datasources = conns), "ghq.level")
  expect_error(ds.nolimitsFit("m", data = "D", max.rounds = 0, datasources = conns), "max.rounds")
  expect_error(ds.nolimitsFit("m", data = "D", se = NA, datasources = conns), "se")
  expect_error(ds.nolimitsFit("m", data = "D", checks = "yes", datasources = conns), "checks")
})

test_that("a bad datasources argument is caught in every entry point", {
  expect_error(ds.nolimitsFit("m", data = "D", datasources = list()), "DSConnection")
  expect_error(ds.nolimitsStatus(datasources = list(1)), "DSConnection")
  expect_error(ds.nolimitsModels(datasources = "nope"), "DSConnection")
})

test_that("agree() names both disagreeing sites", {
  base <- list(names = c("v", "cli"), theta0 = c(1, 2), scale = c(1, 1))
  same <- list(alpha = base, bravo = base, charlie = base)
  expect_identical(agree(same), base)

  drift <- same
  drift$charlie$theta0 <- c(1, 2.0000001)
  err <- expect_error(agree(drift))
  msg <- conditionMessage(err)
  expect_match(msg, "alpha")
  expect_match(msg, "charlie")
  expect_match(msg, "not running the same model")
  expect_match(msg, "theta0")

  drift2 <- same
  drift2$bravo$names <- c("v", "CL")
  expect_error(agree(drift2), "not running the same model")

  drift3 <- same
  drift3$bravo$scale <- c(1, 2)
  expect_error(agree(drift3), "not running the same model")
})

test_that(".agree_model passes on matching hashes and names both on mismatch", {
  h1 <- strrep("a", 64L)
  h2 <- strrep("b", 64L)
  same <- list(alpha = list(model.hash = h1),
               bravo = list(model.hash = h1),
               charlie = list(model.hash = h1))
  expect_identical(.agree_model(same), h1)

  drift <- list(alpha = list(model.hash = h1), bravo = list(model.hash = h2))
  err <- expect_error(.agree_model(drift), "model hash mismatch")
  msg <- conditionMessage(err)
  expect_match(msg, "alpha")
  expect_match(msg, "bravo")
  # Both hashes in short form (first 12 hex).
  expect_match(msg, substr(h1, 1L, 12L))
  expect_match(msg, substr(h2, 1L, 12L))

  # A missing hash (an old server) aborts with the provisioning hint.
  none <- list(alpha = list(theta0 = 1), bravo = list(theta0 = 1))
  expect_error(.agree_model(none), "no model hash")
})

test_that("two real byte-different same-name models trip the agreement check", {
  # DSLite serves every site from ONE R process with ONE global options(), so two
  # sites cannot hold different registered files under the same name (the last
  # server's dsNoLimits.modelDir wins globally). The mismatch is therefore
  # exercised with the REAL server hash of two byte-different sources that differ
  # only by a numerically-inert `+ 0.0` - agree() on names/theta0/scale would
  # pass, and only the canonical model hash separates them.
  skip_if_not_installed("dsNoLimits")
  h <- function(src) dsNoLimits:::.nlds_model_hash(src)

  model.b <- sub("* t), sigma)", "* t) + 0.0, sigma)", NLDS_TEST_MODEL,
                 fixed = TRUE)
  ha <- h(NLDS_TEST_MODEL)
  hb <- h(model.b)
  expect_false(identical(ha, hb))
  # A comment-only edit is the SAME registered model.
  expect_identical(ha, h(paste0(NLDS_TEST_MODEL, "\n# owner note")))

  # The fit's agreement check on the two sites' infos aborts, naming both.
  infos <- list(a = list(model.hash = ha), b = list(model.hash = hb))
  err <- expect_error(.agree_model(infos), "model hash mismatch")
  msg <- conditionMessage(err)
  expect_match(msg, "a")
  expect_match(msg, "b")
  expect_match(msg, substr(ha, 1L, 12L))
  expect_match(msg, substr(hb, 1L, 12L))
})

test_that("agree() tolerates BLAS-level scale drift but not real disagreement", {
  base <- list(names = c("v", "cli"), theta0 = c(1, 2), scale = c(1, 1))
  # A scale that differs only at 1e-12 (heterogeneous BLAS) agrees, and site 1's
  # scale is the one returned.
  tiny <- base; tiny$scale <- c(1, 1 + 1e-12)
  ref <- agree(list(alpha = base, bravo = tiny))
  expect_identical(ref$scale, base$scale)
  # A scale that differs at 1e-3 is a real disagreement and still aborts.
  big <- base; big$scale <- c(1, 1 + 1e-3)
  err <- expect_error(agree(list(alpha = base, bravo = big)),
                      "not running the same model")
  expect_match(conditionMessage(err), "scale")
  # names and theta0 remain EXACT: a 1e-7 theta0 drift still aborts.
  th <- base; th$theta0 <- c(1, 2 + 1e-7)
  expect_error(agree(list(alpha = base, bravo = th)), "not running the same model")
})

test_that(".check_complete names the absent sites", {
  expect_error(
    dsNoLimitsClient:::.check_complete(list(a = 1), c("a", "b"), "status"),
    "b"
  )
  expect_silent(dsNoLimitsClient:::.check_complete(list(a = 1, b = 2), c("a", "b"), "status"))
})

test_that("confint validates its own arguments", {
  fit <- structure(
    list(theta.natural = c(v = 1, cl = 2),
         theta.transformed = c(v = 0, cl = 0.7),
         vcov.natural = matrix(c(0.04, 0, 0, 0.09), 2, 2,
                               dimnames = list(c("v", "cl"), c("v", "cl"))),
         vcov.transformed = NULL),
    class = "ds_nolimits_fit"
  )
  expect_error(confint(fit, level = 1.5), "level")
  expect_error(confint(fit, parm = "nope"), "unknown parameter")
  expect_error(confint(fit, scale = "transformed"), "se = TRUE")

  ci <- confint(fit)
  expect_equal(ci[, 1], c(v = 1, cl = 2) - stats::qnorm(0.975) * c(0.2, 0.3))
  expect_identical(colnames(ci), c("2.5 %", "97.5 %"))
})
