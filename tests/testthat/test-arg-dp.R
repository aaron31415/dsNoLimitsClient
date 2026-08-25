## The differentially private client path, driven with a stub instead of DSI.

dp_theta0 <- c(0.5, -1, 2)
dp_s <- c(1, 2, 0.5)
dp_sites <- c("alpha", "bravo")

# A stub site returning a fixed preconditioned gradient plus a noise draw, and
# recording every call it was asked to make.
dp_stub <- function(counter, grads, noisy = FALSE, p = 3L) {
  function(cally) {
    counter$n <- counter$n + 1L
    counter$calls <- c(counter$calls, list(cally))
    theta <- dsNoLimitsClient:::.ds_num_decode(cally[[3L]])
    counter$thetas <- c(counter$thetas, list(theta))
    lapply(stats::setNames(seq_along(dp_sites), dp_sites), function(k) {
      g <- grads[[k]](theta)
      if (noisy) g <- g + stats::rnorm(p, 0, 1e-6)
      list(grad = g, releases = counter$n, epsilon = 0.5 * counter$n,
           delta = 1e-5, remaining.budget = 10 - 0.5 * counter$n,
           n.subjects = 10L)
    })
  }
}

test_that("the accountant copies in the two packages agree", {
  skip_if_not_installed("dsNoLimits")
  cases <- list(list(1, 1e-5), list(rep(1, 50), 1e-5), list(rep(4, 50), 1e-5),
                list(4, 1e-6), list(c(1, 3, 0.7, 2, 5), 1e-4))
  for (cs in cases) {
    expect_identical(dsNoLimitsClient:::.dp_epsilon(cs[[1L]], cs[[2L]]),
                     dsNoLimits:::.dp_epsilon(cs[[1L]], cs[[2L]]))
  }
})

test_that("the client accountant matches the closed forms", {
  expect_equal(dsNoLimitsClient:::.dp_epsilon(1, 1e-5), 5.2985, tolerance = 1e-3)
  expect_equal(dsNoLimitsClient:::.dp_epsilon(rep(1, 50), 1e-5), 58.9308,
               tolerance = 1e-3)
  expect_equal(dsNoLimitsClient:::.dp_epsilon(rep(4, 50), 1e-5), 10.0452,
               tolerance = 1e-3)
  expect_equal(dsNoLimitsClient:::.dp_epsilon(4, 1e-6), 1.3454, tolerance = 1e-3)
  expect_equal(dsNoLimitsClient:::.dp_epsilon(numeric(0), 1e-5), 0)
})

test_that("the private round sums the sites and never rescales the gradient", {
  counter <- new.env(); counter$n <- 0L; counter$calls <- list()
  grads <- list(function(th) c(1, 2, 3), function(th) c(0.5, -1, 0))
  spec <- list(newobj = "prep", estimator = "laplace", ghq.level = 1L, p = 3L,
               sites = dp_sites, aggregate_fun = dp_stub(counter, grads))
  r <- dsNoLimitsClient:::.dp_round(dp_theta0, 1L, spec)
  # A PLAIN unweighted sum: the servers already clipped and noised in the
  # preconditioned coordinates, so nothing is multiplied by s here.
  expect_equal(r$grad, c(1.5, 1, 3))
  expect_named(r$per.site, dp_sites)
  expect_identical(counter$n, 1L)
})

test_that("a site returning an unusable private gradient aborts the fit", {
  counter <- new.env(); counter$n <- 0L; counter$calls <- list()
  spec <- list(newobj = "prep", estimator = "laplace", ghq.level = 1L, p = 3L,
               sites = dp_sites,
               aggregate_fun = dp_stub(counter,
                                       list(function(th) c(1, 2, 3),
                                            function(th) c(1, NaN, 3))))
  expect_error(dsNoLimitsClient:::.dp_round(dp_theta0, 1L, spec),
               "unusable private gradient")
  # A missing site is an absent element, never a NULL one.
  spec$sites <- c(dp_sites, "charlie")
  expect_error(dsNoLimitsClient:::.dp_round(dp_theta0, 1L, spec), "charlie")
})

test_that("Adam spends exactly one release per round and climbs the objective", {
  counter <- new.env(); counter$n <- 0L; counter$calls <- list()
  counter$thetas <- list()
  # Two quadratic sites: the summed gradient in z vanishes at one point.
  mu <- list(alpha = c(1, 1, 1), bravo = c(-1, 3, 0))
  grads <- lapply(mu, function(m) function(th) dp_s * (-(th - m)))
  spec <- list(newobj = "prep", estimator = "laplace", ghq.level = 1L, p = 3L,
               sites = dp_sites, aggregate_fun = dp_stub(counter, grads))
  run <- dsNoLimitsClient:::.dp_optimize(3L, dp_theta0, dp_s, 0.05, 40L, spec)
  expect_identical(counter$n, 40L)
  expect_identical(run$rounds, 40L)
  # The summed optimum is the average of the two targets; Adam ascends towards
  # it, which is only true if the sign of the step is right.
  target <- (mu$alpha + mu$bravo) / 2
  theta <- dp_theta0 + dp_s * run$z
  expect_lt(max(abs(theta - target)), max(abs(dp_theta0 - target)))
  # Every call is the DP method, with the round number counting up from 1.
  expect_true(all(vapply(counter$calls,
                         function(cl) identical(as.character(cl[[1L]]),
                                                "nolimitsObjGradDpDS"),
                         logical(1))))
  expect_equal(vapply(counter$calls, function(cl) cl[[6L]], numeric(1)),
               as.numeric(seq_len(40L)))
})

test_that("the DP report is one row per site of the servers' own numbers", {
  per.site <- list(
    alpha = list(releases = 12L, epsilon = 3.5, delta = 1e-5,
                 remaining.budget = 6.5),
    bravo = list(releases = 4L, epsilon = 1.1, delta = 1e-6,
                 remaining.budget = 8.9))
  rep <- dsNoLimitsClient:::.dp_report(per.site)
  expect_equal(rep$site, c("alpha", "bravo"))
  expect_equal(rep$releases, c(12L, 4L))
  expect_equal(rep$epsilon, c(3.5, 1.1))
  expect_equal(rep$remaining.budget, c(6.5, 8.9))
})

test_that("dp arguments are validated before any network call", {
  expect_error(ds.nolimitsFit("m", data = "D", dp = "yes"),
               "must be .*TRUE.* or .*FALSE")
  for (est in c("pooled", "map")) {
    expect_error(ds.nolimitsFit("m", data = "D", dp = TRUE, estimator = est),
                 "no per-subject term")
  }
  expect_error(ds.nolimitsFit("m", data = "D", dp = TRUE, dp.rounds = 0),
               "dp.rounds")
  expect_error(ds.nolimitsFit("m", data = "D", dp = TRUE, dp.lr = -1),
               "dp.lr")
  expect_error(ds.nolimitsFit("m", data = "D", dp = TRUE, dp.lr = c(1, 2)),
               "dp.lr")
  # se is not an error under dp, it is a message and a silent downgrade; the
  # abort that follows is the missing connection list, which proves the
  # validation ran first.
  expect_message(
    try(ds.nolimitsFit("m", data = "D", dp = TRUE, se = TRUE,
                       datasources = list()), silent = TRUE),
    "is ignored under")
})

test_that("both DP call shapes are legal on both wire grammars", {
  theta <- dsNoLimitsClient:::.ds_num_encode(
    c(0.5, -1, 2, 1e-300, .Machine$double.xmax))
  for (est in c("laplace", "focei", "ghq", "mle")) {
    for (lvl in 1:9) {
      cally <- call("nolimitsObjGradDpDS", "nolimits.prep", theta, est,
                    as.numeric(lvl), 37)
      expect_silent(ds_wire_check(cally, "aggregate", DSNOLIMITS_AGGREGATE))
    }
  }
  expect_silent(ds_wire_check(call("nolimitsNaturalDS", "nolimits.prep", theta),
                              "aggregate", DSNOLIMITS_AGGREGATE))
})
