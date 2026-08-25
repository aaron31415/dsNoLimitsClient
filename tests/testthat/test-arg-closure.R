## The federated closure is factored out of ds.nolimitsFit so it can be driven
## with a counting stub instead of DSI. optim() asks for fn and gr separately at
## the same point; if the memo ever breaks, every fit silently doubles its
## number of round trips.

theta0 <- c(0.5, -1, 2)
s <- c(1, 2, 0.5)
sites <- c("alpha", "bravo")

# A stub site: value = -sum((theta - target)^2)/2, gradient = -(theta - target),
# so the summed federated objective has a unique known optimum.
stub_agg <- function(targets, counter, finite = TRUE) {
  function(cally) {
    counter$n <- counter$n + 1L
    counter$calls <- c(counter$calls, list(cally))
    theta <- dsNoLimitsClient:::.ds_num_decode(cally[[3L]])
    lapply(targets, function(mu) {
      list(value = -sum((theta - mu)^2) / 2,
           grad = -(theta - mu),
           natural = exp(theta),
           finite = finite,
           n.subjects = 10L)
    })
  }
}

targets <- list(alpha = c(1, 1, 1), bravo = c(-1, 3, 0))

test_that("fn and gr at the same z cost ONE federated round", {
  counter <- new.env(); counter$n <- 0L; counter$calls <- list()
  fed <- dsNoLimitsClient:::.fed_closure(theta0, s, "prep", "laplace", 1L, sites,
                                         stub_agg(targets, counter))
  z <- c(0.1, 0.2, 0.3)
  a <- fed$eval(z)$fn
  b <- fed$eval(z)$gr
  expect_identical(counter$n, 1L)
  expect_identical(fed$rounds(), 1L)

  fed$eval(z + 1e-8)
  expect_identical(counter$n, 2L)
  fed$eval(z)                        # memo holds only the most recent point
  expect_identical(counter$n, 3L)
  expect_type(a, "double")
  expect_length(b, 3L)
})

test_that("the closure sums sites and applies the chain rule and negation", {
  counter <- new.env(); counter$n <- 0L; counter$calls <- list()
  fed <- dsNoLimitsClient:::.fed_closure(theta0, s, "prep", "laplace", 1L, sites,
                                         stub_agg(targets, counter))
  z <- c(0.1, -0.2, 0.3)
  theta <- theta0 + s * z
  r <- fed$eval(z)

  value <- sum(vapply(targets, function(mu) -sum((theta - mu)^2) / 2, numeric(1)))
  grad <- Reduce(`+`, lapply(targets, function(mu) -(theta - mu)))
  expect_equal(r$value, value)
  expect_equal(r$fn, -value)
  expect_equal(r$gr, -(s * grad))
  expect_true(r$finite)
})

test_that("a non-finite site makes L-BFGS-B backtrack instead of dying", {
  counter <- new.env(); counter$n <- 0L; counter$calls <- list()
  fed <- dsNoLimitsClient:::.fed_closure(theta0, s, "prep", "laplace", 1L, sites,
                                         stub_agg(targets, counter, finite = FALSE))
  r <- fed$eval(c(0, 0, 0))
  expect_false(r$finite)
  expect_identical(r$fn, 1e12)
  expect_identical(r$gr, rep(0, 3L))
})

test_that("an absent site aborts rather than pooling k-1 sites", {
  counter <- new.env(); counter$n <- 0L; counter$calls <- list()
  short <- function(cally) stub_agg(targets, counter)(cally)["alpha"]
  fed <- dsNoLimitsClient:::.fed_closure(theta0, s, "prep", "laplace", 1L, sites, short)
  expect_error(fed$eval(c(0, 0, 0)), "bravo")
})

test_that("the emitted objgrad call carries an incrementing round counter", {
  counter <- new.env(); counter$n <- 0L; counter$calls <- list()
  fed <- dsNoLimitsClient:::.fed_closure(theta0, s, "myprep", "ghq", 3L, sites,
                                         stub_agg(targets, counter))
  fed$eval(c(0, 0, 0))
  fed$probe(theta0)
  expect_identical(counter$n, 2L)
  first <- counter$calls[[1L]]
  expect_identical(as.character(first[[1L]]), "nolimitsObjGradDS")
  expect_identical(first[[2L]], "myprep")
  expect_identical(first[[4L]], "ghq")
  expect_identical(first[[5L]], 3)
  expect_identical(first[[6L]], 1)
  expect_identical(counter$calls[[2L]][[6L]], 2)
  # Never 3L: a deparsed "3L" is not guaranteed to lex on Opal.
  expect_false(is.integer(first[[5L]]))
})

test_that("optim drives the closure to the analytic optimum in few rounds", {
  counter <- new.env(); counter$n <- 0L; counter$calls <- list()
  fed <- dsNoLimitsClient:::.fed_closure(theta0, s, "prep", "laplace", 1L, sites,
                                         stub_agg(targets, counter))
  opt <- stats::optim(rep(0, 3), fn = function(z) fed$eval(z)$fn,
                      gr = function(z) fed$eval(z)$gr, method = "L-BFGS-B",
                      control = list(maxit = 100, factr = 1e7, pgtol = 1e-5))
  expect_identical(opt$convergence, 0L)
  expect_equal(theta0 + s * opt$par, (targets$alpha + targets$bravo) / 2,
               tolerance = 1e-6)
  # One round per evaluation, not two: this is the memo doing its job.
  expect_identical(counter$n, fed$rounds())
  expect_lt(counter$n, 30L)
})

## ----------------------------------------------------------------- map ------
## The prior-once combination. A plain sum over the two stub sites would count
## the log-prior twice, which is exactly what these tests fail on.

stub_map_agg <- function(targets, counter, prior.value = -2.5,
                         prior.grad = c(0.1, -0.2, 0.3), bad = character(0)) {
  base <- stub_agg(targets, counter)
  function(cally) {
    res <- base(cally)
    for (nm in names(res)) {
      bump <- if (nm %in% bad) 1e-12 else 0
      res[[nm]]$prior.value <- prior.value + bump
      res[[nm]]$prior.grad <- prior.grad + bump
    }
    res
  }
}

test_that("the map log-prior is added ONCE, not once per site", {
  counter <- new.env(); counter$n <- 0L; counter$calls <- list()
  fed <- dsNoLimitsClient:::.fed_closure(theta0, s, "prep", "map", 1L, sites,
                                         stub_map_agg(targets, counter))
  z <- c(0.1, -0.2, 0.3)
  theta <- theta0 + s * z
  r <- fed$eval(z)

  lik <- sum(vapply(targets, function(mu) -sum((theta - mu)^2) / 2, numeric(1)))
  lik.grad <- Reduce(`+`, lapply(targets, function(mu) -(theta - mu)))
  expect_equal(r$value, lik + (-2.5))
  expect_equal(r$grad, lik.grad + c(0.1, -0.2, 0.3))
  expect_equal(r$fn, -r$value)
  expect_equal(r$gr, -(s * r$grad))

  # Sanity: two sites, so a double-counted prior would be 2 * -2.5.
  expect_false(isTRUE(all.equal(r$value, lik + 2 * (-2.5))))
})

test_that("no prior block is added for the additive estimators", {
  counter <- new.env(); counter$n <- 0L; counter$calls <- list()
  for (est in c("laplace", "focei", "ghq", "pooled", "mle")) {
    fed <- dsNoLimitsClient:::.fed_closure(theta0, s, "prep", est, 1L, sites,
                                           stub_map_agg(targets, counter))
    r <- fed$eval(c(0.1, -0.2, 0.3))
    theta <- theta0 + s * c(0.1, -0.2, 0.3)
    expect_equal(r$value,
                 sum(vapply(targets, function(mu) -sum((theta - mu)^2) / 2, numeric(1))),
                 info = est)
  }
})

test_that("sites that disagree on the prior block abort, naming both", {
  counter <- new.env(); counter$n <- 0L; counter$calls <- list()
  fed <- dsNoLimitsClient:::.fed_closure(theta0, s, "prep", "map", 1L, sites,
                                         stub_map_agg(targets, counter, bad = "bravo"))
  e <- tryCatch(fed$eval(c(0, 0, 0)), error = function(e) e)
  expect_s3_class(e, "error")
  msg <- paste(conditionMessage(e), collapse = " ")
  expect_match(msg, "alpha")
  expect_match(msg, "bravo")
  expect_match(msg, "not running the same model")
})

test_that("a missing prior block on a map fit aborts with a version hint", {
  counter <- new.env(); counter$n <- 0L; counter$calls <- list()
  fed <- dsNoLimitsClient:::.fed_closure(theta0, s, "prep", "map", 1L, sites,
                                         stub_agg(targets, counter))
  expect_error(fed$eval(c(0, 0, 0)), "log-prior block")
})
