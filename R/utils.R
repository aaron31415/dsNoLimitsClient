## Internal helpers. Nothing here is exported.

#' @title Resolve the connection list
#' @param datasources a list of [DSI::DSConnection-class] objects, or NULL.
#' @return the resolved list of connections.
#' @keywords internal
#' @noRd
.setDatasources <- function(datasources) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  ok <- is.list(datasources) && length(datasources) > 0L &&
    !is.null(names(datasources)) && all(nzchar(names(datasources))) &&
    all(vapply(datasources, function(d) methods::is(d, "DSConnection"), logical(1)))
  if (!ok) {
    cli::cli_abort(
      "{.arg datasources} must be a named non-empty list of {.cls DSConnection} objects"
    )
  }
  datasources
}

#' @title One aggregate call over all servers, naming any that failed
#' @description Reports server *names*, never server messages; the analyst gets
#'   the text from [DSI::datashield.errors()].
#' @param datasources the connection list.
#' @param cally a `call` object.
#' @param what a short description used in the warning.
#' @return the named list of per-server results, re-indexed by name.
#' @keywords internal
#' @noRd
.agg <- function(datasources, cally, what) {
  failed <- character(0)
  res <- DSI::datashield.aggregate(
    datasources, cally,
    error = function(server, message) failed <<- c(failed, server),
    errors.print = FALSE
  )
  if (length(failed)) {
    warning(what, " failed on: ", paste(sort(failed), collapse = ", "),
            ". Run DSI::datashield.errors() for details.", call. = FALSE)
  }
  # A failed server is DELETED from the result list, never NULL: intersect()
  # keeps us from inventing an <NA> element here.
  res[intersect(setdiff(names(datasources), failed), names(res))]
}

#' @title Abort when a result is missing sites
#' @param res the list returned by `.agg()`.
#' @param sites the expected server names.
#' @param what a short description used in the error.
#' @return `res`, invisibly, or aborts.
#' @keywords internal
#' @noRd
.check_complete <- function(res, sites, what) {
  absent <- setdiff(sites, names(res))
  if (length(absent)) {
    cli::cli_abort(c(
      "x" = "No {what} returned by: {.val {absent}}",
      "i" = "see {.run DSI::datashield.errors()}"
    ))
  }
  invisible(res)
}

#' @title Require a server-side package of at least a given version
#' @description Replaces `DSI::datashield.pkg_check()`, which is broken:
#'   `version_status` is a matrix and `pkg_check` subsets it with `$`.
#' @param datasources the connection list.
#' @param name the server package name.
#' @param version the minimum version, as a character string.
#' @return invisible TRUE, or aborts.
#' @keywords internal
#' @noRd
.require_server_pkg <- function(datasources, name, version) {
  status <- DSI::datashield.pkg_status(datasources)$version_status
  if (!name %in% rownames(status)) {
    stop("Server package '", name, "' is not installed on any server (>= ",
         version, " required).", call. = FALSE)
  }
  for (server in names(datasources)) {
    # Both dimensions at once: a one-column matrix drops the names otherwise.
    v <- status[name, server]
    if (is.na(v)) {
      stop("Server package '", name, "' is missing on server '", server,
           "' (>= ", version, " required).", call. = FALSE)
    }
    if (numeric_version(v) < numeric_version(version)) {
      stop("Server package '", name, "' on server '", server, "' is ",
           v, " < ", version, ".", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' @title Check that every site prepared the same model
#' @description `names` and `theta0` are literal model-derived quantities, so
#'   exact equality is the right test - do not weaken it. `scale` is only a
#'   preconditioning scale (it affects step sizing, never the optimum), and a
#'   heterogeneous-BLAS node can perturb it at the last bits, so it is compared
#'   with a tiny relative tolerance; site 1's scale is used on agreement. A real
#'   disagreement beyond tolerance still aborts, naming both sites.
#' @param infos the named list returned by `nolimitsPrepInfoDS`.
#' @return the first site's info list, or aborts.
#' @keywords internal
#' @noRd
agree <- function(infos) {
  ref.name <- names(infos)[1L]
  ref <- infos[[1L]]
  for (nm in names(infos)[-1L]) {
    for (f in c("names", "theta0")) {
      if (!identical(infos[[nm]][[f]], ref[[f]])) {
        cli::cli_abort(c(
          "x" = "Sites {.val {ref.name}} and {.val {nm}} disagree on {.field {f}}.",
          "i" = "the sites are not running the same model"
        ))
      }
    }
    if (!isTRUE(all.equal(infos[[nm]][["scale"]], ref[["scale"]],
                          tolerance = 1e-8))) {
      cli::cli_abort(c(
        "x" = "Sites {.val {ref.name}} and {.val {nm}} disagree on {.field scale}.",
        "i" = "the sites are not running the same model"
      ))
    }
  }
  ref
}

#' @title Check that every site resolved the byte-identical registered model
#' @description Each site's `model.hash` is the SHA-256 of the canonical source
#'   of the model it actually prepared, so a strict-mode fit is only sound if
#'   every site's hash is the same. This is ADDITIVE to `agree()`: two different
#'   model files can share `names`, `theta0` and `scale` yet differ in source,
#'   and only the hash catches that. Exact equality is the right test - the hash
#'   is entirely model-derived.
#' @param infos the named list returned by `nolimitsPrepInfoDS`.
#' @return the agreed model hash (character(1)), or aborts naming both sites.
#' @keywords internal
#' @noRd
.agree_model <- function(infos) {
  ref.name <- names(infos)[1L]
  ref <- as.character(infos[[1L]][["model.hash"]])
  if (length(ref) != 1L || is.na(ref) || !nzchar(ref)) {
    min.version <- .DSNOLIMITS_MIN_VERSION
    cli::cli_abort(c(
      "x" = "Site {.val {ref.name}} returned no model hash.",
      "i" = "the model-agreement check needs {.pkg dsNoLimits} >= {.val {min.version}}"
    ))
  }
  for (nm in names(infos)[-1L]) {
    h <- as.character(infos[[nm]][["model.hash"]])
    if (length(h) != 1L || is.na(h)) h <- "NA"
    if (!identical(h, ref)) {
      sref <- substr(ref, 1L, 12L)
      sh <- substr(h, 1L, 12L)
      cli::cli_abort(
        "the sites did not agree on the same model (model hash mismatch: {ref.name} {sref} vs {nm} {sh})"
      )
    }
  }
  ref
}

#' @title Check that every site agrees on the MCEM q1/q2 partition
#' @description The partition is model-derived (names of the observation-side and
#'   the random-effect-distribution parameters), so it must be identical across
#'   sites; exact equality is the right test, the same rationale as `agree()`. On
#'   agreement, maps each part's names to indices in the shared parameter order.
#' @param infos the named list returned by `nolimitsPrepInfoDS`.
#' @param par.names the agreed parameter names.
#' @param estimator `"mcem"` or `"saem"`; SAEM additionally agrees and routes the
#'   closed-form / numerical split.
#' @return a list with `q1.names`, `q2.names`, `q1i`, `q2i`, and, for SAEM,
#'   `cf.names`, `num.names` and `num.by.part` (the numerical names of each part,
#'   in part order), or aborts.
#' @keywords internal
#' @noRd
.agree_em <- function(infos, par.names, estimator = "mcem") {
  ref.name <- names(infos)[1L]
  ref <- infos[[1L]]
  min.version <- .DSNOLIMITS_MIN_VERSION
  is.saem <- identical(estimator, "saem")
  needed <- if (is.saem) c("q1.names", "q2.names", "cf.names", "num.names") else
    c("q1.names", "q2.names")
  if (!isTRUE(ref$em.available) || any(vapply(needed, function(f) is.null(ref[[f]]), logical(1)))) {
    cli::cli_abort(c(
      "x" = "Site {.val {ref.name}} does not provide the {estimator} partition.",
      "i" = "{.code estimator = \"{estimator}\"} needs {.pkg dsNoLimits} >= {.val {min.version}} on NoLimits >= 0.2.7"
    ))
  }
  for (nm in names(infos)[-1L]) {
    for (f in needed) {
      if (!identical(infos[[nm]][[f]], ref[[f]])) {
        cli::cli_abort(c(
          "x" = "Sites {.val {ref.name}} and {.val {nm}} disagree on {.field {f}}.",
          "i" = "the sites are not running the same model"
        ))
      }
    }
  }
  q1 <- as.character(ref$q1.names)
  q2 <- as.character(ref$q2.names)
  bad <- setdiff(c(q1, q2), par.names)
  if (length(bad)) {
    cli::cli_abort("the {estimator} partition refers to unknown parameter{?s}: {.val {bad}}")
  }
  out <- list(q1.names = q1, q2.names = q2,
              q1i = match(q1, par.names), q2i = match(q2, par.names))
  if (is.saem) {
    out$cf.names <- as.character(ref$cf.names)
    out$num.names <- as.character(ref$num.names)
    # Each numerical param is optimized with the Q kernel of the part it lives in.
    out$num.by.part <- list(
      q1 = intersect(out$num.names, q1),
      q2 = intersect(out$num.names, q2))
  }
  out
}

#' @title Run the nested federated MCEM loop
#' @description Per outer iteration: (1) a LOCAL, stateful E-step on every site
#'   (each caches its own subjects' fixed draws; the draws never leave the site,
#'   so this is a per-site assign, the one place a per-site argument - the site
#'   id - forces a short loop over connections); (2) a federated L-BFGS-B over the
#'   q1 coordinates and (3) another over the q2 coordinates, each objective eval
#'   summing the sites' per-part Q value and gradient at the fixed draws. Fixed
#'   outer budget, no convergence test; MCEM tolerates approximate M-steps, so
#'   each inner optimisation is capped at `maxfun` iterations.
#' @param p,theta0,s the parameter count, transformed-scale start and scale.
#' @param q1i,q2i integer indices of each part in the parameter vector.
#' @param datasources the connection list.
#' @param sites the server names.
#' @param site.ids a named integer vector, the stable per-site id.
#' @param spec a list with `newobj`, `outer.iters`, `sample.schedule`, `maxfun`,
#'   `seed` and `agg_fun`.
#' @return a list with `theta.t` (the transformed-scale optimum) and `rounds`.
#' @keywords internal
#' @noRd
.em_run <- function(p, theta0, s, q1i, q2i, datasources, sites, site.ids, spec) {
  theta.t <- theta0
  rounds <- 0L

  estep <- function(outer) {
    for (k in seq_along(sites)) {
      DSI::datashield.assign.expr(
        datasources[k], spec$newobj,
        call("nolimitsEmEStepDS", spec$newobj, .ds_num_encode(theta.t),
             as.numeric(outer), as.numeric(spec$seed),
             as.numeric(site.ids[[sites[k]]]), as.numeric(spec$sample.schedule)))
    }
    rounds <<- rounds + 1L
  }

  mstep <- function(part, idxs) {
    if (!length(idxs)) return(invisible())
    ss <- s[idxs]
    x0 <- theta.t[idxs]
    memo_z <- NULL
    memo <- NULL
    ev <- function(z) {
      if (!is.null(memo_z) && identical(z, memo_z)) return(memo)
      cand <- theta.t
      cand[idxs] <- x0 + ss * z
      rounds <<- rounds + 1L
      res <- spec$agg_fun(call("nolimitsEmQDS", spec$newobj,
                               .ds_num_encode(cand), part, as.numeric(rounds)))
      .check_complete(res, sites, "Q-function")
      finite <- all(vapply(res, function(r) isTRUE(r$finite), logical(1)))
      memo <<- if (finite) {
        Q <- sum(vapply(res, function(r) as.numeric(r$value), numeric(1)))
        g <- Reduce(`+`, lapply(res, function(r) as.numeric(r$grad)))
        # Ascent problem sent to a minimiser: negate value and gradient. The
        # gradient is over the part's free axes only.
        list(fn = -Q, gr = -(ss * g))
      } else {
        list(fn = 1e12, gr = rep(0, length(idxs)))
      }
      memo_z <<- z
      memo
    }
    opt <- stats::optim(rep(0, length(idxs)),
                        fn = function(z) ev(z)$fn, gr = function(z) ev(z)$gr,
                        method = "L-BFGS-B",
                        control = list(maxit = spec$maxfun, factr = 1e7, pgtol = 1e-5))
    theta.t[idxs] <<- x0 + ss * opt$par
    invisible()
  }

  for (outer in seq_len(spec$outer.iters)) {
    estep(outer)
    mstep("q1", q1i)
    mstep("q2", q2i)
  }
  list(theta.t = theta.t, rounds = rounds)
}

#' @title Run the nested federated SAEM loop (hybrid M-step)
#' @description Per outer iteration: (1) a LOCAL, stateful E-step on every site
#'   (reused from MCEM); (2) a STATS round - every site returns its DE-NORMALIZED
#'   per-subject-additive sufficient statistics and the client sums the flat
#'   vectors coordinate-wise; (3) a CLOSED-FORM round on the COORDINATOR site
#'   (lowest site name) - a stateful `saem_closed_form_mstep` at this iteration's
#'   gamma, threaded through an assign that updates the coordinator's cache and an
#'   aggregate that reads the update back; (4) a FEDERATED L-BFGS-B over the
#'   NUMERICAL params of each part (the MCEM Q kernel restricted to the
#'   non-closed-form names). Fixed outer budget, no convergence test.
#' @param p,theta0,s the parameter count, transformed-scale start and scale.
#' @param par.names the agreed parameter names.
#' @param em the `.agree_em(estimator = "saem")` result (q1/q2 + cf/num split).
#' @param datasources the connection list; `coord` the coordinator site name.
#' @param sites the server names; `site.ids` the stable per-site id.
#' @param spec a list with `newobj`, `outer.iters`, `sample.schedule`, `maxfun`
#'   and `seed`.
#' @return a list with `theta.t` (the transformed-scale optimum) and `rounds`.
#' @keywords internal
#' @noRd
.saem_run <- function(p, theta0, s, par.names, em, datasources, sites, site.ids,
                      coord, spec) {
  theta.t <- theta0
  rounds <- 0L
  part.idx <- list(q1 = em$q1i, q2 = em$q2i)

  estep <- function(outer) {
    for (k in seq_along(sites)) {
      DSI::datashield.assign.expr(
        datasources[k], spec$newobj,
        call("nolimitsEmEStepDS", spec$newobj, .ds_num_encode(theta.t),
             as.numeric(outer), as.numeric(spec$seed),
             as.numeric(site.ids[[sites[k]]]), as.numeric(spec$sample.schedule)))
    }
    rounds <<- rounds + 1L
  }

  # Numerical M-step over one part restricted to its numerical free names.
  num_mstep <- function(part, free.names) {
    if (!length(free.names)) return(invisible())
    idxs <- match(free.names, par.names)
    ss <- s[idxs]
    x0 <- theta.t[idxs]
    fn.str <- paste(free.names, collapse = ",")
    memo_z <- NULL
    memo <- NULL
    ev <- function(z) {
      if (!is.null(memo_z) && identical(z, memo_z)) return(memo)
      cand <- theta.t
      cand[idxs] <- x0 + ss * z
      rounds <<- rounds + 1L
      res <- .agg(datasources, call("nolimitsEmQDS", spec$newobj,
                                    .ds_num_encode(cand), part, as.numeric(rounds),
                                    fn.str), "Q-function")
      .check_complete(res, sites, "Q-function")
      finite <- all(vapply(res, function(r) isTRUE(r$finite), logical(1)))
      memo <<- if (finite) {
        Q <- sum(vapply(res, function(r) as.numeric(r$value), numeric(1)))
        g <- Reduce(`+`, lapply(res, function(r) as.numeric(r$grad)))
        list(fn = -Q, gr = -(ss * g))
      } else {
        list(fn = 1e12, gr = rep(0, length(idxs)))
      }
      memo_z <<- z
      memo
    }
    opt <- stats::optim(rep(0, length(idxs)),
                        fn = function(z) ev(z)$fn, gr = function(z) ev(z)$gr,
                        method = "L-BFGS-B",
                        control = list(maxit = spec$maxfun, factr = 1e7, pgtol = 1e-5))
    theta.t[idxs] <<- x0 + ss * opt$par
    invisible()
  }

  for (outer in seq_len(spec$outer.iters)) {
    ## (1) E-step: local, per-site (draws never leave a site).
    estep(outer)
    ## (2) Stats: sum the DE-NORMALIZED flat sufficient statistics coordinate-wise.
    stats <- .agg(datasources, call("nolimitsSaemStatsDS", spec$newobj,
                                    .ds_num_encode(theta.t), as.numeric(outer)),
                  "sufficient statistics")
    .check_complete(stats, sites, "sufficient statistics")
    summed <- Reduce(`+`, lapply(stats, function(r) as.numeric(r$stats)))
    rounds <- rounds + 1L
    ## (3) Closed-form M-step on the COORDINATOR: assign updates the cache,
    ## aggregate reads the transformed-scale updates back.
    DSI::datashield.assign.expr(
      datasources[coord], spec$newobj,
      call("nolimitsSaemMstepDS", spec$newobj, .ds_num_encode(theta.t),
           .ds_num_encode(summed), as.numeric(outer),
           as.numeric(spec$outer.iters)))
    rounds <- rounds + 1L
    upd <- .agg(datasources[coord], call("nolimitsSaemUpdatesDS", spec$newobj),
                "closed-form update")
    .check_complete(upd, coord, "closed-form update")
    u <- upd[[coord]]
    unames <- as.character(u$names)
    uvals <- as.numeric(u$values)
    if (length(unames)) theta.t[match(unames, par.names)] <- uvals
    rounds <- rounds + 1L
    ## (4) Numerical M-step over the non-closed-form params of each part.
    for (part in c("q1", "q2")) {
      num_mstep(part, em$num.by.part[[part]])
    }
  }
  list(theta.t = theta.t, rounds = rounds)
}

#' @title One differentially private EM M-step release, summed over sites
#' @description Calls `nolimitsEmQDpDS` for one part (optionally restricted to a
#'   subset of free names) and sums the per-site NOISED, clipped gradients. The
#'   servers clip AND noise in the preconditioned coordinates, so the returned
#'   gradient is already `dQ/dz`; no objective value is released.
#' @param agg the aggregate closure over the connections.
#' @param sites the server names; `newobj` the prep symbol; `enc` the encoded
#'   theta; `part`/`fn.str` the part and (possibly empty) free-name subset;
#'   `round` the fit round; `k` the expected gradient length.
#' @return a list with `grad` (the summed noised gradient) and `per.site`.
#' @keywords internal
#' @noRd
.dp_em_round <- function(agg, sites, newobj, enc, part, fn.str, round, k) {
  res <- agg(call("nolimitsEmQDpDS", newobj, enc, part, as.numeric(round), fn.str))
  .check_complete(res, sites, "differentially private gradient")
  grads <- lapply(res, function(r) as.numeric(r$grad))
  ok <- vapply(grads, function(g) length(g) == k && all(is.finite(g)), logical(1))
  if (!all(ok)) {
    cli::cli_abort(c(
      "x" = "Site{?s} {.val {names(res)[!ok]}} returned an unusable private gradient.",
      "i" = "see {.run DSI::datashield.errors()}"
    ))
  }
  list(grad = Reduce(`+`, grads), per.site = res)
}

#' @title Fixed-schedule DP-Adam over one M-step part
#' @description The private replacement for the L-BFGS-B M-step: `dp.rounds` Adam
#'   steps on the noised part gradient, no line search (a line search re-evaluates
#'   the objective, which is not released under DP and would spend budget). Ascent
#'   on `z`; the step is not multiplied by `s` again because the servers already
#'   clipped and noised in the preconditioned coordinates.
#' @param theta.t the current transformed-scale vector; `idxs` the part's indices;
#'   `ss` the preconditioning scale at those indices; `part`/`fn.str` the part and
#'   free-name subset; `spec` carries `newobj`, `dp.rounds`, `dp.lr`, `dp_agg`.
#' @return a list with the updated `theta.t` and the final round's `last` result.
#' @keywords internal
#' @noRd
.dp_adam_part <- function(theta.t, idxs, ss, part, fn.str, spec, sites) {
  x0 <- theta.t[idxs]
  z <- rep(0, length(idxs)); m <- z; v <- z
  last <- NULL
  for (t in seq_len(spec$dp.rounds)) {
    cand <- theta.t
    cand[idxs] <- x0 + ss * z
    r <- .dp_em_round(spec$dp_agg, sites, spec$newobj, .ds_num_encode(cand),
                      part, fn.str, t, length(idxs))
    g <- r$grad
    m <- 0.9 * m + 0.1 * g
    v <- 0.999 * v + 0.001 * g * g
    z <- z + spec$dp.lr * (m / (1 - 0.9^t)) / (sqrt(v / (1 - 0.999^t)) + 1e-8)
    last <- r
  }
  theta.t[idxs] <- x0 + ss * z
  list(theta.t = theta.t, last = last)
}

#' @title Run the differentially private federated MCEM loop
#' @description As `.em_run()`, but each M-step part is optimised by DP-Adam over
#'   the servers' per-subject-clipped, noised part gradients (`nolimitsEmQDpDS`).
#'   The E-step is unchanged and releases nothing. The privacy composition is
#'   `outer.iters * (number of non-empty parts) * dp.rounds` releases.
#' @return a list with `theta.t`, `rounds`, `releases` and `last` (the final DP
#'   round's per-site results, for the accountant report).
#' @keywords internal
#' @noRd
.em_run_dp <- function(p, theta0, s, q1i, q2i, datasources, sites, site.ids, spec) {
  theta.t <- theta0
  releases <- 0L
  last <- NULL

  estep <- function(outer) {
    for (k in seq_along(sites)) {
      DSI::datashield.assign.expr(
        datasources[k], spec$newobj,
        call("nolimitsEmEStepDS", spec$newobj, .ds_num_encode(theta.t),
             as.numeric(outer), as.numeric(spec$seed),
             as.numeric(site.ids[[sites[k]]]), as.numeric(spec$sample.schedule)))
    }
  }

  parts <- Filter(function(pr) length(pr$idxs) > 0L,
                  list(list(part = "q1", idxs = q1i),
                       list(part = "q2", idxs = q2i)))
  for (outer in seq_len(spec$outer.iters)) {
    estep(outer)
    for (pr in parts) {
      out <- .dp_adam_part(theta.t, pr$idxs, s[pr$idxs], pr$part, "", spec, sites)
      theta.t <- out$theta.t
      last <- out$last
      releases <- releases + spec$dp.rounds
    }
  }
  expected <- spec$outer.iters * length(parts) * spec$dp.rounds
  if (releases != expected) {
    cli::cli_abort("internal: DP MCEM release count {releases} != composition {expected}")
  }
  list(theta.t = theta.t, rounds = releases, releases = releases, last = last)
}

#' @title Run the differentially private federated SAEM loop (hybrid M-step)
#' @description As `.saem_run()`, but (i) the stats round releases per-subject-
#'   clipped, noised de-normalized statistics (`nolimitsSaemStatsDpDS`) and (ii)
#'   the numerical M-step of each part is DP-Adam over the noised part gradient.
#'   The coordinator's closed-form M-step runs on the NOISED summed stats and is
#'   pure post-processing: it releases nothing and charges no budget. The privacy
#'   composition is `outer.iters * (1 stats + numerical_parts * dp.rounds)`.
#' @return a list with `theta.t`, `rounds`, `releases` and `last`.
#' @keywords internal
#' @noRd
.saem_run_dp <- function(p, theta0, s, par.names, em, datasources, sites, site.ids,
                         coord, spec) {
  theta.t <- theta0
  releases <- 0L
  last <- NULL

  estep <- function(outer) {
    for (k in seq_along(sites)) {
      DSI::datashield.assign.expr(
        datasources[k], spec$newobj,
        call("nolimitsEmEStepDS", spec$newobj, .ds_num_encode(theta.t),
             as.numeric(outer), as.numeric(spec$seed),
             as.numeric(site.ids[[sites[k]]]), as.numeric(spec$sample.schedule)))
    }
  }

  num.parts <- Filter(function(part) length(em$num.by.part[[part]]) > 0L,
                      c("q1", "q2"))
  for (outer in seq_len(spec$outer.iters)) {
    ## (1) E-step: local, per-site (draws never leave a site).
    estep(outer)
    ## (2) DP stats round: sum the NOISED de-normalized stats coordinate-wise.
    res <- .agg(datasources, call("nolimitsSaemStatsDpDS", spec$newobj,
                                  .ds_num_encode(theta.t), as.numeric(outer)),
                "differentially private sufficient statistics")
    .check_complete(res, sites, "differentially private sufficient statistics")
    summed <- Reduce(`+`, lapply(res, function(r) as.numeric(r$stats)))
    releases <- releases + 1L
    last <- list(per.site = res)
    ## (3) Closed-form M-step on the COORDINATOR: post-processing on noised stats,
    ## no release, no budget charge.
    DSI::datashield.assign.expr(
      datasources[coord], spec$newobj,
      call("nolimitsSaemMstepDS", spec$newobj, .ds_num_encode(theta.t),
           .ds_num_encode(summed), as.numeric(outer),
           as.numeric(spec$outer.iters)))
    upd <- .agg(datasources[coord], call("nolimitsSaemUpdatesDS", spec$newobj),
                "closed-form update")
    .check_complete(upd, coord, "closed-form update")
    u <- upd[[coord]]
    unames <- as.character(u$names)
    uvals <- as.numeric(u$values)
    if (length(unames)) theta.t[match(unames, par.names)] <- uvals
    ## (4) DP numerical M-step over the non-closed-form params of each part.
    for (part in num.parts) {
      free.names <- em$num.by.part[[part]]
      idxs <- match(free.names, par.names)
      out <- .dp_adam_part(theta.t, idxs, s[idxs], part,
                           paste(free.names, collapse = ","), spec, sites)
      theta.t <- out$theta.t
      last <- out$last
      releases <- releases + spec$dp.rounds
    }
  }
  expected <- spec$outer.iters * (1L + length(num.parts) * spec$dp.rounds)
  if (releases != expected) {
    cli::cli_abort("internal: DP SAEM release count {releases} != composition {expected}")
  }
  list(theta.t = theta.t, rounds = releases, releases = releases, last = last)
}

## Wire codecs. These are LOCAL COPIES of the identically-named helpers in the
## server package dsNoLimits; test-arg-codec.R pins them to the same reference
## vector so the pair cannot drift. Copies, not an Imports: the two packages are
## installed on different machines.

#' @title Encode a double vector exactly
#' @param x a numeric vector.
#' @return a single comma-joined hexadecimal-float string.
#' @keywords internal
#' @noRd
.ds_num_encode <- function(x) paste0(sprintf("%a", x), collapse = ",")

#' @title Decode a double vector encoded by `.ds_num_encode()`
#' @param s the string.
#' @return a numeric vector.
#' @keywords internal
#' @noRd
.ds_num_decode <- function(s) as.numeric(unlist(strsplit(s, ",", fixed = TRUE)))

#' @title base64url-encode arbitrary text for the restricted parser
#' @param x a single character string.
#' @return the encoded string, over `[A-Za-z0-9_-]`.
#' @keywords internal
#' @noRd
.ds_encode <- function(x) {
  # Pinned identical to the server copy: an empty string must not encode to "",
  # which is an Opal lexer error, so it becomes NULL just as it does server-side.
  if (!nzchar(x)) return(NULL)
  # base64_enc line-wraps anything over 64 bytes; the newlines must go before
  # the "=" strip, or the payload picks up characters no grammar accepts.
  b64 <- gsub("[[:space:]]", "", jsonlite::base64_enc(charToRaw(enc2utf8(x))))
  gsub("=", "", chartr("+/", "-_", b64), fixed = TRUE)
}

#' @title Decode a string encoded by `.ds_encode()`
#' @param x the encoded string.
#' @return the original text.
#' @keywords internal
#' @noRd
.ds_decode <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]", "", x)
  pad <- strrep("=", (4L - nchar(x) %% 4L) %% 4L)
  rawToChar(jsonlite::base64_dec(chartr("-_", "+/", paste0(x, pad))))
}

## Client-side argument validation. UX only, never security: the server
## revalidates everything.

#' @title Validate a single syntactic R name
#' @param x the value supplied.
#' @param arg the argument name.
#' @return invisible TRUE, or aborts.
#' @keywords internal
#' @noRd
.check_name <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !grepl("^[A-Za-z][A-Za-z0-9._]*$", x)) {
    cli::cli_abort(
      "{.arg {arg}} must be a single syntactic R name given as a character string"
    )
  }
  invisible(TRUE)
}

#' @title Validate a single whole number in a range
#' @param x the value supplied.
#' @param arg the argument name.
#' @param lo,hi the inclusive bounds.
#' @return the value as an integer, or aborts.
#' @keywords internal
#' @noRd
.check_count <- function(x, arg, lo, hi) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x != trunc(x) ||
      x < lo || x > hi) {
    cli::cli_abort("{.arg {arg}} must be a single whole number in {lo}:{hi}")
  }
  as.integer(x)
}

#' @title Validate a single TRUE/FALSE
#' @param x the value supplied.
#' @param arg the argument name.
#' @return the value, or aborts.
#' @keywords internal
#' @noRd
.check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    cli::cli_abort("{.arg {arg}} must be {.code TRUE} or {.code FALSE}")
  }
  x
}

#' @title Check that every site returned the same MAP log-prior block
#' @description The log-prior is data-free: it is a function of the model's
#'   declared priors and its parameter transform, evaluated at the parameter
#'   vector this client just sent. So it must be bit-identical across sites, and
#'   exact equality is the right test - the same rationale as `agree()`.
#' @param res the per-site result list.
#' @param p the parameter count.
#' @return a list with `value` and `grad`, or aborts.
#' @keywords internal
#' @noRd
.prior_block <- function(res, p) {
  ref.name <- names(res)[1L]
  ref <- res[[ref.name]]
  min.version <- .DSNOLIMITS_MIN_VERSION
  if (!is.numeric(ref[["prior.value"]]) || length(ref[["prior.value"]]) != 1L ||
      !is.numeric(ref[["prior.grad"]]) || length(ref[["prior.grad"]]) != p) {
    cli::cli_abort(c(
      "x" = "Site {.val {ref.name}} returned no usable MAP log-prior block.",
      "i" = "{.code estimator = \"map\"} needs {.pkg dsNoLimits} >= {.val {min.version}}"
    ))
  }
  for (nm in names(res)[-1L]) {
    for (f in c("prior.value", "prior.grad")) {
      if (!identical(res[[nm]][[f]], ref[[f]])) {
        cli::cli_abort(c(
          "x" = "Sites {.val {ref.name}} and {.val {nm}} disagree on {.field {f}}.",
          "i" = "the sites are not running the same model"
        ))
      }
    }
  }
  list(value = as.numeric(ref[["prior.value"]]),
       grad = as.numeric(ref[["prior.grad"]]))
}

#' @title One federated objective/gradient round
#' @description For every estimator but `"map"` the objective and its gradient
#'   are plain unweighted sums over sites. For `"map"` the sites return the
#'   likelihood-only part and the log-prior separately, and the prior is added
#'   ONCE: summing it with the rest would count it once per site.
#' @param theta the transformed-scale parameter vector.
#' @param round the round number, sent to the server as metadata.
#' @param spec a list carrying `newobj`, `estimator`, `ghq.level`, `p`, `sites`
#'   and `aggregate_fun`.
#' @return a list with `finite`, `value`, `grad`, `natural` and `per.site`.
#' @keywords internal
#' @noRd
.fed_round <- function(theta, round, spec) {
  cally <- call("nolimitsObjGradDS", spec$newobj, .ds_num_encode(theta),
                spec$estimator, as.numeric(spec$ghq.level), as.numeric(round))
  res <- spec$aggregate_fun(cally)
  .check_complete(res, spec$sites, "objective and gradient")
  finite <- all(vapply(res, function(r) isTRUE(r$finite), logical(1)))
  if (!finite) {
    return(list(finite = FALSE, value = NA_real_, grad = rep(0, spec$p),
                natural = rep(NA_real_, spec$p), per.site = res))
  }
  prior <- if (identical(spec$estimator, "map")) {
    .prior_block(res, spec$p)
  } else {
    list(value = 0, grad = rep(0, spec$p))
  }
  list(
    finite = TRUE,
    value = sum(vapply(res, function(r) as.numeric(r$value), numeric(1))) +
      prior$value,
    grad = Reduce(`+`, lapply(res, function(r) as.numeric(r$grad))) + prior$grad,
    natural = as.numeric(res[[1L]]$natural),
    per.site = res
  )
}

#' @title Build the memoised federated objective on the preconditioned scale
#' @description `optim()` calls `fn` and `gr` separately at the same point; the
#'   `identical(z, last_z)` memo keeps that to ONE federated round. Factored out
#'   of `ds.nolimitsFit()` so it is unit-testable without DSI.
#' @param theta0 the transformed-scale starting values.
#' @param s the preconditioning scale.
#' @param newobj the server-side prep symbol.
#' @param estimator one of "laplace", "focei", "ghq", "pooled", "mle", "map".
#' @param ghq.level the quadrature level.
#' @param sites the expected server names.
#' @param aggregate_fun a function of one `call` object returning the named list
#'   of per-site results.
#' @return a list of functions: `eval(z)`, `probe(theta)` and `rounds()`.
#' @keywords internal
#' @noRd
.fed_closure <- function(theta0, s, newobj, estimator, ghq.level, sites,
                         aggregate_fun) {
  p <- length(theta0)
  spec <- list(newobj = newobj, estimator = estimator, ghq.level = ghq.level,
               p = p, sites = sites, aggregate_fun = aggregate_fun)
  rounds <- 0L
  last_z <- NULL
  last <- NULL

  probe <- function(theta) {
    rounds <<- rounds + 1L
    .fed_round(theta, rounds, spec)
  }

  list(
    eval = function(z) {
      if (!is.null(last_z) && identical(z, last_z)) return(last)
      r <- probe(theta0 + s * z)
      # A non-finite objective is a legitimate estimator result at a probe
      # point, not a site failure: make L-BFGS-B backtrack instead of dying.
      last <<- if (r$finite) {
        c(r, list(fn = -r$value, gr = -(s * r$grad)))
      } else {
        c(r, list(fn = 1e12, gr = rep(0, p)))
      }
      last_z <<- z
      last
    },
    probe = probe,
    rounds = function() rounds
  )
}

#' @title sqrt() that yields NA rather than NaN on a negative variance
#' @param v a numeric vector.
#' @return a numeric vector.
#' @keywords internal
#' @noRd
.safe_sqrt <- function(v) ifelse(is.finite(v) & v >= 0, sqrt(v), NA_real_)

#' @title Scalar field of a server result, or NA
#' @param x a list.
#' @param field the field name.
#' @param empty the value to use when the field is absent.
#' @return a length-one vector.
#' @keywords internal
#' @noRd
.field1 <- function(x, field, empty) {
  v <- x[[field]]
  if (is.null(v) || length(v) != 1L) empty else v
}

## Differential privacy. The client neither sets nor can influence any parameter
## of the mechanism: clip, noise multiplier, budget, delta and the round cap are
## all server options set by the data owner. The client only asks for rounds,
## sums what comes back, and reports what the servers said the fit cost.

# Renyi-DP order grid. This is a LOCAL COPY of the accountant in dsNoLimits;
# test-arg-dp.R pins the two copies to the same reference values.
.DP_ALPHAS <- sort(unique(c(seq(1.01, 10.0, length.out = 900), 10:512)))

#' @title Cumulative (epsilon, delta) for a list of Gaussian releases
#' @description One release of the Gaussian mechanism at sensitivity C and noise
#'   standard deviation `sigma * C` is `(alpha, alpha / (2 sigma^2))`-RDP for
#'   every `alpha > 1`, and RDP composes by addition, so a heterogeneous history
#'   costs `alpha * sum_i 1 / (2 sigma_i^2)`. Converting with the standard tail
#'   bound and minimising over alpha gives epsilon.
#' @param sigmas numeric vector of per-release noise multipliers.
#' @param delta the target delta.
#' @return the epsilon.
#' @keywords internal
#' @noRd
.dp_epsilon <- function(sigmas, delta) {
  sigmas <- as.numeric(sigmas)
  if (!length(sigmas)) return(0)
  if (any(!is.finite(sigmas)) || any(sigmas <= 0) ||
      length(delta) != 1L || !is.finite(delta) || delta <= 0 || delta >= 1) {
    cli::cli_abort("invalid arguments to the privacy accountant")
  }
  k <- sum(1 / (2 * sigmas^2))
  min(.DP_ALPHAS * k + log(1 / delta) / (.DP_ALPHAS - 1))
}

#' @title One differentially private federated round
#' @description The servers return only a clipped, noised gradient sum - no
#'   objective value - so there is nothing to check for finiteness beyond the
#'   gradient itself, and nothing to backtrack on. The site sums are added
#'   unweighted, exactly as in the non-private path.
#' @param theta the transformed-scale parameter vector.
#' @param round the round number.
#' @param spec the same list `.fed_round()` takes.
#' @return a list with `grad` and `per.site`.
#' @keywords internal
#' @noRd
.dp_round <- function(theta, round, spec) {
  cally <- call("nolimitsObjGradDpDS", spec$newobj, .ds_num_encode(theta),
                spec$estimator, as.numeric(spec$ghq.level), as.numeric(round))
  res <- spec$aggregate_fun(cally)
  .check_complete(res, spec$sites, "differentially private gradient")
  grads <- lapply(res, function(r) as.numeric(r$grad))
  ok <- vapply(grads, function(g) length(g) == spec$p && all(is.finite(g)),
               logical(1))
  if (!all(ok)) {
    cli::cli_abort(c(
      "x" = "Site{?s} {.val {names(res)[!ok]}} returned an unusable private gradient.",
      "i" = "see {.run DSI::datashield.errors()}"
    ))
  }
  list(grad = Reduce(`+`, grads), per.site = res)
}

#' @title Fixed-schedule Adam ascent on differentially private gradients
#' @description Not L-BFGS-B: a line search re-evaluates the objective to test a
#'   step, which on a noisy gradient is meaningless and spends privacy budget,
#'   and no objective value is released under DP at all. Adam spends exactly
#'   `rounds` releases, which is what makes the round count the budget knob.
#'   There is no convergence test, because every quantity that could gate one is
#'   noisy.
#'
#'   The server clips - and therefore noises - in the PRECONDITIONED
#'   coordinates, so the returned gradient is already `dL/dz`. The step below
#'   must not multiply by `s` a second time.
#' @param p the parameter count.
#' @param theta0 the transformed-scale starting values.
#' @param s the preconditioning scale.
#' @param lr the Adam learning rate.
#' @param rounds how many releases to spend.
#' @param spec the round specification.
#' @return a list with `z`, `rounds` and `last` (the final round's per-site
#'   results).
#' @keywords internal
#' @noRd
.dp_optimize <- function(p, theta0, s, lr, rounds, spec) {
  z <- rep(0, p)
  m <- rep(0, p)
  v <- rep(0, p)
  last <- NULL
  for (t in seq_len(rounds)) {
    r <- .dp_round(theta0 + s * z, t, spec)
    g <- r$grad
    m <- 0.9 * m + 0.1 * g
    v <- 0.999 * v + 0.001 * g * g
    # Ascent: the servers report the gradient of a LOG-likelihood.
    z <- z + lr * (m / (1 - 0.9^t)) / (sqrt(v / (1 - 0.999^t)) + 1e-8)
    last <- r
  }
  list(z = z, rounds = rounds, last = last)
}

#' @title The privacy accounting the servers reported for this fit
#' @param per.site the final DP round's per-site result list.
#' @return a data frame with one row per site, plus the fit-level attributes.
#' @keywords internal
#' @noRd
.dp_report <- function(per.site) {
  num <- function(x, f) as.numeric(.field1(x, f, NA_real_))
  data.frame(
    site = names(per.site),
    releases = vapply(per.site, function(x) as.integer(.field1(x, "releases", NA_integer_)), integer(1)),
    epsilon = vapply(per.site, num, numeric(1), "epsilon"),
    delta = vapply(per.site, num, numeric(1), "delta"),
    remaining.budget = vapply(per.site, num, numeric(1), "remaining.budget"),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
