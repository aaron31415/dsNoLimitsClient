# dsNoLimitsClient

Client-side DataSHIELD package for federated nonlinear mixed-effects (NLME)
analysis of longitudinal data. It fits a NLME model to data held on several
DataSHIELD nodes by driving a federated optimiser over the per-site marginal
log-likelihoods and gradients served by the companion package `dsNoLimits`,
giving estimates and Wald standard errors identical to a pooled analysis without
any individual-level data leaving the servers.

This is the package an analyst installs. Each node must have the server package
[`dsNoLimits`](https://github.com/manuhuth/dsNoLimits) installed and its methods
permitted.

## Usage

The runnable end-to-end examples are under `inst/examples/`.
[`warfarin_dslite.R`](inst/examples/warfarin_dslite.R) fits a one-compartment
oral PK model to the warfarin data split across three in-process DSLite sites
and cross-checks against the pooled fit;
[`warfarin_em_dslite.R`](inst/examples/warfarin_em_dslite.R) does the same for
the SAEM and MCEM estimators. Run one with:

```r
Rscript inst/examples/warfarin_dslite.R
```

The core of a fit, once you have logged into your connections, is a single call:

```r
fit <- ds.nolimitsFit(
  model      = model,        # a registered name, or a model source string
  data       = "D",          # the assigned data frame symbol on each node
  id.col     = "ID",
  time.col   = "t",
  estimator  = "laplace",
  model.mode = "string"
)
print(fit)
summary(fit)
coef(fit)   # estimates; vcov(fit) and confint(fit) also work
```

The fit object carries the estimates, the log-likelihood, Wald standard errors,
per-site contributions, round count and versions.

## Estimators and options

`estimator` is one of `laplace`, `focei`, `ghq`, `pooled`, `mle`, `map`, `mcem`
or `saem`. The marginal-likelihood estimators use random effects; `mle` and
`map` require a model with no random effects; `mcem` and `saem` are
sampling-based EM estimators and report no standard errors, since their
objective is stochastic. Wald standard errors are on by default (`se = TRUE`)
where they apply. Passing `dp = TRUE` runs the fit against the servers'
differentially private gradients with a fixed-schedule optimiser and reports the
spent privacy budget; the privacy parameters themselves are set by the data
owner, not the analyst.

## Federated equals pooled

For the deterministic estimators, a federated fit reproduces a pooled NoLimits
fit to better than 1e-9 relative on the objective; the warfarin ODE end-to-end
example agrees to 1.7e-9. This holds because the marginal log-likelihood is a
sum of independent per-subject terms, so summing site contributions is exact
rather than approximate. MCEM and SAEM agree with their pooled counterparts up
to Monte-Carlo error.

`dsNoLimitsClient` and `dsNoLimits` are versioned in lockstep; install matching
versions on the client and the server.

## Status

Development version, not yet released. It is currently tested against pooled
reference fits under DSLite (in-process); testing and deployment on a real
DataSHIELD server will follow. See [NEWS.md](NEWS.md).

## License

MIT. See [LICENSE](LICENSE).
