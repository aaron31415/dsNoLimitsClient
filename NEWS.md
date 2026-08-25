# dsNoLimitsClient 0.0.5.9000

- `ds.nolimitsFit()` gains `estimator = "mcem"` and `"saem"`: the sampling-based
  EM estimators, driven as a nested loop (a local stateful E-step per site each
  outer iteration, then a federated M-step). SAEM routes its closed-form update
  to a coordinator site. New `em.*` knobs control the outer budget, sample
  schedule and M-step cap; the fit reproduces a pooled `fit_model(MCEM())` /
  `fit_model(SAEM())` up to Monte-Carlo noise. No standard errors and `loglik` /
  `converged` are `NA` (stochastic, fixed budget).
- `dp = TRUE` extends to the EM estimators: the M-steps become fixed-schedule
  DP-Adam over per-subject clipped, noised gradients (MCEM) and sufficient
  statistics (SAEM), the most privacy-expensive path. Release composition is
  `outer * 2 * dp.rounds` (MCEM) and `outer * (1 + numerical_parts * dp.rounds)`
  (SAEM); epsilon is reported from the accountant.

# dsNoLimitsClient 0.0.4.9000

- Cross-site model-agreement check: a strict-mode fit verifies every site
  resolved the byte-identical registered model, aborting on a model-hash
  mismatch and naming the disagreeing sites. The agreed hash is recorded on the
  fit (`fit$model.hash`) and printed in short form. `ds.nolimitsModels` now
  surfaces each model's canonical source hash per site.

# dsNoLimitsClient 0.0.3.9000

Initial development version. Not yet released.

- `ds.nolimitsFit`: federated NLME fit driving a client-side optimiser over the
  per-site marginal log-likelihood and gradient served by `dsNoLimits`.
- Estimators: `laplace`, `focei`, `ghq`, `pooled`, `mle` and `map`.
- Wald standard errors via finite differences of the federated gradient
  (`se = TRUE`), on the transformed and, where defined, natural scale.
- Two model profiles: a registered server-side model name, or an
  analyst-supplied model source string (`model.mode = "string"`).
- Optional differentially private fit (`dp = TRUE`): fixed-schedule optimiser
  against the servers' clipped, noised gradient, reporting the spent privacy
  budget. Privacy parameters are owner-set.
- S3 methods `print`, `summary`, `coef`, `vcov` and `confint` on the fit object.
- Runnable DSLite example in `inst/examples/warfarin_dslite.R`.

Requires the `dsNoLimits` server package installed and permitted on each node.
