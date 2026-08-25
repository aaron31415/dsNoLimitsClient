# Federated one-compartment oral PK fit of the warfarin data over three DSLite
# sites, cross-checked against the same model fitted on the pooled data.
# Run from the dsNoLimitsClient repository root:
#   Rscript inst/examples/warfarin_dslite.R

# Install the three packages from their sibling source directories (with their
# CRAN dependencies) if they are not already available. NoLimitsR additionally
# needs its Julia environment set up once: NoLimitsR::install_nolimits().
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
srcs <- c(NoLimitsR = "../NoLimitsR", dsNoLimits = "../dsNoLimits",
          dsNoLimitsClient = ".")
for (nm in names(srcs)) {
  if (!requireNamespace(nm, quietly = TRUE)) {
    remotes::install_local(srcs[[nm]], dependencies = TRUE, upgrade = "never")
  }
}

library(DSI)
library(DSLite)
library(dsNoLimitsClient)

model <- "
@fixedEffects begin
    ka       = RealNumber(1.0)
    cl       = RealNumber(0.13)
    v        = RealNumber(8.0)
    omega_ka = RealNumber(0.4, scale=:log)
    omega_cl = RealNumber(0.3, scale=:log)
    omega_v  = RealNumber(0.2, scale=:log)
    sigma    = RealNumber(0.5, scale=:log)
end

@covariates begin
    t    = Covariate()
    Dose = ConstantCovariate(constant_on=:ID)
end

@randomEffects begin
    eta_ka = RandomEffect(LogNormal(0.0, omega_ka); column=:ID)
    eta_cl = RandomEffect(LogNormal(0.0, omega_cl); column=:ID)
    eta_v  = RandomEffect(LogNormal(0.0, omega_v);  column=:ID)
end

@preDifferentialEquation begin
    kai = ka * eta_ka
    cli = cl * eta_cl
    vi  = v * eta_v
end

@DifferentialEquation begin
    D(depot)   ~ -kai * depot
    D(central) ~ kai * depot - (cli / vi) * central
end

@initialDE begin
    depot   = Dose
    central = 0.0
end

@formulas begin
    cp = central(t) / vi
    conc ~ Normal(cp, sigma)
end
"

w <- nlmixr2data::warfarin
dose <- stats::setNames(w$amt[w$evid == 1], as.character(w$id[w$evid == 1]))
obs <- w[w$dvid == "cp" & w$evid == 0, ]
d <- data.frame(ID = as.character(obs$id), t = obs$time,
                Dose = unname(dose[as.character(obs$id)]), conc = obs$dv)
cat(sprintf("warfarin: %d subjects, %d observations\n",
            length(unique(d$ID)), nrow(d)))

# Round-robin over SUBJECTS, so no subject is split across sites.
ids <- unique(d$ID)
tables <- lapply(0:2, function(k) d[d$ID %in% ids[seq_along(ids) %% 3 == k], ])
names(tables) <- paste0("site", 1:3)

dslite.server <- DSLite::newDSLiteServer(
  tables = tables,
  config = DSLite::defaultDSConfiguration(include = c("dsBase", "dsNoLimits")),
  strict = TRUE,
  home = tempfile("dslite-")
)
# p = 7 against ~10 subjects per site: the saturation gate p <= nfilter.glm *
# n_subjects needs more than the 0.33 default. In production the data owner
# sets this, not the analyst.
dslite.server$option("nfilter.glm", 0.8)
# The model is sent as source, which is a permissive-profile method.
dslite.server$option("datashield.privacyControlLevel", '"permissive"')

b <- DSI::newDSLoginBuilder(.silent = TRUE)
for (nm in names(tables)) {
  b$append(server = nm, url = "dslite.server", table = nm, driver = "DSLiteDriver")
}
conns <- DSI::datashield.login(b$build(), assign = TRUE, symbol = "D")

fit <- ds.nolimitsFit(model = model, data = "D", id.col = "ID", time.col = "t",
                      estimator = "laplace", model.mode = "string")
print(fit)
print(summary(fit))
DSI::datashield.logout(conns)

nl <- NoLimitsR::nolimits()
dm <- nl$DataModel(NoLimitsR::nl_model(model), NoLimitsR::nl_data(d),
                   primary_id = NoLimitsR::jl_sym("ID"),
                   time_col = NoLimitsR::jl_sym("t"))
pooled <- nl$fit_model(dm, nl$Laplace())
as_vec <- NoLimitsR::nl_eval("as_vec(x) = Vector{Float64}(x)")
theta <- as.numeric(as_vec(nl$get_params(pooled, scale = NoLimitsR::jl_sym("untransformed"))))
# fit_model minimises the negative marginal log-likelihood.
loglik <- -as.numeric(nl$get_objective(pooled))

cat(sprintf("\nloglik federated %.8f  pooled %.8f  rel.diff %.3e\n",
            fit$loglik, loglik, abs(fit$loglik - loglik) / abs(loglik)))
print(data.frame(federated = coef(fit), pooled = theta,
                 abs.diff = abs(coef(fit) - theta)))
