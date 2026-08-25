library(DSI)
library(DSLite)

# The tiers that touch Julia skip themselves when it is not provisioned; see
# nlc_skip_no_julia() in helper-nolimits.R. The MCEM/SAEM tiers need the EM
# dev-API primitives, which ship in NoLimits >= 0.2.7; NoLimitsR resolves them
# from its default shared env.
