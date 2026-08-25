# The canonical DSLite harness (testing.md, "The harness, complete"): empty
# method tables plus explicit registration, so the config is identical under
# devtools::test() and R CMD check and no other installed DS package can shadow
# a method name. strict = TRUE still validates because the tables are non-NULL;
# `config = list()` would skip the parse-and-check step entirely.
#
# DSLite runs every "site" in ONE R process, and therefore over ONE Julia
# session. That is safe here: the fixture model source is byte-identical on
# every site so its re-evaluation into Main.NoLimitsSession is idempotent, and
# each site's prep cache is a separate object bound to the analyst's symbol in
# that site's own DSLite session environment.
ds_lite_conns <- function(tables,
                          aggregate = character(0),
                          assign    = character(0),
                          options   = list(),
                          symbol    = "D",
                          env       = parent.frame()) {
  srv <- DSLite::newDSLiteServer(
    tables = tables,
    config = list(AggregateMethods = data.frame(),
                  AssignMethods    = data.frame(),
                  Options          = list()),
    strict = TRUE,
    home   = tempfile("dslite-")
  )
  for (nm in names(aggregate)) srv$aggregateMethod(nm, aggregate[[nm]])
  for (nm in names(assign))    srv$assignMethod(nm, assign[[nm]])
  for (nm in names(options))   srv$option(nm, options[[nm]])

  assign("dslite.server", srv, envir = env)
  options(datashield.env = env)

  b <- DSI::newDSLoginBuilder(.silent = TRUE)
  for (nm in names(tables)) {
    b$append(server = nm, url = "dslite.server", table = nm, driver = "DSLiteDriver")
  }
  DSI::datashield.login(b$build(), assign = TRUE, symbol = symbol)
}

# DSLite applies server options by pasting `name=value` into options() and
# parsing it, so a character value needs its own inner quotes.
.q <- function(x) paste0('"', x, '"')

# The data owner's DP settings, shared by every DP test file (each test file has
# its own environment, so this must live in a helper to be visible to all).
dp_options <- function(dir, ..., sigma = 1, budget = 1e6, delta = 1e-5) {
  ds_test_options(
    dsNoLimits.dp.ledgerDir = .q(dir),
    default.dsNoLimits.dp.ledgerDir = .q(dir),
    dsNoLimits.dp.clip = 20,
    default.dsNoLimits.dp.clip = 20,
    dsNoLimits.dp.clipMode = .q("per-group"),
    default.dsNoLimits.dp.clipMode = .q("per-group"),
    dsNoLimits.dp.noiseMultiplier = sigma,
    default.dsNoLimits.dp.noiseMultiplier = sigma,
    dsNoLimits.dp.epsilonBudget = budget,
    default.dsNoLimits.dp.epsilonBudget = budget,
    dsNoLimits.dp.delta = delta,
    default.dsNoLimits.dp.delta = delta,
    dsNoLimits.dp.maxT = 200,
    default.dsNoLimits.dp.maxT = 200,
    ...
  )
}

# Disclosure settings for the harness, both spellings of every threshold: a
# server that reads only the live name and a server that reads only the
# `default.`-prefixed one must both see the same value.
#
# nfilter.glm is raised to 0.8 for the CSV fixtures. The saturation gate is
# `p <= nfilter.glm * n_subjects`; with 6 subjects per site and p = 4 the
# canonical 0.33 allows 1.98 parameters and fails. The seeded 30-subject
# generator below is large enough for the 0.33 default, and one equivalence
# test uses it at that default deliberately.
ds_test_options <- function(control = "permissive", ..., nfilter.glm = 0.8) {
  base <- list(
    nfilter.tab = 3, nfilter.subset = 3, nfilter.glm = nfilter.glm,
    nfilter.string = 80, nfilter.stringShort = 20, nfilter.kNN = 3,
    nfilter.levels.density = 0.33, nfilter.levels.max = 40,
    nfilter.noise = 0.25, nfilter.privacy.old = 5
  )
  opts <- c(base, stats::setNames(base, paste0("default.", names(base))))
  opts$datashield.privacyLevel <- 5
  opts$datashield.privacyControlLevel <- .q(control)
  opts$default.datashield.privacyControlLevel <- .q(control)
  c(opts, list(...))
}

# A registry directory holding the fixture model under a registered name, plus
# the two option spellings that point the server at it.
ds_model_dir_options <- function(opts, env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  writeLines(NLDS_TEST_MODEL, file.path(dir, "onecomp_iv.jl"))
  c(opts, list(dsNoLimits.modelDir = .q(dir),
               default.dsNoLimits.modelDir = .q(dir)))
}
