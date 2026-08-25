# Tier 4: the REAL inst/DATASHIELD manifest, as parsed from the installed
# server package. A corrupt DCF file (a blank line, a trailing comma) is
# invisible to every other test in the suite, because the rest of the harness
# registers its methods explicitly.

test_that("the installed dsNoLimits manifest declares exactly the intended methods", {
  skip_if_not_installed("dsNoLimits")
  cfg <- DSLite::defaultDSConfiguration(include = "dsNoLimits")

  expect_setequal(cfg$AggregateMethods$name, DSNOLIMITS_AGGREGATE)
  expect_setequal(cfg$AssignMethods$name, DSNOLIMITS_ASSIGN)
  expect_true(all(grepl("^dsNoLimits::", cfg$AggregateMethods$value)))
  expect_true(all(grepl("^dsNoLimits::", cfg$AssignMethods$value)))
  expect_equal(cfg$AggregateMethods$value,
               paste0("dsNoLimits::", cfg$AggregateMethods$name))
  expect_equal(cfg$AssignMethods$value,
               paste0("dsNoLimits::", cfg$AssignMethods$name))

  # Option values arrive as strings, and every one of ours must be declared.
  expect_true(all(c("default.dsNoLimits.modelDir", "default.dsNoLimits.maxRounds",
                    "default.dsNoLimits.dp.ledgerDir", "default.dsNoLimits.dp.clip",
                    "default.dsNoLimits.dp.clipMode",
                    "default.dsNoLimits.dp.noiseMultiplier",
                    "default.dsNoLimits.dp.epsilonBudget",
                    "default.dsNoLimits.dp.delta", "default.dsNoLimits.dp.maxT") %in%
                    names(cfg$Options)))
  expect_equal(cfg$Options$default.dsNoLimits.maxRounds, "500")
  # A character option must carry its own inner quotes: DSLite and Opal both
  # PARSE the value, so an unquoted string would be read as a symbol.
  expect_equal(cfg$Options$default.dsNoLimits.dp.clipMode, "\"per-group\"")
  expect_equal(cfg$Options$default.dsNoLimits.dp.ledgerDir, "\"\"")
  # Unset by default: DP mode refuses until the data owner provisions a ledger.
  expect_equal(eval(parse(text = cfg$Options$default.dsNoLimits.dp.ledgerDir)), "")
})

test_that("a fit through the real manifest reaches the same methods", {
  skip_if_not_installed("dsNoLimits")
  cfg <- DSLite::defaultDSConfiguration(include = "dsNoLimits")
  srv <- DSLite::newDSLiteServer(tables = list(m1 = nlc_site("a")),
                                 config = cfg, strict = TRUE,
                                 home = tempfile("dslite-"))
  assign("dslite.server", srv, envir = environment())
  options(datashield.env = environment())
  b <- DSI::newDSLoginBuilder(.silent = TRUE)
  b$append(server = "m1", url = "dslite.server", table = "m1", driver = "DSLiteDriver")
  conns <- DSI::datashield.login(b$build(), assign = TRUE, symbol = "D")
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  res <- DSI::datashield.aggregate(conns, call("nolimitsStatusDS"))
  expect_named(res, "m1")
  expect_true("dsnolimits.version" %in% names(res$m1))

  # The server-package version guard (ds.nolimitsFit(checks = TRUE)) resolves a
  # package name only from a config built by defaultDSConfiguration(): DSLite
  # derives the package from the method table, and an explicitly registered
  # method carries no package field. This is therefore the only harness in the
  # suite that can exercise it.
  expect_true(dsNoLimitsClient:::.require_server_pkg(
    conns, "dsNoLimits", dsNoLimitsClient:::.DSNOLIMITS_MIN_VERSION))
  expect_error(
    dsNoLimitsClient:::.require_server_pkg(conns, "dsNoLimits", "99.0.0"),
    "< 99.0.0", fixed = TRUE
  )

  # strict = TRUE is really on: an unpermitted root function is refused.
  expect_error(DSI::datashield.aggregate(conns, quote(mean(D$conc))),
               "There are some DataSHIELD errors", fixed = TRUE)
})
