## Every call shape this client can emit, checked against both grammars before
## it can ever reach a server. DSLite's own validator only inspects the root
## node, so ds_wire_check() does the whitelist and charset work itself.

skip_if_not_installed("DSLite")
skip_if_not_installed("rly")

payload12 <- dsNoLimitsClient:::.ds_num_encode(
  c(0, -0, pi, -1 / 3, 1e-300, 1e300, .Machine$double.xmax,
    .Machine$double.eps, -2.5, 7, 1234567.891011, -1e-17)
)

model.src <- paste(
  "@Model begin",
  "  @fixedEffects begin",
  "    v = 8.0; cli = 0.13; sigma = 0.5",
  "  end",
  "  @formulas begin",
  "    conc ~ Normal(Dose/v * exp(-(cli/v)*t), sigma)",
  "  end",
  "end",
  "# a comment with / and = and 'quotes' and \"double quotes\" and 100%",
  sep = "\n"
)
model.b64 <- dsNoLimitsClient:::.ds_encode(model.src)

test_that("aggregate call shapes are legal and use only permitted methods", {
  expect_silent(ds_wire_check(call("nolimitsStatusDS"), "aggregate", DSNOLIMITS_AGGREGATE))
  expect_silent(ds_wire_check(call("nolimitsModelsDS"), "aggregate", DSNOLIMITS_AGGREGATE))
  expect_silent(ds_wire_check(call("nolimitsPrepInfoDS", "nolimits.prep"),
                              "aggregate", DSNOLIMITS_AGGREGATE))
  for (est in c("laplace", "focei", "ghq", "pooled", "mle", "map")) {
    for (lvl in 1:9) {
      cally <- call("nolimitsObjGradDS", "nolimits.prep", payload12,
                    est, as.numeric(lvl), 1)
      expect_silent(ds_wire_check(cally, "aggregate", DSNOLIMITS_AGGREGATE))
    }
  }
  # A large round number and a long payload still lex.
  expect_silent(ds_wire_check(
    call("nolimitsObjGradDS", "prep", dsNoLimitsClient:::.ds_num_encode(rnorm(200)),
         "laplace", 1, 99999),
    "aggregate", DSNOLIMITS_AGGREGATE
  ))
})

test_that("assign call shapes are legal and use only permitted methods", {
  expect_silent(ds_wire_check(call("nolimitsPrepDS", "D", "warfarin_pk", "ID", "t"),
                              "assign", DSNOLIMITS_ASSIGN))
  expect_silent(ds_wire_check(call("nolimitsPrepDS", "D", "one-compartment.iv_v2",
                                   "SUBJID", "TIME"),
                              "assign", DSNOLIMITS_ASSIGN))
  expect_silent(ds_wire_check(call("nolimitsPrepStringDS", "D", model.b64, "ID", "t"),
                              "assign", DSNOLIMITS_ASSIGN))
})

test_that("the base64url model payload survives the round trip it was checked as", {
  expect_match(model.b64, "^[A-Za-z0-9_-]+$")
  expect_identical(dsNoLimitsClient:::.ds_decode(model.b64), model.src)
})

test_that("ds_wire_check itself rejects the shapes we must never emit", {
  expect_error(ds_wire_check(call("dsNoLimits::nolimitsStatusDS"), "aggregate"), "backtick")
  expect_error(ds_wire_check(call("nolimitsPrepInfoDS", ""), "aggregate"), "empty string")
  expect_error(ds_wire_check(call("nolimitsPrepInfoDS", "D$x y"), "aggregate"),
               "illegal character")
  expect_error(ds_wire_check(call("nolimitsObjGradDS", "prep",
                                  as.call(list(as.name("nolimitsPrepInfoDS"), "prep"))),
                             "aggregate"), "nested DS calls")
  expect_error(ds_wire_check(call("someOtherDS", "prep"), "aggregate",
                             DSNOLIMITS_AGGREGATE), "not permitted methods")
})
