# Tier 4: the client functions against DSLite servers running the real
# dsNoLimits server functions.

test_that("ds.nolimitsStatus returns one row per site", {
  skip_if_not_installed("dsNoLimits")
  conns <- ds_lite_conns(nlc_sites(), aggregate = DS_AGG,
                         options = ds_test_options())
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  out <- ds.nolimitsStatus(datasources = conns)

  expect_s3_class(out, "data.frame")
  expect_equal(out$site, c("a", "b", "c"))
  expect_named(out, c("site", "julia.ok", "julia.version", "nolimits.version",
                      "nolimitsr.version", "dsnolimits.version", "model.dir.set"))
  expect_type(out$julia.ok, "logical")
  # No registry configured in this harness.
  expect_false(any(out$model.dir.set))
})

test_that("ds.nolimitsModels lists the registered model on every site", {
  skip_if_not_installed("dsNoLimits")
  opts <- ds_model_dir_options(ds_test_options())
  conns <- ds_lite_conns(nlc_sites(), aggregate = DS_AGG, options = opts)
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  out <- ds.nolimitsModels(datasources = conns)

  expect_named(out, c("a", "b", "c"))
  expect_s3_class(out$a, "data.frame")
  expect_named(out$a, c("model", "hash"))
  expect_equal(out$a$model, "onecomp_iv")
  expect_match(out$a$hash, "^[0-9a-f]{64}$")
  expect_equal(attr(out, "common"), "onecomp_iv")
  # Same registered model on every site: the hash is identical across nodes.
  expect_identical(out$a$hash, out$b$hash)
  expect_identical(out$a$hash, out$c$hash)

  expect_true(ds.nolimitsStatus(datasources = conns)$model.dir.set[1L])
})

test_that("ds.nolimitsModels warns when no registry is configured", {
  skip_if_not_installed("dsNoLimits")
  conns <- ds_lite_conns(nlc_sites("a"), aggregate = DS_AGG,
                         options = ds_test_options())
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  expect_warning(out <- ds.nolimitsModels(datasources = conns),
                 "no model registry directory is configured")
  expect_equal(nrow(out$a), 0L)
  expect_named(out$a, c("model", "hash"))
})

test_that("the prep assign creates the symbol and prep info round-trips", {
  nlc_skip_no_julia()
  conns <- ds_lite_conns(nlc_sites(), aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options())
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  DSI::datashield.assign(
    conns, "nolimits.prep",
    call("nolimitsPrepStringDS", "D", dsNoLimitsClient:::.ds_encode(NLDS_TEST_MODEL),
         "ID", "t")
  )

  syms <- DSI::datashield.symbols(conns)
  expect_setequal(syms$a, c("D", "nolimits.prep"))
  expect_setequal(syms$c, c("D", "nolimits.prep"))

  info <- DSI::datashield.aggregate(conns, call("nolimitsPrepInfoDS", "nolimits.prep"))
  expect_named(info, c("a", "b", "c"))

  # Everything released by prep info must be a plain R value: no Julia proxy
  # can survive a serialize()/unserialize() round trip.
  expect_identical(unserialize(serialize(info, NULL)), info)

  expect_type(info$a$names, "character")
  expect_type(info$a$theta0, "double")
  expect_type(info$a$scale, "double")
  expect_equal(info$a$names, c("cl", "v", "omega", "sigma"))
  expect_equal(info$a$p, 4L)
  expect_equal(info$a$n.subjects, 6L)
  expect_equal(info$a$n.obs, 36L)

  # theta0 / names / scale are model-derived, so they are identical across
  # sites holding different data. That is what makes agree() sound.
  expect_identical(info$b$theta0, info$a$theta0)
  expect_identical(info$b$scale, info$a$scale)
  expect_identical(info$b$names, info$a$names)

  og <- DSI::datashield.aggregate(
    conns, call("nolimitsObjGradDS", "nolimits.prep",
                dsNoLimitsClient:::.ds_num_encode(info$a$theta0),
                "laplace", 1, 1)
  )
  expect_setequal(names(og$a), c("value", "grad", "natural", "finite", "n.subjects"))
  expect_length(og$a$grad, 4L)
  expect_true(og$a$finite)
  expect_identical(unserialize(serialize(og, NULL)), og)
  # Nothing released grows with n: three length-p vectors and three scalars.
  expect_true(all(vapply(og, function(r) length(unlist(r)), integer(1)) == 11L))
})

test_that("an unregistered model name surfaces through datashield.errors", {
  skip_if_not_installed("dsNoLimits")
  opts <- ds_model_dir_options(ds_test_options())
  conns <- ds_lite_conns(nlc_sites("a"), aggregate = DS_AGG, assign = DS_ASN,
                         options = opts)
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  expect_error(
    DSI::datashield.assign(conns, "nolimits.prep",
                           call("nolimitsPrepDS", "D", "no.such.model", "ID", "t")),
    "There are some DataSHIELD errors", fixed = TRUE
  )
  # A substring, never the exact string: Opal reshapes the message around it.
  expect_match(DSI::datashield.errors()$a,
               "is not a model registered on this server", fixed = TRUE)
})

test_that("the saturation gate fires at the canonical nfilter.glm", {
  nlc_skip_no_julia()
  # 4 parameters, 6 subjects: 4 > 0.33 * 6.
  conns <- ds_lite_conns(nlc_sites("a"), aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options(nfilter.glm = 0.33))
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  expect_error(
    DSI::datashield.assign(
      conns, "nolimits.prep",
      call("nolimitsPrepStringDS", "D",
           dsNoLimitsClient:::.ds_encode(NLDS_TEST_MODEL), "ID", "t")),
    "There are some DataSHIELD errors", fixed = TRUE
  )
  expect_match(DSI::datashield.errors()$a,
               "too many parameters for the number of subjects", fixed = TRUE)
  expect_false("nolimits.prep" %in% DSI::datashield.symbols(conns)$a)
})

test_that("string-mode prep is blocked when the site is not permissive", {
  skip_if_not_installed("dsNoLimits")
  conns <- ds_lite_conns(nlc_sites("a"), aggregate = DS_AGG, assign = DS_ASN,
                         options = ds_test_options(control = ""))
  on.exit(DSI::datashield.logout(conns), add = TRUE)

  expect_error(
    DSI::datashield.assign(
      conns, "nolimits.prep",
      call("nolimitsPrepStringDS", "D",
           dsNoLimitsClient:::.ds_encode(NLDS_TEST_MODEL), "ID", "t")),
    "There are some DataSHIELD errors", fixed = TRUE
  )
  expect_match(DSI::datashield.errors()$a, "BLOCKED", fixed = TRUE)
  expect_false("nolimits.prep" %in% DSI::datashield.symbols(conns)$a)
})
