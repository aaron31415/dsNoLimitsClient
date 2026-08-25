## The wire codecs are duplicated in dsNoLimits (server) and dsNoLimitsClient.
## Duplication is deliberate - the packages live on different machines - so this
## file is the only thing stopping the two copies drifting apart.

# Nasty but transmissible doubles. Denormals are DELIBERATELY excluded: the
# smallest subnormal 5e-324 (0x1p-1074) does not survive sprintf("%a") on all
# platforms and is documented as out of contract.
REFERENCE_DOUBLES <- c(
  0, -0, 1, -1, pi, exp(1), 1 / 3, 2 / 3, 1e-300, 1e300,
  .Machine$double.xmax, -.Machine$double.xmax,
  .Machine$double.xmin, .Machine$double.eps,
  Inf, -Inf, NaN, 1234567.891011, -1e-17
)

test_that(".ds_num_encode/.ds_num_decode round-trip exactly", {
  s <- dsNoLimitsClient:::.ds_num_encode(REFERENCE_DOUBLES)
  back <- dsNoLimitsClient:::.ds_num_decode(s)
  expect_identical(back, REFERENCE_DOUBLES)
  # Signed zero survives too, which identical() does not check by default.
  expect_identical(1 / back[2L], -Inf)
})

test_that("the numeric payload is legal inside a wire string literal", {
  s <- dsNoLimitsClient:::.ds_num_encode(REFERENCE_DOUBLES)
  expect_match(s, "^[A-Za-z0-9_$+,.:;-]+$")
})

test_that("NA_real_ is not silently transmissible", {
  s <- dsNoLimitsClient:::.ds_num_encode(NA_real_)
  expect_identical(s, "NA")
  expect_warning(dsNoLimitsClient:::.ds_num_decode(s))
})

test_that(".ds_encode/.ds_decode round-trip arbitrary text over base64url", {
  sources <- c(
    "conc ~ Normal(Dose/v * exp(-(cli/v)*t), sigma)",
    "@Model begin\n  @fixedEffects begin\n    v = 8.0\n  end\nend\n",
    "a <- \"quoted\"; b = 'single' # comment / slash = equals + plus % pct",
    "unicode: äöü µg/L"
  )
  for (src in sources) {
    enc <- dsNoLimitsClient:::.ds_encode(src)
    expect_match(enc, "^[A-Za-z0-9_-]+$")
    expect_identical(dsNoLimitsClient:::.ds_decode(enc), enc2utf8(src))
  }
})

test_that("the client and server copies of the codecs agree", {
  skip_if_not_installed("dsNoLimits")
  expect_identical(
    dsNoLimits:::.ds_num_encode(REFERENCE_DOUBLES),
    dsNoLimitsClient:::.ds_num_encode(REFERENCE_DOUBLES)
  )
  wire <- dsNoLimitsClient:::.ds_num_encode(REFERENCE_DOUBLES)
  expect_identical(dsNoLimits:::.ds_num_decode(wire), REFERENCE_DOUBLES)

  src <- "conc ~ Normal(Dose/v, sigma)\n# with a comment"
  expect_identical(dsNoLimits:::.ds_encode(src), dsNoLimitsClient:::.ds_encode(src))
  expect_identical(
    dsNoLimits:::.ds_decode(dsNoLimitsClient:::.ds_encode(src)),
    src
  )
  # The empty string must NOT encode to "" (an Opal lexer error): both copies
  # guard it to NULL, and the two copies still agree.
  expect_null(dsNoLimits:::.ds_encode(""))
  expect_null(dsNoLimitsClient:::.ds_encode(""))
  expect_identical(dsNoLimits:::.ds_encode(""), dsNoLimitsClient:::.ds_encode(""))
})

test_that(".ds_encode('') returns NULL in the client copy", {
  expect_null(dsNoLimitsClient:::.ds_encode(""))
})
