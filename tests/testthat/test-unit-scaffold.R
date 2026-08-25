test_that("the package is installed and its version matches the server package", {
  expect_true(nzchar(as.character(utils::packageVersion("dsNoLimitsClient"))))
  skip_if_not_installed("dsNoLimits")
  expect_identical(
    utils::packageVersion("dsNoLimitsClient"),
    utils::packageVersion("dsNoLimits")
  )
})
