# A minimal concrete DSConnection so client-side validation can be exercised
# with a `datasources` argument that passes .setDatasources(). Any DSI call on
# it fails with "unable to find an inherited method", which is unmistakably
# different from a validation abort - that is what makes the Tier 3 tests
# meaningful: a validation error that fires means no request was attempted.
methods::setClass("FakeDSConnection", contains = "DSConnection")

fake_datasources <- function(n = 2L) {
  conns <- replicate(n, methods::new("FakeDSConnection"), simplify = FALSE)
  names(conns) <- paste0("site", seq_len(n))
  conns
}
