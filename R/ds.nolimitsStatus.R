#' @title Per-site readiness of the NoLimits server provisioning
#' @description Reports, for every connected server, whether Julia and NoLimits
#'   are usable and whether a model registry directory has been configured. This
#'   is the diagnostic to run first when a fit fails with a provisioning error;
#'   it touches no data and never errors on a badly provisioned site.
#' @details One aggregate round trip.
#'
#'   Server function called: `nolimitsStatusDS`.
#' @param datasources a list of [DSI::DSConnection-class] objects obtained after
#'   login. If not specified the default set of connections is used: see
#'   [DSI::datashield.connections_default()].
#' @return a data frame with one row per server and columns `site`, `julia.ok`,
#'   `julia.version`, `nolimits.version`, `nolimitsr.version`,
#'   `dsnolimits.version` and `model.dir.set`, printed and returned invisibly.
#' @author Manuel Huth
#' @examples
#' \dontrun{
#' ds.nolimitsStatus(datasources = conns)
#' }
#' @export
ds.nolimitsStatus <- function(datasources = NULL) {
  datasources <- .setDatasources(datasources)

  res <- .agg(datasources, call("nolimitsStatusDS"), "status")
  .check_complete(res, names(datasources), "status")

  out <- data.frame(
    site = names(res),
    julia.ok = vapply(res, .field1, logical(1), "julia.ok", NA),
    julia.version = vapply(res, function(x) as.character(.field1(x, "julia.version", NA_character_)), character(1)),
    nolimits.version = vapply(res, function(x) as.character(.field1(x, "nolimits.version", NA_character_)), character(1)),
    nolimitsr.version = vapply(res, function(x) as.character(.field1(x, "nolimitsr.version", NA_character_)), character(1)),
    dsnolimits.version = vapply(res, function(x) as.character(.field1(x, "dsnolimits.version", NA_character_)), character(1)),
    model.dir.set = vapply(res, .field1, logical(1), "model.dir.set", NA),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  print(out)
  invisible(out)
}
