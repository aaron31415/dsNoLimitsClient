#' @title List the models registered on each server
#' @description Lists the pre-registered model names each data owner has made
#'   available for `model.mode = "registered"` fits. A model that is not listed
#'   by *every* site cannot be fitted federatively.
#' @details One aggregate round trip. Sites that have configured no registry
#'   directory return an empty character vector and are reported with a warning.
#'
#'   Server function called: `nolimitsModelsDS`.
#' @param datasources a list of [DSI::DSConnection-class] objects obtained after
#'   login. If not specified the default set of connections is used: see
#'   [DSI::datashield.connections_default()].
#' @return a named list with one data frame per server, each with a `model`
#'   column (the registered model names) and a `hash` column (the SHA-256 of each
#'   model's canonical source). The list carries an attribute `"common"` holding
#'   the model names available on every server. An auditor compares the hashes
#'   across servers to confirm a given model is byte-identical everywhere.
#' @author Manuel Huth
#' @examples
#' \dontrun{
#' ds.nolimitsModels(datasources = conns)
#' }
#' @export
ds.nolimitsModels <- function(datasources = NULL) {
  datasources <- .setDatasources(datasources)

  res <- .agg(datasources, call("nolimitsModelsDS"), "model registry listing")
  .check_complete(res, names(datasources), "model registry listing")

  unset <- names(res)[!vapply(res, function(x) isTRUE(.field1(x, "model.dir.set", FALSE)), logical(1))]
  if (length(unset)) {
    warning("no model registry directory is configured on: ",
            paste(sort(unset), collapse = ", "), call. = FALSE)
  }

  out <- lapply(res, function(x) {
    models <- as.character(x$models)
    hashes <- if (is.null(x$hashes)) rep(NA_character_, length(models)) else
      as.character(x$hashes)
    if (length(hashes) != length(models)) {
      hashes <- rep(NA_character_, length(models))
    }
    data.frame(model = models, hash = hashes,
               row.names = NULL, stringsAsFactors = FALSE)
  })
  attr(out, "common") <- Reduce(intersect, lapply(out, function(d) d$model))
  out
}
