# Pre-flight legality check for every call this client can emit, from
# wire-protocol.md section 12. It is stricter than DSLite's own validator,
# which recurses on the wrong object and therefore only ever checks the root
# node of the expression.

# Legal characters inside a transmitted string literal: the intersection of
# DSLite's t_STRING/t_DBLSTRING and datashield4j v2's R_STRING/R_DBLSTRING.
.WIRE_STRING_OK <- "^[A-Za-z0-9_$+,.:;-]+$"
.is_empty_arg <- function(a) is.symbol(a) && !nzchar(as.character(a))

ds_wire_check <- function(cally, type = c("aggregate", "assign"), methods = NULL) {
  type <- match.arg(type)
  stopifnot(is.language(cally) || is.character(cally))
  text <- if (is.language(cally)) paste(deparse(cally), collapse = "\n") else cally

  has_ns <- function(e) is.call(e) &&
    (identical(e[[1L]], quote(`::`)) || any(vapply(as.list(e), has_ns, logical(1))))
  if (grepl("`", text, fixed = TRUE)) stop("backtick in expression: ", text)
  if (if (is.language(cally)) has_ns(cally) else grepl("::", text, fixed = TRUE)) {
    stop("namespaced call cannot cross the parser: ", text)
  }

  if (is.language(cally)) {
    strs <- character(0)
    collect <- function(e) {
      if (is.call(e)) {
        for (a in as.list(e)[-1L]) if (!.is_empty_arg(a)) collect(a)
      } else if (is.character(e)) {
        strs <<- c(strs, e)
      }
    }
    collect(cally)
    for (s in strs) {
      if (!nzchar(s)) stop("empty string argument is a lexer error on Opal")
      if (!grepl(.WIRE_STRING_OK, s)) stop("illegal character in string literal: ", s)
    }
  }

  root <- DSLite::testParse(text)

  if (type == "aggregate" && !inherits(root, "FunctionNode")) {
    stop("an aggregate expression must invoke a function at the root: ", text)
  }

  found <- character(0)
  visitor <- list(visit = function(node) {
    if (inherits(node, "FunctionNode")) found <<- c(found, node$name)
  })
  root$accept(visitor)
  found <- unique(found)

  if (!is.null(methods)) {
    absent <- setdiff(found, methods)
    if (length(absent)) stop("not permitted methods: ", paste(absent, collapse = ", "))
  }
  if (sum(grepl("DS[0-9]*$", found)) > 1L) {
    stop("nested DS calls fail on DSLite: ", paste(found, collapse = ", "))
  }

  invisible(found)
}

# The methods dsNoLimits declares in inst/DATASHIELD.
DSNOLIMITS_AGGREGATE <- c("nolimitsStatusDS", "nolimitsModelsDS",
                          "nolimitsPrepInfoDS", "nolimitsObjGradDS",
                          "nolimitsObjGradDpDS", "nolimitsNaturalDS",
                          "nolimitsEmQDS", "nolimitsEmQDpDS",
                          "nolimitsSaemStatsDS", "nolimitsSaemStatsDpDS",
                          "nolimitsSaemUpdatesDS")
DSNOLIMITS_ASSIGN <- c("nolimitsPrepDS", "nolimitsPrepStringDS",
                       "nolimitsEmEStepDS", "nolimitsSaemMstepDS")

# The same two vectors as DSLite method registrations, for helper-dslite.R.
DS_AGG <- stats::setNames(paste0("dsNoLimits::", DSNOLIMITS_AGGREGATE),
                          DSNOLIMITS_AGGREGATE)
DS_ASN <- stats::setNames(paste0("dsNoLimits::", DSNOLIMITS_ASSIGN),
                          DSNOLIMITS_ASSIGN)
