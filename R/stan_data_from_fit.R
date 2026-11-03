#' Find the Stan data list that belongs to a fit
#'
#' Shared by the extractors and diagnostics that take `fit, stan_data =
#' NULL`. In order:
#'
#' 1. an explicitly supplied `stan_data` is returned as-is;
#' 2. otherwise, `attr(fit, "stan_data")`, which [fit_panel()]/
#'    [fit_dyad_ts()] attach to the fit;
#' 3. otherwise, if `fit` is a character vector of CmdStan CSV paths, the
#'    single `*_stan_data.rds` file in their common directory (written by
#'    `fit_*()` when given an `output_dir`), read with a `message()`
#'    naming the file;
#' 4. otherwise an error, or `NULL` when `required = FALSE`.
#'
#' Stan data lists saved before 0.12.0 are returned unchanged (their
#' `dyad_ids` keep `actor_a`/`actor_b`).
#'
#' @param fit A fit object, or a character vector of CSV paths.
#' @param stan_data The caller's `stan_data` argument.
#' @param caller The calling function's name, for messages.
#' @param required If `FALSE`, return `NULL` instead of erroring when
#'   nothing is found (for [icc_curves()], whose `stan_data` is optional).
#' @return A Stan data list, or `NULL` (only when `required = FALSE`).
#' @keywords internal
.stan_data_from_fit <- function(fit, stan_data, caller, required = TRUE) {
  if (!is.null(stan_data)) {
    return(stan_data)
  }
  if (!is.character(fit)) {
    attached <- attr(fit, "stan_data", exact = TRUE)
    if (!is.null(attached)) {
      return(attached)
    }
    problem <- "`fit` carries no \"stan_data\" attribute"
  } else {
    dirs <- unique(normalizePath(dirname(fit), mustWork = FALSE))
    if (length(dirs) != 1) {
      problem <- "the CSV files are not all in one directory"
    } else {
      found <- list.files(dirs, pattern = "_stan_data\\.rds$", full.names = TRUE)
      if (length(found) == 1) {
        message(caller, ": reading `stan_data` from ", found, ".")
        return(readRDS(found))
      }
      problem <- if (length(found) == 0) {
        sprintf("there is no *_stan_data.rds file next to the CSV files (in %s)", dirs)
      } else {
        sprintf(
          "there are %d *_stan_data.rds files next to the CSV files (%s), so which one belongs to them is ambiguous",
          length(found), paste(basename(found), collapse = ", ")
        )
      }
    }
  }
  if (!required) {
    return(NULL)
  }
  stop(
    caller, ": `stan_data` was not supplied and could not be found: ", problem, ". ",
    "Pass `stan_data`, the list from prepare_stan_data() that the model was fit to. ",
    "Fits from fit_panel()/fit_dyad_ts() carry it as attr(fit, \"stan_data\"), and ",
    "save it as <output_basename>_stan_data.rds in their `output_dir`.",
    call. = FALSE
  )
}

#' Event classes of the Stan data attached to a fit
#'
#' For [extract_alpha()]/[extract_mu_intercept()]'s `event_classes = NULL`
#' default: `attr(attr(fit, "stan_data"), "event_classes")`, or `NULL`
#' (unlabelled output) when `fit` is CSV paths or carries no Stan data.
#'
#' @param fit A fit object, or a character vector of CSV paths.
#' @return A character vector, or `NULL`.
#' @keywords internal
.attached_event_classes <- function(fit) {
  if (is.character(fit)) {
    return(NULL)
  }
  attr(attr(fit, "stan_data", exact = TRUE), "event_classes", exact = TRUE)
}
