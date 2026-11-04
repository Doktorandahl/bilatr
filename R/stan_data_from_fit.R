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
#' @param fit A fit object, or a character vector of CSV paths.
#' @param stan_data The caller's `stan_data` argument.
#' @param caller The calling function's name, for messages.
#' @param required If `FALSE`, return `NULL` instead of erroring when
#'   nothing is found (for [icc_curves()], whose `stan_data` is optional).
#' @return A Stan data list, or `NULL` (only when `required = FALSE`).
#' @noRd
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
#' (unlabeled output) when `fit` is CSV paths or carries no Stan data.
#'
#' @param fit A fit object, or a character vector of CSV paths.
#' @return A character vector, or `NULL`.
#' @noRd
.attached_event_classes <- function(fit) {
  if (is.character(fit)) {
    return(NULL)
  }
  attr(attr(fit, "stan_data", exact = TRUE), "event_classes", exact = TRUE)
}

#' Which model a fit was made with
#'
#' Shared by the post-fit functions whose result depends on the model
#' ([diagnose_convergence()], [diagnose_and_extract_bilatr()],
#' [icc_curves()], [check_compositional_residuals()]). In order:
#'
#' 1. an explicit `stan_model` (canonicalized);
#' 2. `attr(fit, "stan_model")`, set by [fit_panel()]/[fit_dyad_ts()];
#' 3. `attr(stan_data, "stan_model")`, set on the Stan data the fit
#'    carries and on its `*_stan_data.rds` sidecar;
#' 4. for CSV paths, the `model = <name>` line of the first CSV's header,
#'    mapped through the registry (CmdStan names a model after its
#'    `.stan` file plus `_model`, e.g. `bilatr_alphanorm_gamma_model`);
#' 5. otherwise `"stable"`, with a `message()`.
#'
#' An explicit `stan_model` that contradicts a recorded one (2-4) is an
#' error naming both.
#'
#' @param fit A fit object, a `draws` object, or CSV paths.
#' @param stan_data The resolved Stan data list, or `NULL`.
#' @param stan_model The caller's `stan_model` argument, or `NULL`.
#' @param caller The calling function's name, for messages.
#' @return A canonical registered model name.
#' @noRd
.resolve_stan_model_for_fit <- function(fit, stan_data, stan_model, caller) {
  recorded <- NULL
  source <- NULL
  if (!is.character(fit) && !is.null(attr(fit, "stan_model", exact = TRUE))) {
    recorded <- attr(fit, "stan_model", exact = TRUE)
    source <- "the fit"
  } else if (!is.null(attr(stan_data, "stan_model", exact = TRUE))) {
    recorded <- attr(stan_data, "stan_model", exact = TRUE)
    source <- "`stan_data`"
  } else if (is.character(fit) && length(fit) > 0) {
    recorded <- .stan_model_from_csv_header(fit[1])
    source <- "the CmdStan CSV header"
  }
  if (!is.null(recorded)) {
    recorded <- .canonical_stan_model(recorded)
  }

  if (!is.null(stan_model)) {
    stan_model <- .canonical_stan_model(stan_model)
    if (!is.null(recorded) && !identical(recorded, stan_model)) {
      stop(
        caller, ": `stan_model = \"", stan_model, "\"` contradicts the model recorded in ",
        source, " (\"", recorded, "\"). Leave `stan_model` unset to use the recorded model.",
        call. = FALSE
      )
    }
    return(stan_model)
  }
  if (!is.null(recorded)) {
    return(recorded)
  }
  message(
    caller, ": no model is recorded on the fit, so `stan_model = \"stable\"` is assumed. ",
    "Pass `stan_model` if the fit used another model."
  )
  "stable"
}

#' Registered model name from a CmdStan CSV's `model = ` header line
#'
#' @param csv_file One CmdStan CSV path.
#' @return A registered model name, or `NULL` if the header has no
#'   `model` line, names a model that is not registered, or cannot be
#'   read.
#' @noRd
.stan_model_from_csv_header <- function(csv_file) {
  cmdstan_name <- tryCatch(
    .scan_stan_csv_header_and_skip(csv_file)$model,
    error = function(e) NA_character_
  )
  if (is.na(cmdstan_name)) {
    return(NULL)
  }
  files <- vapply(.bilatr_stan_models, `[[`, "", "file")
  registered <- paste0(tools::file_path_sans_ext(files), "_model")
  hit <- names(.bilatr_stan_models)[registered == cmdstan_name]
  if (length(hit) == 1) hit else NULL
}
