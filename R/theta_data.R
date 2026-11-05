# Download, cache and verify the estimated theta trajectories, which are
# published as versioned releases of the public repository
# doktorandahl/bilatr-data (archived on Zenodo) rather than shipped in the
# package. Documented in ?theta_trajectories.

# ---------------------------------------------------------------------------
# Release coordinates. Each data release adds its version DOI and its Zenodo
# file URL here (from the record's `files` entries at
# https://zenodo.org/api/records/<id>); NA means "not available".
# ---------------------------------------------------------------------------

#' Data release coordinates
#'
#' `.BILATR_DATA_DEFAULT_VERSION` is the data version this package release
#' is pinned to, so results do not change silently when a new data release
#' appears. `.BILATR_DATA_ZENODO_URL` holds, per version, the URL of the
#' Zenodo copy of the release: either a directory of the release files
#' (ending in `/`) or a `.zip` archive of the repository at that tag (what
#' Zenodo's GitHub integration stores). `NA` means not available yet.
#' @noRd
.BILATR_DATA_DEFAULT_VERSION <- "v2026.11"
.BILATR_DATA_GITHUB_URL <- "https://github.com/doktorandahl/bilatr-data/releases/download/"
.BILATR_DATA_API_URL <- "https://api.github.com/repos/doktorandahl/bilatr-data/releases?per_page=100"
.BILATR_DATA_CONCEPT_DOI <- "10.5281/zenodo.23951818" # all versions
.BILATR_DATA_VERSION_DOI <- c(v2026.11 = "10.5281/zenodo.23951819")
.BILATR_DATA_ZENODO_URL <- c(
  v2026.11 = "https://zenodo.org/records/23951819/files/Doktorandahl/bilatr-data-v2026.11.zip"
)

#' Download (once) the estimated theta trajectories
#'
#' Returns one table of the estimated bilateral-relations trajectories
#' described in [theta_trajectories]. The first call for a table and data
#' version **downloads** it from the bilatr data release (a few MB),
#' verifies its checksum and saves it in a local cache, with a message
#' naming the version, the size and the cache path. Every later call reads
#' the cached copy, silently and without using the internet.
#'
#' The data are published as versioned releases of the public repository
#' <https://github.com/doktorandahl/bilatr-data>, archived on Zenodo, under
#' the CC BY 4.0 license. Each package release is pinned to one data
#' version (`version = NULL`), so the same code gives the same data until
#' you ask for another version. [bilatr_theta_versions()] lists the
#' versions.
#'
#' Files are cached in `tools::R_user_dir("bilatr", "data")`, one
#' directory per version, together with the release's checksum file.
#' Downloads are written to a temporary file and renamed only after the
#' checksum matches, so an interrupted download never leaves a corrupt
#' cache. A cached file is checked against the cached checksum file each
#' time it is read. [bilatr_clear_data_cache()] deletes the cache. Files
#' are downloaded from GitHub first, and from the Zenodo copy of the
#' release if GitHub fails.
#'
#' The checksum is SHA-256 when `tools::sha256sum()` is available (R
#' 4.5.0 or later), and MD5 otherwise; every release carries both
#' checksum files.
#'
#' @param type `"directed"` (one series per sender-receiver pair) or
#'   `"undirected"` (one series per unordered pair).
#' @param version The data version, such as `"v2026.11"`. `NULL` (the
#'   default) uses the version this package release is pinned to.
#'   `"latest"` asks the GitHub API for the newest release, which needs
#'   the internet even when that version is cached.
#' @param refresh If `TRUE`, download again even when a cached copy
#'   exists.
#' @param quiet If `TRUE`, no messages, not even on the first download.
#'
#' @return A tibble with one row per dyad-year and the columns described
#'   in [theta_trajectories].
#'
#' @section Errors:
#' Without an internet connection (and no cached copy), when a file is
#' missing from the release, or when a checksum does not match,
#' `download_bilatr_theta()` stops with an error that names the URL(s) tried and the
#' cache path. It never returns a partial or unverified table.
#'
#' @seealso [theta_trajectories] for what the data are, [theta_sample] for
#'   a small extract that ships with the package.
#' @export
#' @examplesIf interactive()
#' theta <- download_bilatr_theta("directed") # downloads once, then reads the cache
#' subset(theta, dyad %in% c("NOR_SWE", "SWE_NOR") & year >= 2020)
#'
#' # the versions available, and which are cached
#' bilatr_theta_versions()
download_bilatr_theta <- function(type = c("directed", "undirected"), version = NULL,
                         refresh = FALSE, quiet = FALSE) {
  type <- rlang::arg_match(type)
  .check_data_flag(refresh, "refresh")
  .check_data_flag(quiet, "quiet")
  version <- .bilatr_data_resolve_version(version)

  path <- .bilatr_data_fetch(paste0("theta_", type, ".rds"), version, refresh, quiet)
  tryCatch(
    readRDS(path),
    error = function(e) {
      stop("download_bilatr_theta(): could not read the cached file `", path, "` (",
           conditionMessage(e), "). Use `refresh = TRUE` to download it again.",
           call. = FALSE)
    }
  )
}

#' List the available theta data versions
#'
#' Asks the GitHub API which data releases exist, and marks the version
#' this package release is pinned to and the versions already cached.
#' Without an internet connection it lists only the cached versions, with
#' a message.
#'
#' @return A tibble with one row per version, newest first: `version`,
#'   `published` (the release date, `NA` when offline), `default` (the
#'   version [download_bilatr_theta()] uses when `version = NULL`) and `cached`.
#' @seealso [download_bilatr_theta()], [bilatr_clear_data_cache()].
#' @export
#' @examplesIf interactive()
#' bilatr_theta_versions()
bilatr_theta_versions <- function() {
  cached <- .bilatr_data_cached_versions()
  releases <- tryCatch(.bilatr_data_releases(), error = function(e) e)
  if (inherits(releases, "error")) {
    message("bilatr_theta_versions(): ", conditionMessage(releases),
            " Listing cached versions only.")
    out <- tibble::tibble(version = cached, published = as.Date(rep(NA, length(cached))))
  } else {
    out <- releases
    extra <- setdiff(cached, out$version)
    if (length(extra) > 0) {
      out <- dplyr::bind_rows(out, tibble::tibble(version = extra, published = as.Date(NA)))
    }
  }
  out$default <- out$version == .BILATR_DATA_DEFAULT_VERSION
  out$cached <- out$version %in% cached
  out[order(out$version, decreasing = TRUE, method = "radix"), , drop = FALSE]
}

#' Delete cached theta data
#'
#' Deletes the files [download_bilatr_theta()] cached under
#' `tools::R_user_dir("bilatr", "data")`.
#'
#' @param version A data version such as `"v2026.11"`, or `NULL` (the
#'   default) for every cached version.
#' @return The deleted version directories, invisibly.
#' @seealso [download_bilatr_theta()], [bilatr_clear_model_cache()] for compiled
#'   models.
#' @export
#' @examples
#' # see what would be deleted, without deleting anything
#' file.path(tools::R_user_dir("bilatr", "data"), "v2026.11")
#' @examplesIf interactive()
#' bilatr_clear_data_cache("v2026.11")
bilatr_clear_data_cache <- function(version = NULL) {
  root <- .bilatr_data_cache_dir()
  if (is.null(version)) {
    version <- .bilatr_data_cached_versions()
  } else {
    .check_data_version(version)
  }
  dirs <- file.path(root, version)
  dirs <- dirs[dir.exists(dirs)]
  unlink(dirs, recursive = TRUE)
  message("bilatr_clear_data_cache(): deleted ", length(dirs), " cached data version(s)",
          if (length(dirs) > 0) paste0(" (", paste(basename(dirs), collapse = ", "), ")"),
          " from `", root, "`.")
  invisible(dirs)
}

# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------

#' Root of the data cache; `getOption("bilatr.data_cache_dir")` overrides it
#' @noRd
.bilatr_data_cache_dir <- function() {
  getOption("bilatr.data_cache_dir", tools::R_user_dir("bilatr", "data"))
}

#' Base URL of the GitHub release downloads, ending in `/`;
#' `getOption("bilatr.data_base_url")` overrides it (tests use `file://`)
#' @noRd
.bilatr_data_base_url <- function() {
  base <- getOption("bilatr.data_base_url", .BILATR_DATA_GITHUB_URL)
  if (!grepl("/$", base)) base <- paste0(base, "/")
  base
}

#' @noRd
.check_data_flag <- function(x, arg) {
  if (!(is.logical(x) && length(x) == 1 && !is.na(x))) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
}

#' A data version is used in a file path, so only `v` followed by
#' letters, digits, dots, underscores or hyphens is accepted
#' @noRd
.check_data_version <- function(version) {
  if (!(is.character(version) && length(version) == 1 && !is.na(version) &&
    grepl("^v[0-9][0-9A-Za-z._-]*$", version) && !grepl("..", version, fixed = TRUE))) {
    stop("`version` must be a data version such as \"", .BILATR_DATA_DEFAULT_VERSION,
         "\", \"latest\" or NULL.", call. = FALSE)
  }
}

#' Resolve `NULL` and `"latest"` to a concrete version
#' @noRd
.bilatr_data_resolve_version <- function(version) {
  if (is.null(version)) {
    return(.BILATR_DATA_DEFAULT_VERSION)
  }
  if (identical(version, "latest")) {
    releases <- tryCatch(.bilatr_data_releases(), error = function(e) {
      stop("download_bilatr_theta(version = \"latest\"): ", conditionMessage(e),
           " Use `version = NULL` for the pinned version (", .BILATR_DATA_DEFAULT_VERSION,
           "), or name a cached version.", call. = FALSE)
    })
    if (nrow(releases) == 0) {
      stop("download_bilatr_theta(version = \"latest\"): the data repository has no releases.",
           call. = FALSE)
    }
    return(releases$version[[1]])
  }
  .check_data_version(version)
  version
}

#' Fetch the GitHub API's release list as text (one request)
#' @noRd
.bilatr_data_api_get <- function() {
  url <- getOption("bilatr.data_api_url", .BILATR_DATA_API_URL)
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)
  ok <- tryCatch(
    {
      utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
      TRUE
    },
    error = function(e) FALSE,
    warning = function(w) FALSE
  )
  if (!ok || !file.exists(tmp)) {
    stop("Could not list the data releases at <", url, "> (offline, or the GitHub API ",
         "rate limit was reached).", call. = FALSE)
  }
  paste(readLines(tmp, warn = FALSE), collapse = "\n")
}

#' The published `v*` releases, newest first
#'
#' Parses the GitHub API's release list with jsonlite. Drafts, pre-releases
#' and tags that are not data versions are skipped.
#' @return A tibble with `version` and `published` (a `Date`).
#' @noRd
.bilatr_data_releases <- function() {
  json <- .bilatr_data_api_get()
  releases <- tryCatch(jsonlite::fromJSON(json, simplifyVector = FALSE),
                       error = function(e) NULL)
  if (!is.list(releases) || !is.null(names(releases))) {
    # The API reports errors as an object, e.g. {"message": "API rate limit exceeded ..."}
    api_msg <- if (is.list(releases) && is.character(releases$message)) {
      paste0(" (", releases$message[[1]], ")")
    } else {
      ""
    }
    stop("Unexpected response from the GitHub API when listing the data releases",
         api_msg, ".", call. = FALSE)
  }
  string <- function(x) if (is.character(x) && length(x) == 1 && !is.na(x)) x else NULL
  rows <- lapply(releases, function(r) {
    if (!is.list(r)) {
      return(NULL)
    }
    tag <- string(r$tag_name)
    published <- string(r$published_at)
    if (is.null(tag) || is.null(published) || !grepl("^v[0-9]", tag) ||
      !identical(r$prerelease, FALSE) || isTRUE(r$draft)) {
      return(NULL)
    }
    tibble::tibble(version = tag, published = as.Date(substr(published, 1, 10)))
  })
  out <- dplyr::bind_rows(
    tibble::tibble(version = character(0), published = as.Date(character(0))),
    rows
  )
  out[order(out$published, out$version, decreasing = TRUE, method = "radix"), , drop = FALSE]
}

#' Versions present in the cache
#' @noRd
.bilatr_data_cached_versions <- function() {
  root <- .bilatr_data_cache_dir()
  if (!dir.exists(root)) {
    return(character(0))
  }
  v <- list.dirs(root, full.names = FALSE, recursive = FALSE)
  sort(v[grepl("^v[0-9]", v)], method = "radix")
}

#' Checksum algorithm and file: SHA-256 when base R has it, else MD5
#' @noRd
.bilatr_data_checksum <- function() {
  if (exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE)) {
    list(file = "SHA256SUMS", fun = getExportedValue("tools", "sha256sum"))
  } else {
    list(file = "MD5SUMS", fun = tools::md5sum)
  }
}

#' Parse a `SHA256SUMS`/`MD5SUMS` file into a named vector (name = file)
#' @noRd
.bilatr_data_read_sums <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  hash <- sub("\\s.*$", "", lines)
  file <- sub("^\\*", "", sub("^\\S+\\s+", "", lines))
  stats::setNames(tolower(hash), file)
}

#' Where each file of a version can be downloaded from, in order of preference
#' @noRd
.bilatr_data_sources <- function(file, version) {
  out <- list(list(url = paste0(.bilatr_data_base_url(), version, "/", file), archive = FALSE))
  zenodo <- unname(.BILATR_DATA_ZENODO_URL[version])
  if (length(zenodo) == 1 && !is.na(zenodo)) {
    if (grepl("\\.zip$", zenodo)) {
      out <- c(out, list(list(url = zenodo, archive = TRUE)))
    } else {
      if (!grepl("/$", zenodo)) zenodo <- paste0(zenodo, "/")
      out <- c(out, list(list(url = paste0(zenodo, file), archive = FALSE)))
    }
  }
  out
}

#' Download `file` to `dest`, trying each source in turn
#'
#' A `.zip` source is an archive of the whole release; `file` is extracted
#' from it. Returns the URL used, or `NULL` when every source failed.
#' @noRd
.bilatr_data_download <- function(file, version, dest, quiet) {
  for (src in .bilatr_data_sources(file, version)) {
    if (!quiet) message("Downloading ", file, " (", version, ") from <", src$url, ">.")
    target <- if (src$archive) tempfile(fileext = ".zip") else dest
    ok <- tryCatch(
      {
        utils::download.file(src$url, target, mode = "wb", quiet = TRUE)
        file.exists(target)
      },
      error = function(e) FALSE,
      warning = function(w) FALSE
    )
    if (ok && src$archive) {
      ok <- .bilatr_data_unzip_one(target, file, dest)
      unlink(target)
    }
    if (ok) {
      return(src$url)
    }
    unlink(target)
  }
  NULL
}

#' Extract the archive entry whose base name is `file` to `dest`
#' @noRd
.bilatr_data_unzip_one <- function(zip, file, dest) {
  entries <- tryCatch(utils::unzip(zip, list = TRUE)$Name, error = function(e) character(0))
  entry <- entries[basename(entries) == file]
  if (length(entry) != 1) {
    return(FALSE)
  }
  exdir <- tempfile("bilatr_unzip_")
  on.exit(unlink(exdir, recursive = TRUE), add = TRUE)
  utils::unzip(zip, files = entry, exdir = exdir, junkpaths = TRUE)
  file.exists(file.path(exdir, file)) && file.copy(file.path(exdir, file), dest, overwrite = TRUE)
}

#' Return the cached path of `file`, downloading and verifying it if needed
#' @noRd
.bilatr_data_fetch <- function(file, version, refresh, quiet) {
  dir <- file.path(.bilatr_data_cache_dir(), version)
  path <- file.path(dir, file)
  sums <- .bilatr_data_checksum()
  sums_path <- file.path(dir, sums$file)

  if (file.exists(path) && !refresh) {
    if (file.exists(sums_path)) {
      expected <- .bilatr_data_read_sums(sums_path)[file]
      if (!is.na(expected) && !identical(tolower(unname(sums$fun(path))), unname(expected))) {
        stop("download_bilatr_theta(): the cached file `", path, "` does not match its checksum in `",
             sums_path, "`. Use `refresh = TRUE` to download it again.", call. = FALSE)
      }
    }
    return(path)
  }

  old <- options(timeout = max(300, getOption("timeout")))
  on.exit(options(old), add = TRUE)
  if (!dir.exists(dir) && !dir.create(dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("download_bilatr_theta(): could not create the cache directory `", dir, "`.", call. = FALSE)
  }

  # 1. The release's checksum file, written through a temp file.
  sums_tmp <- tempfile(pattern = paste0(".", sums$file, "_"), tmpdir = dir)
  on.exit(unlink(sums_tmp), add = TRUE)
  if (is.null(.bilatr_data_download(sums$file, version, sums_tmp, quiet = TRUE))) {
    stop(.bilatr_data_fail_msg(sums$file, version, path,
                               "could not download the release's checksum file"), call. = FALSE)
  }
  expected <- .bilatr_data_read_sums(sums_tmp)[file]
  if (is.na(expected)) {
    stop(.bilatr_data_fail_msg(file, version, path, paste0(
      "`", file, "` is not listed in the release's ", sums$file, " file"
    )), call. = FALSE)
  }

  # 2. The data file, verified before it is renamed into place.
  tmp <- tempfile(pattern = paste0(".", file, "_"), tmpdir = dir)
  on.exit(unlink(tmp), add = TRUE)
  url <- .bilatr_data_download(file, version, tmp, quiet)
  if (is.null(url)) {
    stop(.bilatr_data_fail_msg(file, version, path, "the download failed"), call. = FALSE)
  }
  actual <- tolower(unname(sums$fun(tmp)))
  if (!identical(actual, unname(expected))) {
    stop("download_bilatr_theta(): the file downloaded from <", url, "> does not match its checksum (",
         sums$file, "), so it was not cached. Nothing was written to `", path, "`.",
         call. = FALSE)
  }
  unlink(c(path, sums_path))
  if (!file.rename(tmp, path) || !file.rename(sums_tmp, sums_path)) {
    unlink(path)
    stop("download_bilatr_theta(): could not write to the cache `", path, "`.", call. = FALSE)
  }
  if (!quiet) {
    message("Downloaded ", file, " (data version ", version, ", ",
            .bilatr_data_format_size(file.size(path)), ") and cached it at `", path,
            "`. Later calls read this copy without downloading.")
  }
  path
}

#' File size for messages, such as "2.3 MB"
#' @noRd
.bilatr_data_format_size <- function(bytes) {
  if (bytes >= 1e6) sprintf("%.1f MB", bytes / 1e6) else sprintf("%.0f kB", max(1, bytes / 1e3))
}

#' Error message naming every URL tried and the cache path
#' @noRd
.bilatr_data_fail_msg <- function(file, version, path, what) {
  urls <- vapply(.bilatr_data_sources(file, version), `[[`, "", "url")
  paste0(
    "download_bilatr_theta(): ", what, " (data version ", version, ").\n",
    "Tried: ", paste0("<", urls, ">", collapse = ", "), ".\n",
    "Cache path: `", path, "`.\n",
    "Check your internet connection, or that the version exists ",
    "(`bilatr_theta_versions()`)."
  )
}
