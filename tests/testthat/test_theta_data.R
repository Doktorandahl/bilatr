# bilatr_theta(), bilatr_theta_versions() and bilatr_clear_data_cache(),
# against a file:// fixture release (no network), with the GitHub API
# mocked.

# A fixture release: <root>/<version>/ with small .rds files, both checksum
# files and a manifest. `corrupt = TRUE` writes wrong checksums.
make_theta_release <- function(root, version = "v2026.11", corrupt = FALSE) {
  dir <- file.path(root, version)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  d <- theta_sample[theta_sample$type == "directed", -1]
  u <- theta_sample[theta_sample$type == "undirected", -1]
  if (version != "v2026.11") d$mean <- d$mean + 1 # tell versions apart
  saveRDS(d, file.path(dir, "theta_directed.rds"), compress = "xz")
  saveRDS(u, file.path(dir, "theta_undirected.rds"), compress = "xz")
  writeLines(sprintf('{"data_version": "%s"}', version), file.path(dir, "manifest.json"))
  files <- c("manifest.json", "theta_directed.rds", "theta_undirected.rds")
  md5 <- unname(tools::md5sum(file.path(dir, files)))
  sha <- if (exists("sha256sum", envir = asNamespace("tools"))) {
    unname(getExportedValue("tools", "sha256sum")(file.path(dir, files)))
  } else {
    strrep("0", 64)
  }
  if (corrupt) {
    md5 <- strrep("0", 32)
    sha <- strrep("0", 64)
  }
  writeLines(paste0(sha, "  ", files), file.path(dir, "SHA256SUMS"))
  writeLines(paste0(md5, "  ", files), file.path(dir, "MD5SUMS"))
  list(dir = dir, directed = d, undirected = u)
}

local_theta_release <- function(..., env = parent.frame()) {
  root <- withr::local_tempdir("bilatr release ", .local_envir = env)
  cache <- withr::local_tempdir("bilatr data cache ", .local_envir = env)
  withr::local_options(
    bilatr.data_base_url = paste0("file://", root, "/"),
    bilatr.data_cache_dir = cache,
    .local_envir = env
  )
  c(make_theta_release(root, ...), root = root, cache = cache)
}

# The GitHub API's release list, as returned for this repository
api_json <- '[
  {"url": "x", "tag_name": "v2027.02-rc", "name": "rc", "draft": false,
   "prerelease": true, "published_at": "2027-02-01T10:00:00Z",
   "assets": [{"name": "theta_directed.rds", "created_at": "2027-02-01T10:00:00Z"}]},
  {"url": "x", "tag_name": "v2027.01", "name": "January 2027", "draft": false,
   "prerelease": false, "published_at": "2027-01-15T10:00:00Z",
   "assets": [{"name": "theta_directed.rds", "created_at": "2027-01-15T10:00:00Z"}]},
  {"url": "x", "tag_name": "notes", "name": "not data", "draft": false,
   "prerelease": false, "published_at": "2026-12-01T10:00:00Z", "assets": []},
  {"url": "x", "tag_name": "v2026.11", "name": "November 2026", "draft": false,
   "prerelease": false, "published_at": "2026-11-04T10:00:00Z",
   "assets": [{"name": "theta_directed.rds", "created_at": "2026-11-04T10:00:00Z"}]}
]'

test_that("theta_sample has the documented shape", {
  expect_s3_class(theta_sample, "tbl_df")
  expect_identical(
    names(theta_sample),
    c("type", "dyad", "year", "mean", "5%", "95%", "filtered_mean", "filtered_5%", "filtered_95%")
  )
  expect_identical(nrow(theta_sample), 1038L)
  expect_identical(as.vector(table(theta_sample$type)), c(690L, 348L))
  expect_length(unique(theta_sample$dyad[theta_sample$type == "directed"]), 20)
  expect_length(unique(theta_sample$dyad[theta_sample$type == "undirected"]), 10)
  expect_true(all(c("USA_CHN", "CHN_USA", "NOR_SWE", "SWE_NOR") %in% theta_sample$dyad))
  expect_false(anyNA(theta_sample))
})

test_that("bilatr_theta() downloads, verifies and caches, then reads the cache", {
  fx <- local_theta_release()
  n <- 0
  real_download <- .bilatr_data_download
  local_mocked_bindings(.bilatr_data_download = function(...) {
    n <<- n + 1
    real_download(...)
  })

  expect_message(x <- bilatr_theta("directed"), "Downloading theta_directed.rds")
  expect_identical(x, fx$directed)
  expect_identical(n, 2) # checksum file + data file
  sums <- .bilatr_data_checksum()$file
  expect_true(file.exists(file.path(fx$cache, "v2026.11", "theta_directed.rds")))
  expect_true(file.exists(file.path(fx$cache, "v2026.11", sums)))
  expect_setequal(list.files(file.path(fx$cache, "v2026.11"), all.files = TRUE, no.. = TRUE),
                  c("theta_directed.rds", sums))

  # second call: from the cache, no download
  expect_no_message(y <- bilatr_theta("directed"))
  expect_identical(y, fx$directed)
  expect_identical(n, 2)

  # refresh = TRUE downloads again
  expect_identical(bilatr_theta("directed", refresh = TRUE, quiet = TRUE), fx$directed)
  expect_identical(n, 4)

  expect_identical(bilatr_theta("undirected", quiet = TRUE), fx$undirected)
})

test_that("MD5SUMS is used where tools::sha256sum() is unavailable", {
  fx <- local_theta_release()
  local_mocked_bindings(.bilatr_data_checksum = function() {
    list(file = "MD5SUMS", fun = tools::md5sum)
  })
  expect_identical(bilatr_theta("directed", quiet = TRUE), fx$directed)
  expect_true(file.exists(file.path(fx$cache, "v2026.11", "MD5SUMS")))
  expect_false(file.exists(file.path(fx$cache, "v2026.11", "SHA256SUMS")))
})

test_that("quiet = TRUE silences the progress messages", {
  local_theta_release()
  expect_no_message(bilatr_theta("directed", quiet = TRUE))
})

test_that("a checksum mismatch errors and leaves nothing in the cache", {
  fx <- local_theta_release(corrupt = TRUE)
  expect_error(bilatr_theta("directed", quiet = TRUE), "does not match its checksum")
  expect_length(list.files(file.path(fx$cache, "v2026.11"), all.files = TRUE, no.. = TRUE), 0)
})

test_that("a corrupted cached file errors and points to refresh = TRUE", {
  fx <- local_theta_release()
  bilatr_theta("directed", quiet = TRUE)
  path <- file.path(fx$cache, "v2026.11", "theta_directed.rds")
  writeLines("garbage", path)
  expect_error(bilatr_theta("directed"), "refresh = TRUE")
  expect_identical(bilatr_theta("directed", refresh = TRUE, quiet = TRUE), fx$directed)
})

test_that("an unreachable URL errors with the URL and the cache path", {
  fx <- local_theta_release()
  withr::local_options(bilatr.data_base_url = paste0("file://", fx$root, "/does-not-exist/"))
  err <- expect_error(bilatr_theta("directed", quiet = TRUE), "could not download")
  expect_match(conditionMessage(err), "does-not-exist/v2026.11/", fixed = TRUE)
  expect_match(conditionMessage(err), file.path(fx$cache, "v2026.11", "theta_directed.rds"),
               fixed = TRUE)
  expect_length(list.files(file.path(fx$cache, "v2026.11"), all.files = TRUE, no.. = TRUE), 0)
})

test_that("a file missing from the release errors informatively", {
  fx <- local_theta_release()
  unlink(file.path(fx$dir, "theta_undirected.rds"))
  expect_error(bilatr_theta("undirected", quiet = TRUE), "the download failed")
})

test_that("the Zenodo archive is used when GitHub fails", {
  fx <- local_theta_release()
  zip <- file.path(fx$root, "bilatr-data-v2026.11.zip")
  withr::with_dir(fx$root, utils::zip(zip, "v2026.11", flags = "-r9Xq"))
  withr::local_options(bilatr.data_base_url = paste0("file://", fx$root, "/does-not-exist/"))
  local_mocked_bindings(.BILATR_DATA_ZENODO_URL = c(v2026.11 = paste0("file://", zip)))
  expect_identical(bilatr_theta("directed", quiet = TRUE), fx$directed)
})

test_that("bilatr_clear_data_cache() deletes one version or all of them", {
  fx <- local_theta_release()
  make_theta_release(fx$root, "v2027.01")
  bilatr_theta("directed", quiet = TRUE)
  bilatr_theta("directed", version = "v2027.01", quiet = TRUE)
  expect_identical(.bilatr_data_cached_versions(), c("v2026.11", "v2027.01"))

  expect_message(gone <- bilatr_clear_data_cache("v2026.11"), "deleted 1")
  expect_identical(basename(gone), "v2026.11")
  expect_identical(.bilatr_data_cached_versions(), "v2027.01")

  expect_message(bilatr_clear_data_cache(), "deleted 1")
  expect_identical(.bilatr_data_cached_versions(), character(0))
  expect_true(dir.exists(fx$cache)) # the root itself is kept
  expect_message(bilatr_clear_data_cache(), "deleted 0")
})

test_that("versions are validated before they are used in a path", {
  local_theta_release()
  expect_error(bilatr_theta(version = "../etc"), "data version")
  expect_error(bilatr_theta(version = c("v1", "v2")), "data version")
  expect_error(bilatr_clear_data_cache("v2026.11/.."), "data version")
  expect_error(bilatr_theta(refresh = NA), "TRUE or FALSE")
  expect_error(bilatr_theta("both"), "must be one of")
})

test_that("version = 'latest' resolves to the newest published v* release", {
  fx <- local_theta_release()
  make_theta_release(fx$root, "v2027.01")
  local_mocked_bindings(.bilatr_data_api_get = function() api_json)
  expect_identical(.bilatr_data_resolve_version("latest"), "v2027.01")
  x <- bilatr_theta("directed", version = "latest", quiet = TRUE)
  expect_identical(x$mean, fx$directed$mean + 1)
  expect_identical(.bilatr_data_resolve_version(NULL), .BILATR_DATA_DEFAULT_VERSION)
})

test_that("version = 'latest' errors clearly when the API cannot be reached", {
  local_theta_release()
  withr::local_options(bilatr.data_api_url = "file:///does-not-exist/releases.json")
  expect_error(bilatr_theta(version = "latest"), "Could not list the data releases")
  local_mocked_bindings(.bilatr_data_api_get = function() "<html>rate limited</html>")
  expect_error(bilatr_theta(version = "latest"), "no releases")
})

test_that("bilatr_theta_versions() marks the default and the cached versions", {
  fx <- local_theta_release()
  bilatr_theta("directed", quiet = TRUE)
  local_mocked_bindings(.bilatr_data_api_get = function() api_json)
  v <- bilatr_theta_versions()
  expect_identical(v$version, c("v2027.01", "v2026.11"))
  expect_identical(v$published, as.Date(c("2027-01-15", "2026-11-04")))
  expect_identical(v$default, c(FALSE, TRUE))
  expect_identical(v$cached, c(FALSE, TRUE))
})

test_that("bilatr_theta_versions() lists only cached versions when offline", {
  fx <- local_theta_release()
  bilatr_theta("directed", quiet = TRUE)
  local_mocked_bindings(.bilatr_data_api_get = function() stop("offline"))
  expect_message(v <- bilatr_theta_versions(), "cached versions only")
  expect_identical(v$version, "v2026.11")
  expect_true(is.na(v$published))
  expect_true(v$cached)
  expect_true(v$default)
})

test_that("the checksum file parser handles both sum formats", {
  f <- withr::local_tempfile()
  writeLines(c("ABC123  theta_directed.rds", "def456 *theta_undirected.rds", ""), f)
  expect_identical(.bilatr_data_read_sums(f),
                   c(theta_directed.rds = "abc123", theta_undirected.rds = "def456"))
})

test_that("the published data release downloads and matches theta_sample", {
  skip_on_cran()
  skip_if_offline("github.com")
  skip_if(is.na(.BILATR_DATA_VERSION_DOI[[.BILATR_DATA_DEFAULT_VERSION]]),
          "the data release coordinates are still pending")
  withr::local_options(bilatr.data_cache_dir = withr::local_tempdir())
  x <- bilatr_theta("directed", quiet = TRUE)
  expect_identical(nrow(x), 355821L)
  s <- theta_sample[theta_sample$type == "directed", -1]
  expect_identical(x[x$dyad %in% s$dyad, ], s)
})
