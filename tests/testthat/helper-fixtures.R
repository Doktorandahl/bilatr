make_fake_events <- function(n = 400, seed = 42, years = 2015:2019) {
  set.seed(seed)
  countries <- c("USA", "CHN", "RUS")
  dates <- format(
    sample(
      seq(as.Date(paste0(min(years), "-01-01")), as.Date(paste0(max(years), "-12-31")), by = "day"),
      n,
      replace = TRUE
    ),
    "%Y%m%d"
  )
  events <- tibble::tibble(
    Actor1CountryCode = sample(countries, n, replace = TRUE),
    Actor2CountryCode = sample(countries, n, replace = TRUE),
    SQLDATE = as.integer(dates),
    EventCode = sample(bilatr::cameo_lookup$CAMEOEVENTCODE, n, replace = TRUE)
  )
  events[events$Actor1CountryCode != events$Actor2CountryCode, ]
}
