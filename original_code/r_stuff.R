library(tidyverse)
library(foreach)
library(httr)

download_and_summarize_gdelt <- function(
  date,
  local_folder = "data/gdelt_data",
  type = c('daily', 'monthly', 'yearly')
) {
  header <- get_gdelt_column_names()

  if (type == 'daily') {
    f <- paste0(
      "http://data.gdeltproject.org/events/",
      str_remove_all(date, "-"),
      ".export.CSV.zip"
    )
    f2 <- paste0(local_folder, "/", str_remove_all(date, "-"), ".zip")
  } else if (type == 'monthly') {
    year_month <- str_remove_all(date, "-") %>% str_sub(1, 6)
    f <- paste0(
      "http://data.gdeltproject.org/events/",
      year_month,
      ".zip"
    )
    f2 <- paste0(local_folder, "/", year_month, ".zip")
  } else if (type == 'yearly') {
    year <- str_remove_all(date, "-") %>% str_sub(1, 4)
    f <- paste0(
      "http://data.gdeltproject.org/events/",
      year,
      ".zip"
    )
    f2 <- paste0(local_folder, "/", year, ".zip")
  } else {
    stop("Invalid type specified. Choose from 'daily', 'monthly', or 'yearly'.")
  }

  url_check <- httr::HEAD(f)

  if (url_check$status_code != 200L) {
    warning("GDELT data not available for ", date)
    return(NULL)
  } else {
    try(download.file(f, destfile = f2, mode = "wb", quiet = TRUE))
  }

  out <- read_delim(
    unz(
      f2,
      unzip(f2, list = TRUE)$Name[1]
    ),
    num_threads = 1,
    col_names = header,
    progress = FALSE
  )

  gs_neg10 <- out %>%
    filter(GoldsteinScale == -10)

  if (nrow(gs_neg10) > 0) {
    saveRDS(
      gs_neg10,
      paste0("data/goldstein_neg10/", str_remove_all(date, "-"), ".rds")
    )
  }

  if (type == 'daily') {
    out <- out %>%
      filter(!is.na(Actor1CountryCode) & !is.na(Actor2CountryCode)) %>%
      group_by(Actor1CountryCode, Actor2CountryCode) %>%
      summarize(
        QC1 = sum(QuadClass == 1),
        QC2 = sum(QuadClass == 2),
        QC3 = sum(QuadClass == 3),
        QC4 = sum(QuadClass == 4),
        total_events = n(),
        mean_tone = mean(AvgTone, na.rm = TRUE),
        prop_positive = mean(AvgTone > 0, na.rm = TRUE),
        prop_negative = mean(AvgTone < 0, na.rm = TRUE),
        mean_goldstein = mean(GoldsteinScale, na.rm = TRUE),
        prop_positive_goldstein = mean(GoldsteinScale > 0, na.rm = TRUE),
        prop_negative_goldstein = mean(GoldsteinScale < 0, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      mutate(date = date)
  } else {
    out <- out %>%
      filter(!is.na(Actor1CountryCode) & !is.na(Actor2CountryCode)) %>%
      group_by(Actor1CountryCode, Actor2CountryCode, SQLDATE) %>%
      summarize(
        QC1 = sum(QuadClass == 1),
        QC2 = sum(QuadClass == 2),
        QC3 = sum(QuadClass == 3),
        QC4 = sum(QuadClass == 4),
        total_events = n(),
        mean_tone = mean(AvgTone, na.rm = TRUE),
        prop_positive = mean(AvgTone > 0, na.rm = TRUE),
        prop_negative = mean(AvgTone < 0, na.rm = TRUE),
        mean_goldstein = mean(GoldsteinScale, na.rm = TRUE),
        prop_positive_goldstein = mean(GoldsteinScale > 0, na.rm = TRUE),
        prop_negative_goldstein = mean(GoldsteinScale < 0, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      mutate(date = ymd(SQLDATE)) %>%
      select(-SQLDATE)
  }

  saveRDS(
    out,
    paste0(local_folder, "/gdelt_", str_remove_all(date, "-"), ".rds")
  )
  file.remove(f2)
}


get_gdelt_column_names <- function() {
  ## Gdelt column names
  c(
    "GLOBALEVENTID",
    "SQLDATE",
    "MonthYear",
    "Year",
    "FractionDate",
    "Actor1Code",
    "Actor1Name",
    "Actor1CountryCode",
    "Actor1KnownGroupCode",
    "Actor1EthnicCode",
    "Actor1Religion1Code",
    "Actor1Religion2Code",
    "Actor1Type1Code",
    "Actor1Type2Code",
    "Actor1Type3Code",
    "Actor2Code",
    "Actor2Name",
    "Actor2CountryCode",
    "Actor2KnownGroupCode",
    "Actor2EthnicCode",
    "Actor2Religion1Code",
    "Actor2Religion2Code",
    "Actor2Type1Code",
    "Actor2Type2Code",
    "Actor2Type3Code",
    "IsRootEvent",
    "EventCode",
    "EventBaseCode",
    "EventRootCode",
    "QuadClass",
    "GoldsteinScale",
    "NumMentions",
    "NumSources",
    "NumArticles",
    "AvgTone",
    "Actor1Geo_Type",
    "Actor1Geo_FullName",
    "Actor1Geo_CountryCode",
    "Actor1Geo_ADM1Code",
    "Actor1Geo_Lat",
    "Actor1Geo_Long",
    "Actor1Geo_FeatureID",
    "Actor2Geo_Type",
    "Actor2Geo_FullName",
    "Actor2Geo_CountryCode",
    "Actor2Geo_ADM1Code",
    "Actor2Geo_Lat",
    "Actor2Geo_Long",
    "Actor2Geo_FeatureID",
    "ActionGeo_Type",
    "ActionGeo_FullName",
    "ActionGeo_CountryCode",
    "ActionGeo_ADM1Code",
    "ActionGeo_Lat",
    "ActionGeo_Long",
    "ActionGeo_FeatureID",
    "DATEADDED"
  )
}

relevant_actors <- function() {
  c("GOV", "MIL", "SPY")
}


download_gdelt_raw_zip <- function(
  date,
  local_folder = "data/gdelt_raw",
  type = c('daily', 'monthly', 'yearly')
) {
  if (type == 'daily') {
    f <- paste0(
      "http://data.gdeltproject.org/events/",
      str_remove_all(date, "-"),
      ".export.CSV.zip"
    )
    f2 <- paste0(local_folder, "/", str_remove_all(date, "-"), ".zip")
  } else if (type == 'monthly') {
    year_month <- str_remove_all(date, "-") %>% str_sub(1, 6)
    f <- paste0(
      "http://data.gdeltproject.org/events/",
      year_month,
      ".zip"
    )
    f2 <- paste0(local_folder, "/", year_month, ".zip")
  } else if (type == 'yearly') {
    year <- str_remove_all(date, "-") %>% str_sub(1, 4)
    f <- paste0(
      "http://data.gdeltproject.org/events/",
      year,
      ".zip"
    )
    f2 <- paste0(local_folder, "/", year, ".zip")
  } else {
    stop("Invalid type specified. Choose from 'daily', 'monthly', or 'yearly'.")
  }

  url_check <- httr::HEAD(f)

  if (url_check$status_code != 200L) {
    warning("GDELT data not available for ", date)
    return(NULL)
  } else {
    try(download.file(f, destfile = f2, mode = "wb", quiet = TRUE))
  }
}


summarize_gdelt <- function(
  file,
  destination_folder = "data/gdelt_summarized",
  save_file = TRUE,
  return_data = FALSE
) {
  header <- get_gdelt_column_names()

  rel_act <- relevant_actors()

  out <- read_delim(
    unz(
      file,
      unzip(file, list = TRUE)$Name[1]
    ),
    num_threads = 1,
    col_names = header,
    progress = FALSE
  )

  dest_file <- str_remove_all(basename(file), ".zip")

  ## Filte relevant events
  out <- out %>%
    filter(
      !is.na(Actor1CountryCode) & !is.na(Actor2CountryCode),
      Actor1CountryCode != Actor2CountryCode,
      ((Actor1Type1Code %in%
        rel_act |
        Actor1Type2Code %in% rel_act |
        Actor1Type3Code %in% rel_act) &
        (Actor2Type1Code %in%
          rel_act |
          Actor2Type2Code %in% rel_act |
          Actor2Type3Code %in% rel_act))
    )

  out <- out %>%
    group_by(Actor1CountryCode, Actor2CountryCode, SQLDATE) %>%
    summarize(
      QC1 = sum(QuadClass == 1),
      QC2 = sum(QuadClass == 2),
      QC3 = sum(QuadClass == 3),
      QC4 = sum(QuadClass == 4),
      total_events = n(),
      mean_tone = mean(AvgTone, na.rm = TRUE),
      prop_positive = mean(AvgTone > 0, na.rm = TRUE),
      prop_negative = mean(AvgTone < 0, na.rm = TRUE),
      mean_goldstein = mean(GoldsteinScale, na.rm = TRUE),
      prop_positive_goldstein = mean(GoldsteinScale > 0, na.rm = TRUE),
      prop_negative_goldstein = mean(GoldsteinScale < 0, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(date = ymd(SQLDATE)) %>%
    select(-SQLDATE)
  if (save_file) {
    saveRDS(
      out,
      paste0(destination_folder, "/gdelt_", dest_file, ".rds")
    )
  }
  if (return_data) {
    return(out)
  }
}


get_actor_combos <- function(file) {
  header <- get_gdelt_column_names()
  rel_act <- relevant_actors()

  out <- read_delim(
    unz(
      file,
      unzip(file, list = TRUE)$Name[1]
    ),
    num_threads = 1,
    col_names = header,
    progress = FALSE,
    show_col_types = FALSE
  )

  out %>%
    filter(
      !is.na(Actor1CountryCode) & !is.na(Actor2CountryCode)
    ) %>%
    select(Actor1Type1Code, Actor2Type1Code) %>%
    count(Actor1Type1Code, Actor2Type1Code) %>%
    mutate(date = str_remove_all(basename(file), ".zip"))
}

summaries_to_dyad_period <- function(
  files,
  resolution = c('monthly', 'yearly'),
  out_name,
  directed = c(TRUE, FALSE)
) {
  all_data <- foreach(f = files, .final = bind_rows) %do%
    {
      readRDS(f)
    }

  if (directed) {
    all_data <- all_data %>%
      mutate(dyad = paste(Actor1CountryCode, Actor2CountryCode, sep = "_"))
  } else {
    all_data <- all_data %>%
      mutate(
        dyad = pmap_chr(
          list(Actor1CountryCode, Actor2CountryCode),
          ~ paste(sort(c(...)), collapse = "_")
        )
      )
  }

  all_data <- all_data %>%
    select(
      dyad,
      date,
      QC1,
      QC2,
      QC3,
      QC4,
      total_events,
      mean_tone,
      prop_positive,
      prop_negative,
      mean_goldstein,
      prop_positive_goldstein,
      prop_negative_goldstein
    )

  if (resolution == 'monthly') {
    all_data <- all_data %>%
      group_by(dyad, year = year(date), month = month(date)) %>%
      summarize(
        QC1 = sum(QC1),
        QC2 = sum(QC2),
        QC3 = sum(QC3),
        QC4 = sum(QC4),
        mean_tone = weighted.mean(mean_tone, w = total_events, na.rm = TRUE),
        mean_goldstein = weighted.mean(
          mean_goldstein,
          w = total_events,
          na.rm = TRUE
        ),
        prop_positive = weighted.mean(
          prop_positive,
          w = total_events,
          na.rm = TRUE
        ),
        prop_negative = weighted.mean(
          prop_negative,
          w = total_events,
          na.rm = TRUE
        ),
        prop_positive_goldstein = weighted.mean(
          prop_positive_goldstein,
          w = total_events,
          na.rm = TRUE
        ),
        prop_negative_goldstein = weighted.mean(
          prop_negative_goldstein,
          w = total_events,
          na.rm = TRUE
        ),
        total_events = sum(total_events)
      ) %>%
      ungroup()
  } else if (resolution == 'yearly') {
    all_data <- all_data %>%
      group_by(dyad, year = year(date)) %>%
      summarize(
        QC1 = sum(QC1),
        QC2 = sum(QC2),
        QC3 = sum(QC3),
        QC4 = sum(QC4),
        mean_tone = weighted.mean(mean_tone, w = total_events, na.rm = TRUE),
        mean_goldstein = weighted.mean(
          mean_goldstein,
          w = total_events,
          na.rm = TRUE
        ),
        prop_positive = weighted.mean(
          prop_positive,
          w = total_events,
          na.rm = TRUE
        ),
        prop_negative = weighted.mean(
          prop_negative,
          w = total_events,
          na.rm = TRUE
        ),
        prop_positive_goldstein = weighted.mean(
          prop_positive_goldstein,
          w = total_events,
          na.rm = TRUE
        ),
        prop_negative_goldstein = weighted.mean(
          prop_negative_goldstein,
          w = total_events,
          na.rm = TRUE
        ),
        total_events = sum(total_events)
      ) %>%
      ungroup()
  }

  saveRDS(all_data, out_name)
}


extract_all_relevant_gdelt <- function(
  file,
  destination_folder = "data/gdelt_filtered",
  save_file = TRUE,
  return_data = FALSE
) {
  header <- get_gdelt_column_names()

  rel_act <- relevant_actors()

  out <- read_delim(
    unz(
      file,
      unzip(file, list = TRUE)$Name[1]
    ),
    num_threads = 1,
    col_names = header,
    progress = FALSE
  )

  dest_file <- str_remove_all(basename(file), ".zip")

  ## Filte relevant events
  out <- out %>%
    filter(
      !is.na(Actor1CountryCode) & !is.na(Actor2CountryCode),
      Actor1CountryCode != Actor2CountryCode,
      ((Actor1Type1Code %in%
        rel_act |
        Actor1Type2Code %in% rel_act |
        Actor1Type3Code %in% rel_act) &
        (Actor2Type1Code %in%
          rel_act |
          Actor2Type2Code %in% rel_act |
          Actor2Type3Code %in% rel_act))
    )

  if (save_file) {
    saveRDS(
      out,
      paste0(destination_folder, "/gdelt_filtered_", dest_file, ".rds")
    )
  }
  if (return_data) {
    return(out)
  }
}


relevents_to_dyad_period <- function(
  data,
  resolution = c('monthly', 'yearly'),
  typeevents = c('quadclass', 'pentaclass'),
  out_name,
  directed = c(TRUE, FALSE)
) {
  if (directed) {
    data <- data %>%
      mutate(dyad = paste(Actor1CountryCode, Actor2CountryCode, sep = "_"))
  } else {
    data <- data %>%
      mutate(
        dyad = pmap_chr(
          list(Actor1CountryCode, Actor2CountryCode),
          ~ paste(sort(c(...)), collapse = "_")
        )
      )
  }

  data <- data %>%
    mutate(date = ymd(SQLDATE))

  if (typeevents == 'quadclass') {
    data <- data %>%
      select(
        dyad,
        date,
        QuadClass
      ) %>%
      rename(event_type = QuadClass)
  } else if (typeevents == 'pentaclass') {
    data <- data %>%
      select(
        dyad,
        date,
        PentaClass
      ) %>%
      rename(event_type = PentaClass)
  } else {
    stop(
      "Invalid typeevents specified. Choose from 'quadclass' or 'pentaclass'."
    )
  }

  if (typeevents == 'quadclass') {
    prefix <- "QC"
  } else if (typeevents == 'pentaclass') {
    prefix <- "PC"
  }

  column_order <- c(
    "dyad",
    "year",
    "month",
    paste0("QC", 1:4),
    paste0("PC", 0:4),
    "total_events"
  )

  if (resolution == 'monthly') {
    data2 <- data %>%
      group_by(dyad, year = year(date), month = month(date)) %>%
      count(event_type) %>%
      pivot_wider(
        names_from = event_type,
        values_from = n,
        values_fill = list(n = 0)
      ) %>%
      ungroup() %>%
      rename_with(~ paste0(prefix, .x), -c(dyad, year, month))

    total_events <- data %>%
      group_by(dyad, year = year(date), month = month(date)) %>%
      count() %>%
      ungroup() %>%
      rename(total_events = n)

    data2 <- data2 %>%
      left_join(total_events, by = c("dyad", "year", "month"))
  } else if (resolution == 'yearly') {
    data2 <- data %>%
      group_by(dyad, year = year(date)) %>%
      count(event_type) %>%
      pivot_wider(
        names_from = event_type,
        values_from = n,
        values_fill = list(n = 0)
      ) %>%
      ungroup() %>%
      rename_with(~ paste0(prefix, .x), -c(dyad, year))

    total_events <- data %>%
      group_by(dyad, year = year(date)) %>%
      count() %>%
      ungroup() %>%
      rename(total_events = n)

    data2 <- data2 %>%
      left_join(total_events, by = c("dyad", "year"))
  }

  data2 <- data2 %>%
    select(any_of(column_order))

  saveRDS(data2, out_name)
}


grouped_events_to_dyad_period <- function(
  data,
  resolution = c('monthly', 'yearly'),
  grouping_var,
  directed = c(TRUE, FALSE),
  reference_category = NULL,
  reference_hostile = NULL
) {
  if (directed) {
    data <- data %>%
      mutate(dyad = paste(Actor1CountryCode, Actor2CountryCode, sep = "_"))
  } else {
    data <- data %>%
      mutate(
        dyad = pmap_chr(
          list(Actor1CountryCode, Actor2CountryCode),
          ~ paste(sort(c(...)), collapse = "_")
        )
      )
  }

  data <- data %>%
    mutate(date = ymd(SQLDATE))

  data <- data %>%
    select(
      dyad,
      date,
      !!sym(grouping_var)
    ) %>%
    rename(event_type = !!sym(grouping_var))

  if (
    is.null(reference_category) ||
      !(reference_category %in% unique(data$event_type)) &&
        (is.null(reference_hostile) ||
          !(reference_hostile %in% unique(data$event_type)))
  ) {
    column_order <- c(
      "dyad",
      "year",
      "month",
      paste0("EventClass_", sort(unique(data$event_type))),
      "total_events"
    )
  } else if (
    is.null(reference_hostile) ||
      !(reference_hostile %in% unique(data$event_type))
  ) {
    column_order <- c(
      "dyad",
      "year",
      "month",
      paste0("EventClass_", reference_category),
      paste0(
        "EventClass_",
        sort(setdiff(unique(data$event_type), reference_category))
      ),
      "total_events"
    )
  } else if (
    is.null(reference_category) ||
      !(reference_category %in% unique(data$event_type))
  ) {
    column_order <- c(
      "dyad",
      "year",
      "month",
      paste0(
        "EventClass_",
        sort(setdiff(unique(data$event_type), reference_hostile))
      ),
      paste0("EventClass_", reference_hostile),
      "total_events"
    )
  } else {
    column_order <- c(
      "dyad",
      "year",
      "month",
      paste0("EventClass_", reference_category),
      paste0(
        "EventClass_",
        sort(setdiff(
          unique(data$event_type),
          c(reference_category, reference_hostile)
        ))
      ),
      paste0("EventClass_", reference_hostile),
      "total_events"
    )
  }

  if (resolution == 'monthly') {
    data2 <- data %>%
      group_by(dyad, year = year(date), month = month(date)) %>%
      count(event_type) %>%
      pivot_wider(
        names_from = event_type,
        values_from = n,
        values_fill = list(n = 0)
      ) %>%
      ungroup() %>%
      rename_with(~ paste0("EventClass_", .x), -c(dyad, year, month))

    total_events <- data %>%
      group_by(dyad, year = year(date), month = month(date)) %>%
      count() %>%
      ungroup() %>%
      rename(total_events = n)

    data2 <- data2 %>%
      left_join(total_events, by = c("dyad", "year", "month"))
  } else if (resolution == 'yearly') {
    data2 <- data %>%
      group_by(dyad, year = year(date)) %>%
      count(event_type) %>%
      pivot_wider(
        names_from = event_type,
        values_from = n,
        values_fill = list(n = 0)
      ) %>%
      ungroup() %>%
      rename_with(~ paste0("EventClass_", .x), -c(dyad, year))

    total_events <- data %>%
      group_by(dyad, year = year(date)) %>%
      count() %>%
      ungroup() %>%
      rename(total_events = n)

    data2 <- data2 %>%
      left_join(total_events, by = c("dyad", "year"))
  }

  data2 <- data2 %>%
    select(any_of(column_order))

  return(data2)
}

library(cmdstanr)
library(posterior)
library(tidyverse)

extract_theta_trajectory <- function(
  obj_path,
  actor1,
  actor2,
  data,
  model_name = NULL,
  data_set = NULL,
  years = 2000:2025
) {
  obj <- readRDS(obj_path)

  obj_data <- obj$summary(
    variables = c("theta"),
    mean,
    ~ quantile(.x, probs = c(0.05, 0.5, 0.95))
  ) %>%
    mutate(
      model = model_name,
      actor1 = actor1,
      actor2 = actor2,
      data_set = data_set,
      time_index = 1:n()
    )

  dyad_data <- make_stan_data_irt(
    data,
    years = years,
    actors = c(actor1, actor2),
    type = model_name,
    data_set = data_set,
    return_for = 'summary'
  ) %>%
    mutate(
      time_index = 1:n(),
      dyad = paste(actor1, actor2, sep = '-'),
      date = ymd(paste(year, month, '01', sep = '-'))
    ) %>%
    right_join(obj_data, by = 'time_index')

  return(dyad_data)
}


extract_theta_hierarchial <- function(
  obj_path,
  data,
  years = 2000:2025,
  dyads = NULL,
  resolution = c('monthly', 'yearly'),
  min_n_events = 1,
  data_set = c('gdelt', 'icews'),
  obj_ids = NULL
) {
  resolution <- match.arg(resolution, choices = c('monthly', 'yearly'))
  data_set <- match.arg(data_set, choices = c('gdelt', 'icews'))

  obj <- readRDS(obj_path)

  if (is.character(data)) {
    data <- readRDS(data)
  }

  obj_data <- obj$summary(
    variables = c("theta"),
    mean,
    ~ quantile(.x, probs = c(0.05, 0.5, 0.95))
  ) %>%
    mutate(variable = str_remove(str_remove(variable, "theta\\["), "\\]")) %>%
    separate(variable, into = c("dyad_id", "time_index"), sep = ",") %>%
    mutate(
      dyad_id = as.integer(dyad_id),
      time_index = as.integer(time_index)
      #model = model_name,
      #data_set = data_set
    )

  dyad_identifiers <- make_dyad_ids(
    data,
    years = years,
    resolution = resolution,
    min_n_events = min_n_events,
    data_set = data_set
  )

  obj_data <- obj_data %>%
    left_join(dyad_identifiers, by = 'dyad_id') %>%
    mutate(dyad2 = paste(sort(str_split_1(dyad, "_")), collapse = '_'))

  if (!is.null(dyads)) {
    obj_data <- obj_data %>%
      filter(dyad2 %in% dyads)
  }
  if (!is.null(obj_ids)) {
    obj_data <- obj_data %>%
      bind_cols(obj_ids)
  }

  return(obj_data)
}


event_type <- 'PentaClass'
reference_category <- '044'
reference_hostile <- '19'
temporal_resolution <- 'yearly'
years <- 1990:2025
model <- 'flexevents'
directed <- 'undirected'
start_year <- 1990
min_n_events <- 100
base_dir <- 'model_output'
cores <- 1
save_alpha <- TRUE

summarize_theta_trajectories <- function(
  event_type,
  reference_category,
  reference_hostile,
  temporal_resolution,
  model,
  directed,
  start_year,
  end_year = 2025,
  min_n_events,
  base_dir = 'model_output',
  cores = 1,
  save_alpha = TRUE,
  data_path = NULL
) {
  dir_short <- ifelse(directed == 'directed', 'dir', 'undir')

  obj_dir <- paste0(
    base_dir,
    '/',
    model,
    '_',
    event_type,
    '_ref',
    reference_category,
    '_c',
    substr(temporal_resolution, 1, 1),
    '_',
    dir_short,
    '_',
    start_year,
    '_',
    min_n_events
  )

  save_file_name <- paste0(
    substr(dir_short, 1, 1),
    'd',
    substr(temporal_resolution, 1, 1),
    '_',
    event_type,
    '_',
    model,
    '_',
    event_type,
    '_ref',
    reference_category,
    '_',
    min_n_events,
    '_',
    start_year,
    '.rds'
  )

  obj <- readRDS(
    list.files(obj_dir, full.names = TRUE) %>%
      str_subset(".rds")
  )

  if (is.null(data_path)) {
    data_path <- "data/gdelt_filtered_full.rds"
  }

  years <- start_year:end_year

  data <- readRDS(data_path)

  dirlogic <- ifelse(directed == 'directed', TRUE, FALSE)

  data <- grouped_events_to_dyad_period(
    data = data,
    resolution = temporal_resolution,
    grouping_var = event_type,
    directed = dirlogic,
    reference_category = reference_category,
    reference_hostile = reference_hostile
  )

  data <- data %>%
    filter(year >= start_year) %>%
    group_by(dyad) %>%
    mutate(n_events = sum(total_events)) %>%
    filter(n_events >= min_n_events) %>%
    ungroup()

  if (temporal_resolution == 'yearly') {
    temp_skeleton <- expand_grid(
      dyad = unique(data$dyad),
      year = years
    )

    dyad_identifiers <- data %>%
      right_join(
        temp_skeleton,
        by = c('dyad', 'year')
      ) %>%
      arrange(dyad, year) %>%
      select(dyad, year) %>%
      group_by(dyad) %>%
      mutate(time_index = row_number(), dyad_id = cur_group_id()) %>%
      ungroup()
  } else if (temporal_resolution == 'monthly') {
    temp_skeleton <- expand_grid(
      dyad = unique(data$dyad),
      year = years,
      month = 1:12
    )
    dyad_identifiers <- data %>%
      right_join(
        temp_skeleton,
        by = c('dyad', 'year', 'month')
      ) %>%
      arrange(dyad, year, month) %>%
      select(dyad, year, month) %>%
      group_by(dyad) %>%
      mutate(time_index = row_number(), dyad_id = cur_group_id()) %>%
      ungroup()
  }

  obj_data <- obj$summary(
    variables = c("theta"),
    mean,
    ~ quantile(.x, probs = c(0.05, 0.5, 0.95)),
    .cores = cores
  ) %>%
    mutate(variable = str_remove(str_remove(variable, "theta\\["), "\\]")) %>%
    separate(variable, into = c("dyad_id", "time_index"), sep = ",") %>%
    mutate(
      dyad_id = as.integer(dyad_id),
      time_index = as.integer(time_index)
      #model = model_name,
      #data_set = data_set
    )

  if (save_alpha) {
    alpha_data <- obj$summary(
      variables = c("alpha"),
      mean,
      ~ quantile(.x, probs = c(0.05, 0.5, 0.95)),
      .cores = cores
    )

    saveRDS(
      alpha_data,
      file = paste0(
        'data/alphas/alpha_',
        save_file_name
      )
    )
  }

  obj_data <- obj_data %>%
    left_join(dyad_identifiers, by = c('dyad_id', 'time_index'))

  attr(obj_data, "metadata") <- list(
    temporal_resolution = temporal_resolution,
    dir_type = dir_short,
    model_type = model,
    min_n_events = min_n_events,
    start_year = start_year
  )

  saveRDS(
    obj_data,
    file = paste0(
      'data/theta_trajectories/theta_',
      save_file_name
    )
  )
}

library(tidyverse)

gdelt_to_stan <- function(
  data,
  years = 2000:2025,
  resolution = c('monthly', 'yearly'),
  directed = c(TRUE, FALSE),
  grouping_var = c('TriClass', 'QuadClass', 'PentClass', 'EventRootCode'),
  reference_category = NULL,
  reference_hostile = NULL,
  min_n_events = 1,
  chunk_size = 100
) {
  resolution <- match.arg(resolution, choices = c('monthly', 'yearly'))

  Anum <- length(unique(data[[grouping_var]]))

  if (!is.null(reference_category)) {
    reference_category <- as.character(reference_category)
    if (!(reference_category %in% unique(data[[grouping_var]]))) {
      warning(
        "Reference category not found in the data. Please check your reference_category parameter. The first category in the data will be used as the reference category."
      )
      reference_category <- NULL
    }
  }

  if (!is.null(reference_hostile)) {
    reference_hostile <- as.character(reference_hostile)
    if (!(reference_hostile %in% unique(data[[grouping_var]]))) {
      warning(
        "Reference hostile category not found in the data. Please check your reference_hostile parameter. Running without reference hostile category. Note: the last category in the data will be used as the reference hostile category if reference_hostile is not specified and will be forced in the opposite direction of the reference category if reference_category is specified, or the first category in the data if reference_category is not specified."
      )
      reference_hostile <- NULL
    }
  }

  data <- grouped_events_to_dyad_period(
    data = data,
    resolution = resolution,
    grouping_var = grouping_var,
    directed = directed,
    reference_category = reference_category,
    reference_hostile = reference_hostile
  )

  if (resolution == 'yearly') {
    temp_skeleton <- expand_grid(
      dyad = unique(data$dyad),
      year = years
    )

    data <- data %>%
      right_join(
        temp_skeleton,
        by = c('dyad', 'year')
      ) %>%
      arrange(dyad, year) %>%
      mutate(across(
        c(starts_with("EventClass"), total_events),
        ~ replace_na(.x, 0)
      )) %>%
      mutate(
        is_obs = if_else(total_events == 0, 0, 1),
      ) %>%
      ungroup()

    T <- length(years)
  } else if (resolution == 'monthly') {
    temp_skeleton <- expand_grid(
      dyad = unique(data$dyad),
      year = years,
      month = 1:12
    )

    data <- data %>%
      right_join(
        temp_skeleton,
        by = c('dyad', 'year', 'month')
      ) %>%
      arrange(dyad, year, month) %>%
      mutate(across(
        c(starts_with("EventClass"), total_events),
        ~ replace_na(.x, 0)
      )) %>%
      mutate(
        is_obs = if_else(total_events == 0, 0, 1),
      ) %>%
      ungroup()

    T <- length(years) * 12
  }

  drop_dyads <- data %>%
    group_by(dyad) %>%
    summarise(total_events = sum(total_events), .groups = 'drop') %>%
    filter(total_events < min_n_events) %>%
    pull(dyad)

  data <- data %>%
    filter(!(dyad %in% drop_dyads))

  obs_matrix <- data %>%
    pull(is_obs) %>%
    matrix(nrow = length(unique(data$dyad)), byrow = TRUE)

  events_list <- data %>%
    group_by(dyad) %>%
    group_split() %>%
    map(~ select(.x, starts_with("EventClass")) %>% as.matrix())

  events_array <- array(
    dim = c(
      length(unique(data$dyad)),
      nrow(data) / length(unique(data$dyad)),
      Anum
    )
  )

  for (i in seq_along(events_list)) {
    events_array[i, , ] <- events_list[[i]]
  }

  data_list <- list(
    D = length(unique(data$dyad)),
    T = T,
    A = Anum,
    is_obs = obs_matrix,
    Y = events_array,
    C = chunk_size
  )

  return(data_list)
}


gdelt_to_stan_weight <- function(
  data,
  years = 2000:2025,
  resolution = c('monthly', 'yearly'),
  directed = c(TRUE, FALSE),
  grouping_var = c('TriClass', 'QuadClass', 'PentClass', 'EventRootCode'),
  reference_category = NULL,
  reference_hostile = NULL,
  min_n_events = 1,
  chunk_size = 100
) {
  resolution <- match.arg(resolution, choices = c('monthly', 'yearly'))

  Anum <- length(unique(data[[grouping_var]]))

  if (!is.null(reference_category)) {
    reference_category <- as.character(reference_category)
    if (!(reference_category %in% unique(data[[grouping_var]]))) {
      warning(
        "Reference category not found in the data. Please check your reference_category parameter. The first category in the data will be used as the reference category."
      )
      reference_category <- NULL
    }
  }

  if (!is.null(reference_hostile)) {
    reference_hostile <- as.character(reference_hostile)
    if (!(reference_hostile %in% unique(data[[grouping_var]]))) {
      warning(
        "Reference hostile category not found in the data. Please check your reference_hostile parameter. Running without reference hostile category. Note: the last category in the data will be used as the reference hostile category if reference_hostile is not specified and will be forced in the opposite direction of the reference category if reference_category is specified, or the first category in the data if reference_category is not specified."
      )
      reference_hostile <- NULL
    }
  }

  data <- grouped_events_to_dyad_period(
    data = data,
    resolution = resolution,
    grouping_var = grouping_var,
    directed = directed,
    reference_category = reference_category,
    reference_hostile = reference_hostile
  )

  if (resolution == 'yearly') {
    temp_skeleton <- expand_grid(
      dyad = unique(data$dyad),
      year = years
    )

    data <- data %>%
      right_join(
        temp_skeleton,
        by = c('dyad', 'year')
      ) %>%
      arrange(dyad, year) %>%
      mutate(across(
        c(starts_with("EventClass"), total_events),
        ~ replace_na(.x, 0)
      )) %>%
      mutate(
        is_obs = if_else(total_events == 0, 0, 1),
      ) %>%
      ungroup()

    T <- length(years)
  } else if (resolution == 'monthly') {
    temp_skeleton <- expand_grid(
      dyad = unique(data$dyad),
      year = years,
      month = 1:12
    )

    data <- data %>%
      right_join(
        temp_skeleton,
        by = c('dyad', 'year', 'month')
      ) %>%
      arrange(dyad, year, month) %>%
      mutate(across(
        c(starts_with("EventClass"), total_events),
        ~ replace_na(.x, 0)
      )) %>%
      mutate(
        is_obs = if_else(total_events == 0, 0, 1),
      ) %>%
      ungroup()

    T <- length(years) * 12
  }

  drop_dyads <- data %>%
    group_by(dyad) %>%
    summarise(total_events = sum(total_events), .groups = 'drop') %>%
    filter(total_events < min_n_events) %>%
    pull(dyad)

  data <- data %>%
    filter(!(dyad %in% drop_dyads))

  obs_matrix <- data %>%
    pull(is_obs) %>%
    matrix(nrow = length(unique(data$dyad)), byrow = TRUE)

  events_list <- data %>%
    group_by(dyad) %>%
    group_split() %>%
    map(~ select(.x, starts_with("EventClass")) %>% as.matrix())

  events_array <- array(
    dim = c(
      length(unique(data$dyad)),
      nrow(data) / length(unique(data$dyad)),
      Anum
    )
  )

  for (i in seq_along(events_list)) {
    events_array[i, , ] <- events_list[[i]]
  }

  # -- Dyad-period weights (addresses dyad and period dominance) --------------

  dyad_totals <- apply(events_array, 1, sum)
  period_totals <- apply(events_array, 2, sum)

  dyad_weight <- 1 / log1p(pmax(dyad_totals, 1))
  dyad_weight <- dyad_weight / mean(dyad_weight)

  period_weight <- 1 / log1p(pmax(period_totals, 1))
  period_weight <- period_weight / mean(period_weight)

  # -- Action type weights (upweights rare, theoretically informative classes) -
  action_means <- apply(events_array, 3, function(x) {
    nonzero <- x[x > 0]
    if (length(nonzero) == 0) {
      return(1)
    } # fallback for all-zero action types
    mean(nonzero)
  })
  action_weight <- 1 / log1p(action_means)
  action_weight <- action_weight / mean(action_weight) # normalise to mean 1

  # -- Update data_list --------------------------------------------------------
  data_list <- list(
    D = length(unique(data$dyad)),
    T = T,
    A = Anum,
    is_obs = obs_matrix,
    Y = events_array,
    C = chunk_size,
    dyad_weight = dyad_weight,
    period_weight = period_weight,
    action_weight = action_weight
  )

  return(data_list)
}

make_stan_data <- function(
  data,
  years = 2000:2025,
  resolution = c('monthly', 'yearly'),
  event_type = c('pentaclass', 'quadclass'),
  data_set = c('gdelt', 'icews'),
  min_n_events = 1,
  chunk_size = 100
) {
  data_set <- match.arg(data_set, choices = c('gdelt', 'icews'))
  if (data_set == 'icews' && max(years) > 2023) {
    warning("ICEWS data only goes up to 2023. Adjusting max year to 2023.")
    years <- min(years):2023
  }

  resolution <- match.arg(resolution, choices = c('monthly', 'yearly'))
  event_type <- match.arg(event_type, choices = c('pentaclass', 'quadclass'))

  if (any(str_detect(colnames(data), "PC")) & event_type == 'quadclass') {
    warning(
      "Data contains pentaclass columns but event_type is set to 'quadclass'. Please check your data and event_type parameter. Running with event_type = 'pentaclass'."
    )
    event_type <- 'pentaclass'
  } else if (
    any(str_detect(colnames(data), "QC")) & event_type == 'pentaclass'
  ) {
    warning(
      "Data contains quadclass columns but event_type is set to 'pentaclass'. Please check your data and event_type parameter. Running with event_type = 'quadclass'."
    )
    event_type <- 'quadclass'
  }

  Anum <- ifelse(event_type == 'quadclass', 4, 5)

  if (resolution == 'yearly') {
    temp_skeleton <- expand_grid(
      dyad = unique(data$dyad),
      year = years
    )

    data <- data %>%
      right_join(
        temp_skeleton,
        by = c('dyad', 'year')
      ) %>%
      arrange(dyad, year) %>%
      replace_na(list(
        QC1 = 0,
        QC2 = 0,
        QC3 = 0,
        QC4 = 0,
        PC0 = 0,
        PC1 = 0,
        PC2 = 0,
        PC3 = 0,
        PC4 = 0,
        total_events = 0
      )) %>%
      mutate(
        is_obs = if_else(total_events == 0, 0, 1),
      ) %>%
      ungroup()

    T <- length(years)
  } else if (resolution == 'monthly') {
    temp_skeleton <- expand_grid(
      dyad = unique(data$dyad),
      year = years,
      month = 1:12
    )

    data <- data %>%
      right_join(
        temp_skeleton,
        by = c('dyad', 'year', 'month')
      ) %>%
      arrange(dyad, year, month) %>%
      replace_na(list(
        QC1 = 0,
        QC2 = 0,
        QC3 = 0,
        QC4 = 0,
        PC0 = 0,
        PC1 = 0,
        PC2 = 0,
        PC3 = 0,
        PC4 = 0,
        total_events = 0
      )) %>%
      mutate(
        is_obs = if_else(total_events == 0, 0, 1),
      ) %>%
      ungroup()

    T <- length(years) * 12
  }

  drop_dyads <- data %>%
    group_by(dyad) %>%
    summarise(total_events = sum(total_events), .groups = 'drop') %>%
    filter(total_events < min_n_events) %>%
    pull(dyad)

  data <- data %>%
    filter(!(dyad %in% drop_dyads))

  obs_matrix <- data %>%
    pull(is_obs) %>%
    matrix(nrow = length(unique(data$dyad)), byrow = TRUE)

  events_list <- data %>%
    group_by(dyad) %>%
    group_split() %>%
    map(~ select(.x, (starts_with("QC") | starts_with("PC"))) %>% as.matrix())

  events_array <- array(
    dim = c(
      length(unique(data$dyad)),
      nrow(data) / length(unique(data$dyad)),
      Anum
    )
  )

  for (i in seq_along(events_list)) {
    events_array[i, , ] <- events_list[[i]]
  }

  data_list <- list(
    D = length(unique(data$dyad)),
    T = T,
    A = Anum,
    is_obs = obs_matrix,
    Y = events_array,
    C = chunk_size
  )

  return(data_list)
}


make_stan_data_wstartvals <- function(
  data,
  years = 2000:2025,
  resolution = c('monthly', 'yearly'),
  event_type = c('pentaclass', 'quadclass'),
  data_set = c('gdelt', 'icews'),
  min_n_events = 1,
  chunk_size = 100,
  starting_values = NULL
) {
  data_set <- match.arg(data_set, choices = c('gdelt', 'icews'))
  if (data_set == 'icews' && max(years) > 2023) {
    warning("ICEWS data only goes up to 2023. Adjusting max year to 2023.")
    years <- min(years):2023
  }

  resolution <- match.arg(resolution, choices = c('monthly', 'yearly'))
  event_type <- match.arg(event_type, choices = c('pentaclass', 'quadclass'))

  if (any(str_detect(colnames(data), "PC")) & event_type == 'quadclass') {
    warning(
      "Data contains pentaclass columns but event_type is set to 'quadclass'. Please check your data and event_type parameter. Running with event_type = 'pentaclass'."
    )
    event_type <- 'pentaclass'
  } else if (
    any(str_detect(colnames(data), "QC")) & event_type == 'pentaclass'
  ) {
    warning(
      "Data contains quadclass columns but event_type is set to 'pentaclass'. Please check your data and event_type parameter. Running with event_type = 'quadclass'."
    )
    event_type <- 'quadclass'
  }

  Anum <- ifelse(event_type == 'quadclass', 4, 5)

  if (resolution == 'yearly') {
    temp_skeleton <- expand_grid(
      dyad = unique(data$dyad),
      year = years
    )

    data <- data %>%
      right_join(
        temp_skeleton,
        by = c('dyad', 'year')
      ) %>%
      arrange(dyad, year) %>%
      replace_na(list(
        QC1 = 0,
        QC2 = 0,
        QC3 = 0,
        QC4 = 0,
        PC0 = 0,
        PC1 = 0,
        PC2 = 0,
        PC3 = 0,
        PC4 = 0,
        total_events = 0
      )) %>%
      mutate(
        is_obs = if_else(total_events == 0, 0, 1),
      ) %>%
      ungroup()

    T <- length(years)
  } else if (resolution == 'monthly') {
    temp_skeleton <- expand_grid(
      dyad = unique(data$dyad),
      year = years,
      month = 1:12
    )

    data <- data %>%
      right_join(
        temp_skeleton,
        by = c('dyad', 'year', 'month')
      ) %>%
      arrange(dyad, year, month) %>%
      replace_na(list(
        QC1 = 0,
        QC2 = 0,
        QC3 = 0,
        QC4 = 0,
        PC0 = 0,
        PC1 = 0,
        PC2 = 0,
        PC3 = 0,
        PC4 = 0,
        total_events = 0
      )) %>%
      mutate(
        is_obs = if_else(total_events == 0, 0, 1),
      ) %>%
      ungroup()

    T <- length(years) * 12
  }

  drop_dyads <- data %>%
    group_by(dyad) %>%
    summarise(total_events = sum(total_events), .groups = 'drop') %>%
    filter(total_events < min_n_events) %>%
    pull(dyad)

  data <- data %>%
    filter(!(dyad %in% drop_dyads))

  obs_matrix <- data %>%
    pull(is_obs) %>%
    matrix(nrow = length(unique(data$dyad)), byrow = TRUE)

  events_list <- data %>%
    group_by(dyad) %>%
    group_split() %>%
    map(~ select(.x, (starts_with("QC") | starts_with("PC"))) %>% as.matrix())

  events_array <- array(
    dim = c(
      length(unique(data$dyad)),
      nrow(data) / length(unique(data$dyad)),
      Anum
    )
  )

  for (i in seq_along(events_list)) {
    events_array[i, , ] <- events_list[[i]]
  }

  if (is.null(starting_values)) {
    starting_values <- list(
      theta0_mu = rep(0, length(unique(data$dyad))),
      theta0_sigma = rep(1, length(unique(data$dyad)))
    )
  } else if (is.numeric(starting_values) && length(starting_values) == 2) {
    starting_values <- list(
      theta0_mu = rep(starting_values[1], length(unique(data$dyad))),
      theta0_sigma = rep(starting_values[2], length(unique(data$dyad)))
    )
  } else if (is.data.frame(starting_values) && ncol(starting_values) == 3) {
    start_vals <- data %>%
      select(dyad) %>%
      left_join(starting_values, by = 'dyad') %>%
      group_by(dyad) %>%
      slice(1) %>%
      replace_na(list(
        theta0_mu = 0,
        theta0_sigma = 1
      )) %>%
      ungroup()
    starting_values <- list(
      theta0_mu = start_vals$theta0_mu,
      theta0_sigma = start_vals$theta0_sigma
    )
  } else if (is.character(starting_values) && starting_values == 'oldver') {
    # This is a special case to replicate the starting values from the old version of the model for comparison purposes. It will be ignored in the make_stan_data function and handled separately in the model fitting function.
    starting_values <- 'oldver'
  } else {
    stop(
      "Invalid format for starting_values. Please provide NULL, a numeric vector of length 2, or a data frame with columns 'dyad', 'theta0_mu', and 'theta0_sigma'."
    )
  }

  if (starting_values == 'oldver') {
    data_list <- list(
      D = length(unique(data$dyad)),
      T = T,
      A = Anum,
      is_obs = obs_matrix,
      Y = events_array,
      C = chunk_size
    )
  } else {
    data_list <- list(
      D = length(unique(data$dyad)),
      T = T,
      A = Anum,
      is_obs = obs_matrix,
      Y = events_array,
      C = chunk_size,
      theta0_mu = starting_values$theta0_mu,
      theta0_sigma = starting_values$theta0_sigma
    )
  }

  return(data_list)
}


make_dyad_ids <- function(
  data,
  years = 2000:2025,
  resolution = c('monthly', 'yearly'),
  data_set = c('gdelt', 'icews'),
  min_n_events = 1,
  dir_type = c('directed', 'undirected')
) {
  if (data_set == 'icews' && max(years) > 2023) {
    warning("ICEWS data only goes up to 2023. Adjusting max year to 2023.")
    years <- min(years):2023
  }

  resolution <- match.arg(resolution, choices = c('monthly', 'yearly'))

  dyad_secondid <- data %>%
    group_by(dyad) %>%
    slice(1) %>%
    mutate(
      a1 = str_sub(dyad, 1, 3),
      a2 = str_sub(dyad, 5, 7),
      dyad2 = case_when(
        a1 < a2 ~ paste(a1, a2, sep = "_"),
        TRUE ~ paste(a2, a1, sep = "_")
      )
    ) %>%
    ungroup() %>%
    select(dyad, dyad2)

  if (resolution == 'yearly') {
    temp_skeleton <- expand_grid(
      dyad = unique(data$dyad),
      year = years
    )

    data <- data %>%
      right_join(
        temp_skeleton,
        by = c('dyad', 'year')
      ) %>%
      arrange(dyad, year)
  } else if (resolution == 'monthly') {
    temp_skeleton <- expand_grid(
      dyad = unique(data$dyad),
      year = years,
      month = 1:12
    )

    data <- data %>%
      right_join(
        temp_skeleton,
        by = c('dyad', 'year', 'month')
      ) %>%
      arrange(dyad, year, month)
  }

  data <- data %>%
    replace_na(list(
      total_events = 0
    )) %>%
    ungroup()

  drop_dyads <- data %>%
    group_by(dyad) %>%
    summarise(total_events = sum(total_events), .groups = 'drop') %>%
    filter(total_events < min_n_events) %>%
    pull(dyad)

  data <- data %>%
    filter(!(dyad %in% drop_dyads)) %>%
    left_join(dyad_secondid, by = 'dyad')

  data <- data %>%
    group_by(dyad) %>%
    mutate(dyad_id = cur_group_id(), time_index = row_number()) %>%
    ungroup() %>%
    select(any_of(c('dyad_id', 'time_index', 'dyad', 'dyad2', 'year', 'month')))

  return(data)
}

estimate_starting_values <- function(
  data,
  calib_years = 1985:1989,
  anchor_in = c('goldstein', 'tone')
) {
  if (length(calib_years) == 1) {
    stop(
      "calib_years should be a range of years, not a single year, for random effects to be estimated"
    )
  }

  anchor_in <- match.arg(anchor_in, choices = c('goldstein', 'tone'))

  data <- data %>%
    filter(year %in% calib_years)

  if (anchor_in == 'goldstein') {
    data <- data %>%
      mutate(dv = mean_goldstein)
  } else if (anchor_in == 'tone') {
    data <- data %>%
      mutate(dv = mean_tone)
  }

  m <- lme4::lmer(
    dv ~ 1 + (1 | dyad),
    data = data,
    weights = log(total_events + 1)
  )
  random_effects <- lme4::ranef(m) %>%
    as.data.frame() %>%
    select(grp, condval, condsd) %>%
    rename(
      dyad = grp,
      theta0_mu = condval,
      theta0_sigma = condsd
    )
  return(random_effects)
}

library(cmdstanr)


model_flex_runner <- function(
  model_path,
  data,
  resolution = c('monthly', 'yearly'),
  grouping_var = c('TriClass', 'QuadClass', 'PentClass', 'EventRootCode'),
  reference_category = NULL,
  reference_hostile = NULL,
  directed = c(TRUE, FALSE),
  years = 2000:2025,
  output_dir,
  chains = 1,
  parallel_chains = 1,
  threads_per_chain = 1,
  iter_warmup = 500,
  iter_sampling = 1000,
  min_n_events = 1,
  chunk_size = 100,
  seed = 123
) {
  stanmod <- cmdstan_model(model_path, cpp_options = list(stan_threads = TRUE))

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  resolution <- match.arg(resolution, choices = c('monthly', 'yearly'))

  if (str_detect(model_path, 'weight')) {
    standata <- gdelt_to_stan_weight(
      data,
      years = years,
      grouping_var = grouping_var,
      reference_category = reference_category,
      reference_hostile = reference_hostile,
      directed = directed,
      resolution = resolution,
      min_n_events = min_n_events,
      chunk_size = chunk_size
    )
  } else {
    standata <- gdelt_to_stan(
      data,
      years = years,
      grouping_var = grouping_var,
      reference_category = reference_category,
      reference_hostile = reference_hostile,
      directed = directed,
      resolution = resolution,
      min_n_events = min_n_events,
      chunk_size = chunk_size
    )
  }

  init_fn <- function() {
    list(
      mu_theta0 = 0,
      sigma_theta0 = 0.5,
      z_theta0 = rep(0, standata$D),
      theta_raw = matrix(0, standata$D, standata$T),
      mu_log_noise = log(0.2),
      sigma_log_noise = 0.3,
      mu_log_phi = 0,
      sigma_log_phi = 0.5,
      phi = rep(1, standata$D),
      process_noise = rep(0.2, standata$D),
      alpha_raw = rep(0, standata$A - 2),
      alpha_hostile = 0.5,
      mu_intercept_raw = rep(0, standata$A - 1)
    )
  }

  obj <- stanmod$sample(
    data = standata,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    threads_per_chain = threads_per_chain,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    output_dir = output_dir,
    init = init_fn,
    refresh = 1,
    show_messages = TRUE,
    show_exceptions = TRUE
  )

  saveRDS(
    obj,
    file.path(
      output_dir,
      paste0(
        basename(output_dir),
        '.rds'
      )
    )
  )
}


# model_comparison_runner <- function(
#   model_path,
#   data,
#   data_set = c('gdelt', 'icews'),
#   actor1,
#   actor2,
#   years = 2000:2025,
#   output_dir,
#   chains = 1,
#   parallel_chains = 1,
#   iter_warmup = 500,
#   iter_sampling = 1000
# ) {
#   stanmod <- cmdstan_model(model_path)

#   obj <- stanmod$sample(
#     data = make_stan_data_irt(
#       data,
#       years = years,
#       actors = c(actor1, actor2),
#       type = basename(model_path) %>%
#         str_remove('_irt.stan'),
#       data_set = data_set
#     ),
#     seed = 123,
#     chains = chains,
#     parallel_chains = parallel_chains,
#     iter_warmup = iter_warmup,
#     iter_sampling = iter_sampling,
#     output_dir = output_dir
#   )

#   saveRDS(
#     obj,
#     file.path(
#       output_dir,
#       paste0(
#         basename(model_path) %>% tools::file_path_sans_ext(),
#         '_',
#         actor1,
#         '_',
#         actor2,
#         '_',
#         data_set,
#         '.rds'
#       )
#     )
#   )
# }

model_hierarchial_runner <- function(
  model_path,
  data,
  data_set = c('gdelt', 'icews'),
  resolution = c('monthly', 'yearly'),
  years = 2000:2025,
  output_dir,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 1,
  iter_warmup = 500,
  iter_sampling = 1000,
  min_n_events = 1,
  chunk_size = 100,
  #estimate_starting_values = FALSE,
  calib_years = 5
) {
  stanmod <- cmdstan_model(model_path, cpp_options = list(stan_threads = TRUE))

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  resolution <- match.arg(resolution, choices = c('monthly', 'yearly'))

  # if (estimate_starting_values == TRUE) {
  #   estimated_starting_values <- estimate_starting_values(
  #     data,
  #     calib_years = ((min(years) - calib_years):(min(years) - 1)),
  #     anchor_in = 'goldstein'
  #   )
  # } else if (estimate_starting_values == 'oldver') {
  #   estimated_starting_values <- 'oldver'
  # } else {
  #   estimated_starting_values <- NULL
  # }

  obj <- stanmod$sample(
    data = make_stan_data(
      data,
      years = years,
      data_set = data_set,
      resolution = resolution,
      min_n_events = min_n_events,
      chunk_size = chunk_size #,
      #starting_values = estimated_starting_values
    ),
    seed = 123,
    chains = chains,
    parallel_chains = parallel_chains,
    threads_per_chain = threads_per_chain,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    output_dir = output_dir,
    refresh = 10
  )

  saveRDS(
    obj,
    file.path(
      output_dir,
      paste0(
        basename(model_path) %>% tools::file_path_sans_ext(),
        data_set,
        '.rds'
      )
    )
  )
}

library(foreach)
library(tidyverse)

cameo_to_quad <- function(cameo_code) {
  qc <- icews::cameo_codes
  qc <- qc %>%
    select(cameo_code, quad_category) %>%
    mutate(
      quad_class = case_when(
        quad_category == "verbal cooperation" ~ 1,
        quad_category == "material cooperation" ~ 2,
        quad_category == "verbal conflict" ~ 3,
        quad_category == "material conflict" ~ 4,
        TRUE ~ NA_real_
      )
    )
  return(qc$quad_class[match(cameo_code, qc$cameo_code)])
}


relevant_actors_icews <- c(
  "Air Force",
  "Army",
  "Cabinet",
  "Coast Guard",
  "Ministry",
  "Security Mi",
  "Execut",
  "Gov",
  "Legisl",
  "Lower House",
  "Milit",
  "Supreme Court",
  "Navy",
  "Police",
  "Unicameral",
  "Upper House"
)


file <- "data/icews/20230108-icews-events.zip"
files <- list.files('data/icews', pattern = 'zip', full.names = TRUE)


all_icews <- foreach(file = files) %do%
  {
    read_delim(
      unz(
        file,
        unzip(file, list = TRUE)$Name[1]
      ),
      num_threads = 1,
      progress = FALSE
    )
  }

all_icews2 <- foreach(i = 1:length(all_icews), .final = bind_rows) %do%
  {
    all_icews[[i]] %>%
      mutate(
        cameo = as.character(`CAMEO Code`),
        cameo = case_when(
          str_length(cameo) < 3 ~ paste0("0", cameo),
          TRUE ~ cameo
        )
      ) %>%
      select(-`CAMEO Code`)
  }

all_icews2 <- all_icews2 %>%
  mutate(
    actor1_list = str_split(`Source Sectors`, ","),
    actor2_list = str_split(`Target Sectors`, ",")
  ) %>%
  filter(`Source Country` != `Target Country`)

unique_actors <- unique(c(
  unlist(all_icews2$actor1_list),
  unlist(all_icews2$actor2_list)
))


all_icews3 <- all_icews2 %>%
  filter(
    str_detect(`Source Sectors`, paste(relevant_actors_icews, collapse = "|")) &
      str_detect(`Target Sectors`, paste(relevant_actors_icews, collapse = "|"))
  )

qc <- icews::cameo_codes
qc <- qc %>%
  select(cameo_code, quad_category) %>%
  mutate(
    quad_class = case_when(
      quad_category == "verbal cooperation" ~ 1,
      quad_category == "material cooperation" ~ 2,
      quad_category == "verbal conflict" ~ 3,
      quad_category == "material conflict" ~ 4,
      TRUE ~ NA_real_
    )
  )

all_icews3 <- all_icews3 %>%
  left_join(qc, by = c("cameo" = "cameo_code"))

all_icews3 <- all_icews3 %>%
  mutate(year = year(`Event Date`), month = month(`Event Date`))

all_icews4 <- all_icews3 %>%
  mutate(
    cameo = case_when(is.na(quad_class) ~ paste0("0", cameo), TRUE ~ cameo)
  ) %>%
  select(-quad_category, -quad_class) %>%
  left_join(qc, by = c("cameo" = "cameo_code"))


icews_country_month_dyads <- all_icews4 %>%
  group_by(`Source Country`, `Target Country`, year, month) %>%
  summarise(
    sQC1 = sum(quad_class == 1),
    sQC2 = sum(quad_class == 2),
    sQC3 = sum(quad_class == 3),
    sQC4 = sum(quad_class == 4),
    mean_intensity = mean(Intensity, na.rm = TRUE),
    mean_positive = mean(Intensity > 0, na.rm = TRUE),
    mean_negative = mean(Intensity < 0, na.rm = TRUE),
    total_events = n(),
    pQC1 = sQC1 / total_events,
    pQC2 = sQC2 / total_events,
    pQC3 = sQC3 / total_events,
    pQC4 = sQC4 / total_events,

    .groups = "drop"
  ) %>%
  rename(source_country = `Source Country`, target_country = `Target Country`)

saveRDS(icews_country_month_dyads, "data/icews_country_month_dyads.rds")

library(dplyr)
library(stringr)

# ── 1. Build the lookup table from your data ──────────────────────────────────
cameo_df <- tibble::tribble(
  ~CAMEOEVENTCODE , ~GoldsteinScore ,
  "01"            ,   0.0           , "010"  ,   0.0 , "011"  ,  -0.1 , "012"  ,  -0.4 , "013"  ,   0.4 ,
  "014"           ,   0.0           , "015"  ,   0.0 , "016"  ,  -2.0 , "017"  ,   0.0 , "018"  ,   3.4 ,
  "019"           ,   3.4           , "02"   ,   3.0 , "020"  ,   3.0 , "021"  ,   3.4 , "0211" ,   3.4 ,
  "0212"          ,   3.4           , "0213" ,   3.4 , "0214" ,   3.4 , "022"  ,   3.2 , "023"  ,   3.4 ,
  "0231"          ,   3.4           , "0232" ,   3.4 , "0233" ,   3.4 , "0234" ,   3.4 , "024"  ,  -0.3 ,
  "0241"          ,  -0.3           , "0242" ,  -0.3 , "0243" ,  -0.3 , "0244" ,  -0.3 , "025"  ,  -0.3 ,
  "0251"          ,  -0.3           , "0252" ,  -0.3 , "0253" ,  -0.3 , "0254" ,  -0.3 , "0255" ,  -0.3 ,
  "0256"          ,  -0.3           , "026"  ,   4.0 , "027"  ,   4.0 , "028"  ,   4.0 , "03"   ,   4.0 ,
  "030"           ,   4.0           , "031"  ,   5.2 , "0311" ,   5.2 , "0312" ,   5.2 , "0313" ,   5.2 ,
  "0314"          ,   5.2           , "032"  ,   4.5 , "033"  ,   5.2 , "0331" ,   5.2 , "0332" ,   5.2 ,
  "0333"          ,   5.2           , "0334" ,   6.0 , "034"  ,   7.0 , "0341" ,   7.0 , "0342" ,   7.0 ,
  "0343"          ,   7.0           , "0344" ,   7.0 , "035"  ,   7.0 , "0351" ,   7.0 , "0352" ,   7.0 ,
  "0353"          ,   7.0           , "0354" ,   7.0 , "0355" ,   7.0 , "0356" ,   7.0 , "036"  ,   4.0 ,
  "037"           ,   5.0           , "038"  ,   7.0 , "039"  ,   5.0 , "04"   ,   1.0 , "040"  ,   1.0 ,
  "041"           ,   1.0           , "042"  ,   1.9 , "043"  ,   2.8 , "044"  ,   2.5 , "045"  ,   5.0 ,
  "046"           ,   7.0           , "05"   ,   3.5 , "050"  ,   3.5 , "051"  ,   3.4 , "052"  ,   3.5 ,
  "053"           ,   3.8           , "054"  ,   6.0 , "055"  ,   7.0 , "056"  ,   7.0 , "057"  ,   8.0 ,
  "06"            ,   6.0           , "060"  ,   6.0 , "061"  ,   6.4 , "062"  ,   7.4 , "063"  ,   7.4 ,
  "064"           ,   7.0           , "07"   ,   7.0 , "070"  ,   7.0 , "071"  ,   7.4 , "072"  ,   8.3 ,
  "073"           ,   7.4           , "074"  ,   8.5 , "075"  ,   7.0 , "08"   ,   5.0 , "080"  ,   5.0 ,
  "081"           ,   5.0           , "0811" ,   5.0 , "0812" ,   5.0 , "0813" ,   5.0 , "0814" ,   5.0 ,
  "082"           ,   5.0           , "083"  ,   5.0 , "0831" ,   5.0 , "0832" ,   5.0 , "0833" ,   5.0 ,
  "0834"          ,   5.0           , "084"  ,   7.0 , "0841" ,   7.0 , "0842" ,   7.0 , "085"  ,   7.0 ,
  "086"           ,   9.0           , "0861" ,   9.0 , "0862" ,   9.0 , "0863" ,   9.0 , "087"  ,   9.0 ,
  "0871"          ,   9.0           , "0872" ,   9.0 , "0873" ,   9.0 , "0874" ,  10.0 , "09"   ,  -2.0 ,
  "090"           ,  -2.0           , "091"  ,  -2.0 , "092"  ,  -2.0 , "093"  ,  -2.0 , "094"  ,  -2.0 ,
  "10"            ,  -5.0           , "100"  ,  -5.0 , "101"  ,  -5.0 , "1011" ,  -5.0 , "1012" ,  -5.0 ,
  "1013"          ,  -5.0           , "1014" ,  -5.0 , "102"  ,  -5.0 , "103"  ,  -5.0 , "1031" ,  -5.0 ,
  "1032"          ,  -5.0           , "1033" ,  -5.0 , "1034" ,  -5.0 , "104"  ,  -5.0 , "1041" ,  -5.0 ,
  "1042"          ,  -5.0           , "1043" ,  -5.0 , "1044" ,  -5.0 , "105"  ,  -5.0 , "1051" ,  -5.0 ,
  "1052"          ,  -5.0           , "1053" ,  -5.0 , "1054" ,  -5.0 , "1055" ,  -5.0 , "1056" ,  -5.0 ,
  "107"           ,  -5.0           , "108"  ,  -5.0 , "11"   ,  -2.0 , "110"  ,  -2.0 , "111"  ,  -2.0 ,
  "112"           ,  -2.0           , "1121" ,  -2.0 , "1122" ,  -2.0 , "1123" ,  -2.0 , "1124" ,  -2.0 ,
  "1125"          ,  -2.0           , "113"  ,  -2.0 , "114"  ,  -2.0 , "115"  ,  -2.0 , "116"  ,  -2.0 ,
  "12"            ,  -4.0           , "120"  ,  -4.0 , "121"  ,  -4.0 , "1211" ,  -4.0 , "1212" ,  -4.0 ,
  "122"           ,  -4.0           , "1221" ,  -4.0 , "1222" ,  -4.0 , "1223" ,  -4.0 , "1224" ,  -4.0 ,
  "123"           ,  -4.0           , "1231" ,  -4.0 , "1232" ,  -4.0 , "1233" ,  -4.0 , "1234" ,  -4.0 ,
  "124"           ,  -4.0           , "1241" ,  -4.0 , "1242" ,  -4.0 , "1243" ,  -4.0 , "1244" ,  -4.0 ,
  "1245"          ,  -4.0           , "1246" ,  -4.0 , "125"  ,  -5.0 , "126"  ,  -5.0 , "127"  ,  -5.0 ,
  "128"           ,  -5.0           , "129"  ,  -5.0 , "13"   ,  -6.0 , "130"  ,  -4.4 , "131"  ,  -5.8 ,
  "1311"          ,  -5.8           , "1312" ,  -5.8 , "1313" ,  -5.8 , "132"  ,  -5.8 , "1321" ,  -5.8 ,
  "1322"          ,  -5.8           , "1323" ,  -5.8 , "1324" ,  -5.8 , "133"  ,  -5.8 , "134"  ,  -5.8 ,
  "135"           ,  -5.8           , "136"  ,  -7.0 , "137"  ,  -7.0 , "138"  ,  -7.0 , "1381" ,  -7.0 ,
  "1382"          ,  -7.0           , "1383" ,  -7.0 , "1384" ,  -7.0 , "1385" ,  -7.0 , "139"  ,  -7.0 ,
  "14"            ,  -6.5           , "140"  ,  -6.5 , "141"  ,  -6.5 , "1411" ,  -6.5 , "1412" ,  -6.5 ,
  "1413"          ,  -6.5           , "1414" ,  -6.5 , "142"  ,  -6.5 , "1421" ,  -6.5 , "1422" ,  -6.5 ,
  "1423"          ,  -6.5           , "1424" ,  -6.5 , "143"  ,  -6.5 , "1431" ,  -6.5 , "1432" ,  -6.5 ,
  "1433"          ,  -6.5           , "1434" ,  -6.5 , "144"  ,  -7.5 , "1441" ,  -7.5 , "1442" ,  -7.5 ,
  "1443"          ,  -7.5           , "1444" ,  -7.5 , "145"  ,  -7.5 , "1451" ,  -7.5 , "1452" ,  -7.5 ,
  "1453"          ,  -7.5           , "1454" ,  -7.5 , "15"   ,  -7.2 , "150"  ,  -7.2 , "151"  ,  -7.2 ,
  "152"           ,  -7.2           , "153"  ,  -7.2 , "154"  ,  -7.2 , "16"   ,  -4.0 , "160"  ,  -4.0 ,
  "161"           ,  -4.0           , "162"  ,  -5.6 , "1621" ,  -5.6 , "1622" ,  -5.6 , "1623" ,  -5.6 ,
  "163"           ,  -8.0           , "164"  ,  -7.0 , "165"  ,  -6.5 , "166"  ,  -7.0 , "1661" ,  -7.0 ,
  "1662"          ,  -7.0           , "1663" ,  -7.0 , "17"   ,  -7.0 , "170"  ,  -7.0 , "171"  ,  -9.2 ,
  "1711"          ,  -9.2           , "1712" ,  -9.2 , "172"  ,  -5.0 , "1721" ,  -5.0 , "1722" ,  -5.0 ,
  "1723"          ,  -5.0           , "1724" ,  -5.0 , "173"  ,  -5.0 , "174"  ,  -5.0 , "175"  ,  -9.0 ,
  "18"            ,  -9.0           , "180"  ,  -9.0 , "181"  ,  -9.0 , "182"  ,  -9.5 , "1821" ,  -9.0 ,
  "1822"          ,  -9.0           , "1823" , -10.0 , "183"  , -10.0 , "1831" , -10.0 , "1832" , -10.0 ,
  "1833"          , -10.0           , "184"  ,  -8.0 , "185"  ,  -8.0 , "186"  , -10.0 , "19"   , -10.0 ,
  "190"           , -10.0           , "191"  ,  -9.5 , "192"  ,  -9.5 , "193"  , -10.0 , "194"  , -10.0 ,
  "195"           , -10.0           , "196"  ,  -9.5 , "20"   , -10.0 , "200"  , -10.0 , "201"  ,  -9.5 ,
  "202"           , -10.0           , "203"  , -10.0 , "204"  , -10.0 , "2041" , -10.0 , "2042" , -10.0
)

# ── 2. Helper: extract the root 2-digit code from any CAMEO code ──────────────
get_root <- function(code) {
  # Codes starting with a leading zero (e.g. "011", "0211") → first 2 chars
  # Codes starting with 1x/2x (e.g. "100", "1121", "20") → first 2 chars
  str_sub(str_pad(as.character(code), width = 2, pad = "0"), 1, 2)
}

# ── 3. Classification functions ───────────────────────────────────────────────
assign_quad <- function(root) {
  root_num <- as.integer(root)
  case_when(
    root_num >= 1 & root_num <= 5 ~ 1L,
    root_num >= 6 & root_num <= 9 ~ 2L,
    root_num >= 10 & root_num <= 14 ~ 3L,
    root_num >= 15 & root_num <= 20 ~ 4L,
    TRUE ~ NA_integer_
  )
}

assign_penta <- function(root, quad) {
  root_num <- as.integer(root)
  case_when(
    root_num %in% c(1, 2) ~ 0L, # 01–02 → Verbal cooperation
    root_num == 14 ~ 4L, # 14 → Protest (override quad=3)
    root_num == 16 ~ 3L, # 16 → Reduce relations (override quad=4)
    TRUE ~ quad # all others inherit QuadClass
  )
}

# ── 4. Add all columns ────────────────────────────────────────────────────────
cameo_df <- cameo_df |>
  mutate(
    root = get_root(CAMEOEVENTCODE),
    QuadClass = assign_quad(root),
    PentaClass = assign_penta(root, QuadClass),
    PentaClass_modified = if_else(
      GoldsteinScore <= 1 & QuadClass == 1,
      0L,
      PentaClass
    )
  ) |>
  select(-root) # drop helper column if not needed


cameo_df <- cameo_df |>
  mutate(
    CAMEOLabel = case_when(
      CAMEOEVENTCODE == "01" ~ "Make Public Statement",
      CAMEOEVENTCODE == "010" ~ "Make statement, not specified below",
      CAMEOEVENTCODE == "011" ~ "Decline comment",
      CAMEOEVENTCODE == "012" ~ "Make pessimistic comment",
      CAMEOEVENTCODE == "013" ~ "Make optimistic comment",
      CAMEOEVENTCODE == "014" ~ "Consider policy option",
      CAMEOEVENTCODE == "015" ~ "Acknowledge or claim responsibility",
      CAMEOEVENTCODE == "016" ~ "Deny responsibility",
      CAMEOEVENTCODE == "017" ~ "Engage in symbolic act",
      CAMEOEVENTCODE == "018" ~ "Make empathetic comment",
      CAMEOEVENTCODE == "019" ~ "Express accord",
      CAMEOEVENTCODE == "02" ~ "Appeal",
      CAMEOEVENTCODE ==
        "020" ~ "Make an appeal or request, not specified below",
      CAMEOEVENTCODE ==
        "021" ~ "Appeal for material cooperation, not specified below",
      CAMEOEVENTCODE == "0211" ~ "Appeal for economic cooperation",
      CAMEOEVENTCODE == "0212" ~ "Appeal for military cooperation",
      CAMEOEVENTCODE == "0213" ~ "Appeal for judicial cooperation",
      CAMEOEVENTCODE == "0214" ~ "Appeal for intelligence",
      CAMEOEVENTCODE == "022" ~ "Appeal for diplomatic cooperation",
      CAMEOEVENTCODE == "023" ~ "Appeal for aid, not specified below",
      CAMEOEVENTCODE == "0231" ~ "Appeal for economic aid",
      CAMEOEVENTCODE == "0232" ~ "Appeal for military aid",
      CAMEOEVENTCODE == "0233" ~ "Appeal for humanitarian aid",
      CAMEOEVENTCODE ==
        "0234" ~ "Appeal for military protection or peacekeeping",
      CAMEOEVENTCODE ==
        "024" ~ "Appeal for political reform, not specified below",
      CAMEOEVENTCODE == "0241" ~ "Appeal for change in leadership",
      CAMEOEVENTCODE == "0242" ~ "Appeal for policy change",
      CAMEOEVENTCODE == "0243" ~ "Appeal for rights",
      CAMEOEVENTCODE == "0244" ~ "Appeal for change in institutions, regime",
      CAMEOEVENTCODE == "025" ~ "Appeal to yield, not specified below",
      CAMEOEVENTCODE ==
        "0251" ~ "Appeal for easing of administrative sanctions",
      CAMEOEVENTCODE == "0252" ~ "Appeal for easing of political dissent",
      CAMEOEVENTCODE == "0253" ~ "Appeal for release of persons or property",
      CAMEOEVENTCODE ==
        "0254" ~ "Appeal for easing of economic sanctions, boycott, or embargo",
      CAMEOEVENTCODE ==
        "0255" ~ "Appeal for target to allow international involvement (non-mediation)",
      CAMEOEVENTCODE ==
        "0256" ~ "Appeal for de-escalation of military engagement",
      CAMEOEVENTCODE == "026" ~ "Appeal to others to meet or negotiate",
      CAMEOEVENTCODE == "027" ~ "Appeal to others to settle dispute",
      CAMEOEVENTCODE == "028" ~ "Appeal to engage in or accept mediation",
      CAMEOEVENTCODE == "03" ~ "Express Intent to Cooperate",
      CAMEOEVENTCODE ==
        "030" ~ "Express intent to cooperate, not specified below",
      CAMEOEVENTCODE ==
        "031" ~ "Express intent to engage in material cooperation, not specified below",
      CAMEOEVENTCODE == "0311" ~ "Express intent to cooperate economically",
      CAMEOEVENTCODE == "0312" ~ "Express intent to cooperate militarily",
      CAMEOEVENTCODE ==
        "0313" ~ "Express intent to cooperate on judicial matters",
      CAMEOEVENTCODE == "0314" ~ "Express intent to cooperate on intelligence",
      CAMEOEVENTCODE ==
        "032" ~ "Express intent to engage in diplomatic cooperation",
      CAMEOEVENTCODE ==
        "033" ~ "Express intent to provide material aid, not specified below",
      CAMEOEVENTCODE == "0331" ~ "Express intent to provide economic aid",
      CAMEOEVENTCODE == "0332" ~ "Express intent to provide military aid",
      CAMEOEVENTCODE == "0333" ~ "Express intent to provide humanitarian aid",
      CAMEOEVENTCODE ==
        "0334" ~ "Express intent to provide military protection or peacekeeping",
      CAMEOEVENTCODE ==
        "034" ~ "Express intent to institute political reform, not specified below",
      CAMEOEVENTCODE == "0341" ~ "Express intent to change leadership",
      CAMEOEVENTCODE == "0342" ~ "Express intent to change policy",
      CAMEOEVENTCODE == "0343" ~ "Express intent to provide rights",
      CAMEOEVENTCODE ==
        "0344" ~ "Express intent to change institutions, regime",
      CAMEOEVENTCODE == "035" ~ "Express intent to yield, not specified below",
      CAMEOEVENTCODE ==
        "0351" ~ "Express intent to ease administrative sanctions",
      CAMEOEVENTCODE == "0352" ~ "Express intent to ease popular dissent",
      CAMEOEVENTCODE ==
        "0353" ~ "Express intent to release persons or property",
      CAMEOEVENTCODE ==
        "0354" ~ "Express intent to ease economic sanctions, boycott, or embargo",
      CAMEOEVENTCODE ==
        "0355" ~ "Express intent to allow international involvement (non-mediation)",
      CAMEOEVENTCODE ==
        "0356" ~ "Express intent to de-escalate military engagement",
      CAMEOEVENTCODE == "036" ~ "Express intent to meet or negotiate",
      CAMEOEVENTCODE == "037" ~ "Express intent to settle dispute",
      CAMEOEVENTCODE == "038" ~ "Express intent to accept mediation",
      CAMEOEVENTCODE == "039" ~ "Express intent to mediate",
      CAMEOEVENTCODE == "04" ~ "Consult",
      CAMEOEVENTCODE == "040" ~ "Consult, not specified below",
      CAMEOEVENTCODE == "041" ~ "Discuss by telephone",
      CAMEOEVENTCODE == "042" ~ "Make a visit",
      CAMEOEVENTCODE == "043" ~ "Host a visit",
      CAMEOEVENTCODE == "044" ~ "Meet at a third location",
      CAMEOEVENTCODE == "045" ~ "Mediate",
      CAMEOEVENTCODE == "046" ~ "Engage in negotiation",
      CAMEOEVENTCODE == "05" ~ "Engage in Diplomatic Cooperation",
      CAMEOEVENTCODE ==
        "050" ~ "Engage in diplomatic cooperation, not specified below",
      CAMEOEVENTCODE == "051" ~ "Praise or endorse",
      CAMEOEVENTCODE == "052" ~ "Defend verbally",
      CAMEOEVENTCODE == "053" ~ "Rally support on behalf of",
      CAMEOEVENTCODE == "054" ~ "Grant diplomatic recognition",
      CAMEOEVENTCODE == "055" ~ "Apologize",
      CAMEOEVENTCODE == "056" ~ "Forgive",
      CAMEOEVENTCODE == "057" ~ "Sign formal agreement",
      CAMEOEVENTCODE == "06" ~ "Engage in Material Cooperation",
      CAMEOEVENTCODE ==
        "060" ~ "Engage in material cooperation, not specified below",
      CAMEOEVENTCODE == "061" ~ "Cooperate economically",
      CAMEOEVENTCODE == "062" ~ "Cooperate militarily",
      CAMEOEVENTCODE == "063" ~ "Engage in judicial cooperation",
      CAMEOEVENTCODE == "064" ~ "Share intelligence or information",
      CAMEOEVENTCODE == "07" ~ "Provide Aid",
      CAMEOEVENTCODE == "070" ~ "Provide aid, not specified below",
      CAMEOEVENTCODE == "071" ~ "Provide economic aid",
      CAMEOEVENTCODE == "072" ~ "Provide military aid",
      CAMEOEVENTCODE == "073" ~ "Provide humanitarian aid",
      CAMEOEVENTCODE == "074" ~ "Provide military protection or peacekeeping",
      CAMEOEVENTCODE == "075" ~ "Grant asylum",
      CAMEOEVENTCODE == "08" ~ "Yield",
      CAMEOEVENTCODE == "080" ~ "Yield, not specified below",
      CAMEOEVENTCODE ==
        "081" ~ "Ease administrative sanctions, not specified below",
      CAMEOEVENTCODE == "0811" ~ "Ease restrictions on political freedoms",
      CAMEOEVENTCODE == "0812" ~ "Ease ban on political parties or politicians",
      CAMEOEVENTCODE == "0813" ~ "Ease curfew",
      CAMEOEVENTCODE == "0814" ~ "Ease state of emergency or martial law",
      CAMEOEVENTCODE == "082" ~ "Ease political dissent",
      CAMEOEVENTCODE ==
        "083" ~ "Accede to requests or demands for political reform, not specified below",
      CAMEOEVENTCODE == "0831" ~ "Accede to demands for change in leadership",
      CAMEOEVENTCODE == "0832" ~ "Accede to demands for change in policy",
      CAMEOEVENTCODE == "0833" ~ "Accede to demands for rights",
      CAMEOEVENTCODE ==
        "0834" ~ "Accede to demands for change in institutions, regime",
      CAMEOEVENTCODE == "084" ~ "Return, release, not specified below",
      CAMEOEVENTCODE == "0841" ~ "Return, release person(s)",
      CAMEOEVENTCODE == "0842" ~ "Return, release property",
      CAMEOEVENTCODE == "085" ~ "Ease economic sanctions, boycott, embargo",
      CAMEOEVENTCODE ==
        "086" ~ "Allow international involvement, not specified below",
      CAMEOEVENTCODE == "0861" ~ "Receive deployment of peacekeepers",
      CAMEOEVENTCODE == "0862" ~ "Receive inspectors",
      CAMEOEVENTCODE == "0863" ~ "Allow humanitarian access",
      CAMEOEVENTCODE == "087" ~ "De-escalate military engagement",
      CAMEOEVENTCODE == "0871" ~ "Declare truce, ceasefire",
      CAMEOEVENTCODE == "0872" ~ "Ease military blockade",
      CAMEOEVENTCODE == "0873" ~ "Demobilize armed forces",
      CAMEOEVENTCODE == "0874" ~ "Retreat or surrender militarily",
      CAMEOEVENTCODE == "09" ~ "Investigate",
      CAMEOEVENTCODE == "090" ~ "Investigate, not specified below",
      CAMEOEVENTCODE == "091" ~ "Investigate crime, corruption",
      CAMEOEVENTCODE == "092" ~ "Investigate human rights abuses",
      CAMEOEVENTCODE == "093" ~ "Investigate military action",
      CAMEOEVENTCODE == "094" ~ "Investigate war crimes",
      CAMEOEVENTCODE == "10" ~ "Demand",
      CAMEOEVENTCODE == "100" ~ "Demand, not specified below",
      CAMEOEVENTCODE ==
        "101" ~ "Demand material cooperation, not specified below",
      CAMEOEVENTCODE == "1011" ~ "Demand economic cooperation",
      CAMEOEVENTCODE == "1012" ~ "Demand military cooperation",
      CAMEOEVENTCODE == "1013" ~ "Demand judicial cooperation",
      CAMEOEVENTCODE == "1014" ~ "Demand intelligence cooperation",
      CAMEOEVENTCODE == "102" ~ "Demand diplomatic cooperation",
      CAMEOEVENTCODE == "103" ~ "Demand material aid, not specified below",
      CAMEOEVENTCODE == "1031" ~ "Demand economic aid",
      CAMEOEVENTCODE == "1032" ~ "Demand military aid",
      CAMEOEVENTCODE == "1033" ~ "Demand humanitarian aid",
      CAMEOEVENTCODE == "1034" ~ "Demand military protection or peacekeeping",
      CAMEOEVENTCODE == "104" ~ "Demand political reform, not specified below",
      CAMEOEVENTCODE == "1041" ~ "Demand change in leadership",
      CAMEOEVENTCODE == "1042" ~ "Demand policy change",
      CAMEOEVENTCODE == "1043" ~ "Demand rights",
      CAMEOEVENTCODE == "1044" ~ "Demand change in institutions, regime",
      CAMEOEVENTCODE ==
        "105" ~ "Demand that target yields, not specified below",
      CAMEOEVENTCODE == "1051" ~ "Demand easing of administrative sanctions",
      CAMEOEVENTCODE == "1052" ~ "Demand easing of political dissent",
      CAMEOEVENTCODE == "1053" ~ "Demand release of persons or property",
      CAMEOEVENTCODE ==
        "1054" ~ "Demand easing of economic sanctions, boycott, or embargo",
      CAMEOEVENTCODE ==
        "1055" ~ "Demand that target allows international involvement (non-mediation)",
      CAMEOEVENTCODE == "1056" ~ "Demand de-escalation of military engagement",
      CAMEOEVENTCODE == "107" ~ "Demand settling of dispute",
      CAMEOEVENTCODE == "108" ~ "Demand mediation",
      CAMEOEVENTCODE == "11" ~ "Disapprove",
      CAMEOEVENTCODE == "110" ~ "Disapprove, not specified below",
      CAMEOEVENTCODE == "111" ~ "Criticize or denounce",
      CAMEOEVENTCODE == "112" ~ "Accuse, not specified below",
      CAMEOEVENTCODE == "1121" ~ "Accuse of crime, corruption",
      CAMEOEVENTCODE == "1122" ~ "Accuse of human rights abuses",
      CAMEOEVENTCODE == "1123" ~ "Accuse of aggression",
      CAMEOEVENTCODE == "1124" ~ "Accuse of war crimes",
      CAMEOEVENTCODE == "1125" ~ "Accuse of espionage, treason",
      CAMEOEVENTCODE == "113" ~ "Rally opposition against",
      CAMEOEVENTCODE == "114" ~ "Complain officially",
      CAMEOEVENTCODE == "115" ~ "Bring lawsuit against",
      CAMEOEVENTCODE == "116" ~ "Find guilty or liable (legally)",
      CAMEOEVENTCODE == "12" ~ "Reject",
      CAMEOEVENTCODE == "120" ~ "Reject, not specified below",
      CAMEOEVENTCODE == "121" ~ "Reject material cooperation",
      CAMEOEVENTCODE == "1211" ~ "Reject economic cooperation",
      CAMEOEVENTCODE == "1212" ~ "Reject military cooperation",
      CAMEOEVENTCODE ==
        "122" ~ "Reject request or demand for material aid, not specified below",
      CAMEOEVENTCODE == "1221" ~ "Reject request for economic aid",
      CAMEOEVENTCODE == "1222" ~ "Reject request for military aid",
      CAMEOEVENTCODE == "1223" ~ "Reject request for humanitarian aid",
      CAMEOEVENTCODE ==
        "1224" ~ "Reject request for military protection or peacekeeping",
      CAMEOEVENTCODE ==
        "123" ~ "Reject request or demand for political reform, not specified below",
      CAMEOEVENTCODE == "1231" ~ "Reject request for change in leadership",
      CAMEOEVENTCODE == "1232" ~ "Reject request for policy change",
      CAMEOEVENTCODE == "1233" ~ "Reject request for rights",
      CAMEOEVENTCODE ==
        "1234" ~ "Reject request for change in institutions, regime",
      CAMEOEVENTCODE == "124" ~ "Refuse to yield, not specified below",
      CAMEOEVENTCODE == "1241" ~ "Refuse to ease administrative sanctions",
      CAMEOEVENTCODE == "1242" ~ "Refuse to ease popular dissent",
      CAMEOEVENTCODE == "1243" ~ "Refuse to release persons or property",
      CAMEOEVENTCODE ==
        "1244" ~ "Refuse to ease economic sanctions, boycott, or embargo",
      CAMEOEVENTCODE ==
        "1245" ~ "Refuse to allow international involvement (non-mediation)",
      CAMEOEVENTCODE == "1246" ~ "Refuse to de-escalate military engagement",
      CAMEOEVENTCODE ==
        "125" ~ "Reject proposal to meet, discuss, or negotiate",
      CAMEOEVENTCODE == "126" ~ "Reject mediation",
      CAMEOEVENTCODE == "127" ~ "Reject plan, agreement to settle dispute",
      CAMEOEVENTCODE == "128" ~ "Defy norms, law",
      CAMEOEVENTCODE == "129" ~ "Veto",
      CAMEOEVENTCODE == "13" ~ "Threaten",
      CAMEOEVENTCODE == "130" ~ "Threaten, not specified below",
      CAMEOEVENTCODE == "131" ~ "Threaten non-force, not specified below",
      CAMEOEVENTCODE == "1311" ~ "Threaten to reduce or stop aid",
      CAMEOEVENTCODE == "1312" ~ "Threaten with sanctions, boycott, embargo",
      CAMEOEVENTCODE == "1313" ~ "Threaten to reduce or break relations",
      CAMEOEVENTCODE ==
        "132" ~ "Threaten with administrative sanctions, not specified below",
      CAMEOEVENTCODE ==
        "1321" ~ "Threaten with restrictions on political freedoms",
      CAMEOEVENTCODE ==
        "1322" ~ "Threaten to ban political parties or politicians",
      CAMEOEVENTCODE == "1323" ~ "Threaten to impose curfew",
      CAMEOEVENTCODE ==
        "1324" ~ "Threaten to impose state of emergency or martial law",
      CAMEOEVENTCODE == "133" ~ "Threaten with political dissent, protest",
      CAMEOEVENTCODE == "134" ~ "Threaten to halt negotiations",
      CAMEOEVENTCODE == "135" ~ "Threaten to halt mediation",
      CAMEOEVENTCODE ==
        "136" ~ "Threaten to halt international involvement (non-mediation)",
      CAMEOEVENTCODE == "137" ~ "Threaten with repression",
      CAMEOEVENTCODE ==
        "138" ~ "Threaten with military force, not specified below",
      CAMEOEVENTCODE == "1381" ~ "Threaten blockade",
      CAMEOEVENTCODE == "1382" ~ "Threaten occupation",
      CAMEOEVENTCODE == "1383" ~ "Threaten unconventional violence",
      CAMEOEVENTCODE == "1384" ~ "Threaten conventional attack",
      CAMEOEVENTCODE == "1385" ~ "Threaten attack with WMD",
      CAMEOEVENTCODE == "139" ~ "Give ultimatum",
      CAMEOEVENTCODE == "14" ~ "Protest",
      CAMEOEVENTCODE ==
        "140" ~ "Engage in political dissent, not specified below",
      CAMEOEVENTCODE == "141" ~ "Demonstrate or rally, not specified below",
      CAMEOEVENTCODE == "1411" ~ "Demonstrate for leadership change",
      CAMEOEVENTCODE == "1412" ~ "Demonstrate for policy change",
      CAMEOEVENTCODE == "1413" ~ "Demonstrate for rights",
      CAMEOEVENTCODE ==
        "1414" ~ "Demonstrate for change in institutions, regime",
      CAMEOEVENTCODE == "142" ~ "Conduct hunger strike, not specified below",
      CAMEOEVENTCODE == "1421" ~ "Conduct hunger strike for leadership change",
      CAMEOEVENTCODE == "1422" ~ "Conduct hunger strike for policy change",
      CAMEOEVENTCODE == "1423" ~ "Conduct hunger strike for rights",
      CAMEOEVENTCODE ==
        "1424" ~ "Conduct hunger strike for change in institutions, regime",
      CAMEOEVENTCODE ==
        "143" ~ "Conduct strike or boycott, not specified below",
      CAMEOEVENTCODE ==
        "1431" ~ "Conduct strike or boycott for leadership change",
      CAMEOEVENTCODE == "1432" ~ "Conduct strike or boycott for policy change",
      CAMEOEVENTCODE == "1433" ~ "Conduct strike or boycott for rights",
      CAMEOEVENTCODE ==
        "1434" ~ "Conduct strike or boycott for change in institutions, regime",
      CAMEOEVENTCODE == "144" ~ "Obstruct passage, block, not specified below",
      CAMEOEVENTCODE == "1441" ~ "Obstruct passage to demand leadership change",
      CAMEOEVENTCODE == "1442" ~ "Obstruct passage to demand policy change",
      CAMEOEVENTCODE == "1443" ~ "Obstruct passage to demand rights",
      CAMEOEVENTCODE ==
        "1444" ~ "Obstruct passage to demand change in institutions, regime",
      CAMEOEVENTCODE == "145" ~ "Protest violently, riot, not specified below",
      CAMEOEVENTCODE ==
        "1451" ~ "Engage in violent protest for leadership change",
      CAMEOEVENTCODE == "1452" ~ "Engage in violent protest for policy change",
      CAMEOEVENTCODE == "1453" ~ "Engage in violent protest for rights",
      CAMEOEVENTCODE ==
        "1454" ~ "Engage in violent protest for change in institutions, regime",
      CAMEOEVENTCODE == "15" ~ "Exhibit Force Posture",
      CAMEOEVENTCODE ==
        "150" ~ "Demonstrate military or police power, not specified below",
      CAMEOEVENTCODE == "151" ~ "Increase police alert status",
      CAMEOEVENTCODE == "152" ~ "Increase military alert status",
      CAMEOEVENTCODE == "153" ~ "Mobilize or increase police power",
      CAMEOEVENTCODE == "154" ~ "Mobilize or increase armed forces",
      CAMEOEVENTCODE == "155" ~ "Mobilize or increase cyber-forces",
      CAMEOEVENTCODE == "16" ~ "Reduce Relations",
      CAMEOEVENTCODE == "160" ~ "Reduce relations, not specified below",
      CAMEOEVENTCODE == "161" ~ "Reduce or break diplomatic relations",
      CAMEOEVENTCODE ==
        "162" ~ "Reduce or stop material aid, not specified below",
      CAMEOEVENTCODE == "1621" ~ "Reduce or stop economic assistance",
      CAMEOEVENTCODE == "1622" ~ "Reduce or stop military assistance",
      CAMEOEVENTCODE == "1623" ~ "Reduce or stop humanitarian assistance",
      CAMEOEVENTCODE == "163" ~ "Impose embargo, boycott, or sanctions",
      CAMEOEVENTCODE == "164" ~ "Halt negotiations",
      CAMEOEVENTCODE == "165" ~ "Halt mediation",
      CAMEOEVENTCODE == "166" ~ "Expel or withdraw, not specified below",
      CAMEOEVENTCODE == "1661" ~ "Expel or withdraw peacekeepers",
      CAMEOEVENTCODE == "1662" ~ "Expel or withdraw inspectors, observers",
      CAMEOEVENTCODE == "1663" ~ "Expel or withdraw aid agencies",
      CAMEOEVENTCODE == "17" ~ "Coerce",
      CAMEOEVENTCODE == "170" ~ "Coerce, not specified below",
      CAMEOEVENTCODE == "171" ~ "Seize or damage property, not specified below",
      CAMEOEVENTCODE == "1711" ~ "Confiscate property",
      CAMEOEVENTCODE == "1712" ~ "Destroy property",
      CAMEOEVENTCODE ==
        "172" ~ "Impose administrative sanctions, not specified below",
      CAMEOEVENTCODE == "1721" ~ "Impose restrictions on political freedoms",
      CAMEOEVENTCODE == "1722" ~ "Ban political parties or politicians",
      CAMEOEVENTCODE == "1723" ~ "Impose curfew",
      CAMEOEVENTCODE == "1724" ~ "Impose state of emergency or martial law",
      CAMEOEVENTCODE == "173" ~ "Arrest, detain, or charge with legal action",
      CAMEOEVENTCODE == "174" ~ "Expel or deport individuals",
      CAMEOEVENTCODE == "175" ~ "Use tactics of violent repression",
      CAMEOEVENTCODE == "176" ~ "Attack cybernetically",
      CAMEOEVENTCODE == "18" ~ "Assault",
      CAMEOEVENTCODE ==
        "180" ~ "Use unconventional violence, not specified below",
      CAMEOEVENTCODE == "181" ~ "Abduct, hijack, or take hostage",
      CAMEOEVENTCODE == "182" ~ "Physically assault, not specified below",
      CAMEOEVENTCODE == "1821" ~ "Sexually assault",
      CAMEOEVENTCODE == "1822" ~ "Torture",
      CAMEOEVENTCODE == "1823" ~ "Kill by physical assault",
      CAMEOEVENTCODE ==
        "183" ~ "Conduct suicide, car, or other non-military bombing, not specified below",
      CAMEOEVENTCODE == "1831" ~ "Carry out suicide bombing",
      CAMEOEVENTCODE == "1832" ~ "Carry out vehicular bombing",
      CAMEOEVENTCODE == "1833" ~ "Carry out roadside bombing",
      CAMEOEVENTCODE == "1834" ~ "Carry out location bombing",
      CAMEOEVENTCODE == "184" ~ "Use as human shield",
      CAMEOEVENTCODE == "185" ~ "Attempt to assassinate",
      CAMEOEVENTCODE == "186" ~ "Assassinate",
      CAMEOEVENTCODE == "19" ~ "Fight",
      CAMEOEVENTCODE ==
        "190" ~ "Use conventional military force, not specified below",
      CAMEOEVENTCODE == "191" ~ "Impose blockade, restrict movement",
      CAMEOEVENTCODE == "192" ~ "Occupy territory",
      CAMEOEVENTCODE == "193" ~ "Fight with small arms and light weapons",
      CAMEOEVENTCODE == "194" ~ "Fight with artillery and tanks",
      CAMEOEVENTCODE == "195" ~ "Employ aerial weapons, not specified below",
      CAMEOEVENTCODE == "1951" ~ "Employ precision-guided aerial munitions",
      CAMEOEVENTCODE == "1952" ~ "Employ remotely piloted aerial munitions",
      CAMEOEVENTCODE == "196" ~ "Violate ceasefire",
      CAMEOEVENTCODE == "20" ~ "Use Unconventional Mass Violence",
      CAMEOEVENTCODE ==
        "200" ~ "Use unconventional mass violence, not specified below",
      CAMEOEVENTCODE == "201" ~ "Engage in mass expulsion",
      CAMEOEVENTCODE == "202" ~ "Engage in mass killings",
      CAMEOEVENTCODE == "203" ~ "Engage in ethnic cleansing",
      CAMEOEVENTCODE ==
        "204" ~ "Use weapons of mass destruction, not specified below",
      CAMEOEVENTCODE ==
        "2041" ~ "Use chemical, biological, or radiological weapons",
      CAMEOEVENTCODE == "2042" ~ "Detonate nuclear weapons",
      TRUE ~ NA_character_
    )
  )


cameo_df <- cameo_df %>%
  select(CAMEOEVENTCODE, CAMEOLabel, everything())


gdelt_cameos <- gdelt_filtered_full %>% group_by(EventCode) %>% count()

cameo_df <- cameo_df %>%
  left_join(
    gdelt_cameos,
    by = c("CAMEOEVENTCODE" = "EventCode")
  )

saveRDS(cameo_df, "data/cameo_df.rds")
