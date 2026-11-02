# trim_theta_before_first_obs.R
#
# Drops dyad-periods before each dyad's first observed interaction
# (is_obs == 1) from an extract_theta() tibble, then slims and saves the
# result as an xz-compressed .rda for data/.
#
# Rationale: the bilatr likelihood is only evaluated where is_obs[d, t] == 1,
# so theta before a dyad's first event is pure random-walk prior smoothed
# backwards from the first observation -- it carries no information from
# that dyad's data.
#
# Requires the stan_data.rds the theta file was extracted from (it carries
# the D x T `is_obs` matrix, indexed by the same dyad_id/time_index).

trim_theta_before_first_obs <- function(theta, stan_data, trim_after_last = FALSE) {
  is_obs <- stan_data$is_obs
  stopifnot(
    is.matrix(is_obs),
    max(theta$dyad_id) <= nrow(is_obs),
    max(theta$time_index) <= ncol(is_obs)
  )

  obs_rows <- which(is_obs == 1L, arr.ind = TRUE) # cols: row = dyad_id, col = time_index
  bounds <- data.frame(
    dyad_id = seq_len(nrow(is_obs)),
    first_t = tapply(obs_rows[, "col"], factor(obs_rows[, "row"], levels = seq_len(nrow(is_obs))), min),
    last_t  = tapply(obs_rows[, "col"], factor(obs_rows[, "row"], levels = seq_len(nrow(is_obs))), max)
  )
  # min_n_events >= 1 guarantees every retained dyad has at least one event
  stopifnot(!anyNA(bounds$first_t))

  i <- match(theta$dyad_id, bounds$dyad_id)
  keep <- theta$time_index >= bounds$first_t[i]
  if (trim_after_last) keep <- keep & theta$time_index <= bounds$last_t[i]

  message(sprintf(
    "kept %s of %s dyad-periods (%.1f%%)",
    format(sum(keep), big.mark = ","), format(length(keep), big.mark = ","), 100 * mean(keep)
  ))
  theta[keep, , drop = FALSE]
}

# Slim columns for shipping: dyad_id/time_index are recoverable from
# dyad/year, dyad2 from dyad; factors and integers compress much better
# than repeated character strings and doubles.
slim_theta <- function(theta, digits = 3) {
  out <- data.frame(
    dyad = factor(theta$dyad),
    year = as.integer(theta$year),
    mean = signif(theta$mean, digits),
    q05  = signif(theta$`5%`, digits),
    q50  = signif(theta$`50%`, digits),
    q95  = signif(theta$`95%`, digits)
  )
  if ("month" %in% names(theta)) out$month <- as.integer(theta$month)
  out[order(out$dyad, out$year), ]
}

# ---- run ----------------------------------------------------------------
theta_path     <- "runscripts/thetatests/stable_gamma_gdelt_bilatr_ERC16NZ_ref5_y_dir_1991_1_wnone_theta.rds"
stan_data_path <- "path/to/stable_gamma_gdelt_bilatr_ERC16NZ_ref5_y_dir_1991_1_wnone_stan_data.rds"
obj_name       <- "theta_dir_smoothed"   # name of the object in data/

theta     <- readRDS(theta_path)
stan_data <- readRDS(stan_data_path)

theta_trim <- slim_theta(trim_theta_before_first_obs(theta, stan_data))

assign(obj_name, theta_trim)
save(list = obj_name, file = file.path("data", paste0(obj_name, ".rda")),
     compress = "xz", compression_level = 9)
