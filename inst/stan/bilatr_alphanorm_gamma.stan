// bilatr model `stable_gamma` (experimental): the `stable` model
// (bilatr_alphanorm.stan) plus country-level category offsets gamma.
// Fit with fit_panel(stan_model = "stable_gamma").
//
// Everything in bilatr_alphanorm.stan's header applies here unchanged: the
// likelihood, the random-walk latent process, the priors, and the
// identification of alpha, theta and mu_intercept, including the
// orientation fold. This header covers what gamma adds.
//
// PURPOSE
//   Some countries use some action classes more (or less) than their
//   dyads' theta implies, whomever they interact with. A dyad-level offset
//   cannot capture this, because each dyad already has a free theta; a
//   country-level offset can, because each country appears in many dyads.
//   gamma[, c] is country c's offset vector over the A action classes, so
//   theta is measured net of each country's baseline mix of event types.
//
// ADDITIONAL DATA
//   n_countries    number of distinct countries in the retained dyads
//   ctry_a[d]      the sender (directed) or first country (undirected)
//   ctry_b[d]      the receiver (directed) or second country (undirected)
//   w_send[d]      1 for directed data; for undirected data, the share of
//                  the pair's events with ctry_a[d] as the acting country
//
// LINEAR PREDICTOR
//   eta[d, t] = alpha * theta[d, t] - mu_eff[d],  mu_eff[d] = mu_intercept + g_d
//   - directed, d = (i -> j): g_d = gamma[, i] (the sender generates the
//     events, so its repertoire applies);
//   - undirected, d = {i, j}: g_d = w_send[d] * gamma[, i]
//     + (1 - w_send[d]) * gamma[, j].
//   One program serves both: w_send = 1 reduces the undirected formula to
//   the directed one exactly, so directed and undirected fits estimate the
//   same gamma, with the same prior and interpretation. g_d is constant
//   within a dyad and is folded into the intercept once per dyad
//   (partial_log_lik_offset()). The models nest: as sigma_gamma -> 0, this
//   program and `stable` coincide, and with n_countries = 1, gamma is
//   identically 0.
//
// ADDITIONAL PARAMETERS AND PRIORS
//   gamma_z[a, c]  ~ std_normal()          (non-centered)
//   sigma_gamma[a] ~ half-normal(0, 0.3)   (per-class spread across countries)
//   The prior scale reflects a dyad-level compositional residual with an
//   RMS of about 0.22 on the linear-predictor scale, found by
//   check_compositional_residuals() in `stable` fits: a country-level
//   effect can only be part of such a residual, so 0.22 is an upper bound
//   on it.
//
// IDENTIFICATION: GAMMA CONSTRAINTS. Three constraints on gamma, enforced
// exactly by construction in transformed parameters:
//   1. Orthogonal to alpha, per country. Shifting gamma[, c] by k * alpha
//      shifts every dyad containing c by exactly what that dyad's own free
//      theta can absorb, so that component is not identified. It is
//      projected out with dot_self(alpha) rather than a hard-coded A
//      (equal only because alpha has RMS 1), so the projection does not
//      depend on alpha's normalization.
//   2. Sum to zero over classes, per country. A constant added to every
//      class of gamma[, c] cancels in softmax().
//   3. Sum to zero over countries, per class. The country-independent part
//      of a class's offset is mu_intercept's.
// gamma is built in four steps: scale gamma_z by sigma_gamma; center each
// class across countries (3); center each country across classes (2);
// project out alpha per country (1). Each later step preserves the
// earlier steps' zero sums, so all three constraints hold simultaneously
// (the package tests check this numerically). The directions of gamma_z
// that the construction projects out are informed by its prior only, so
// gamma_z is left out of the convergence diagnostics, like the
// sign-ambiguous raw parameters.
//
// gamma and the orientation fold. gamma is not multiplied by
// orientation_sign(alpha_raw): the projection
// v - (dot_product(v, alpha) / dot_self(alpha)) * alpha removes the same
// component of v whether alpha or -alpha is used, so gamma, g_d and eta
// are unchanged when the fold flips alpha and theta.
//
// HONEST CAVEAT (mixture approximation). If undirected counts are the
// sum of two directed processes, the exact likelihood is a mixture of
// two softmax compositions, not a single softmax with averaged offsets.
// Averaging gamma_i/gamma_j in the linear predictor (as g_d does above)
// is a first-order approximation, accurate when gamma is small -- which
// the residual check this model is built to test suggests it is. The
// exact two-component mixture breaks the Dirichlet-multinomial form and
// is not worth the cost; this is acknowledged here deliberately, not an
// oversight.
//
// GENERATED QUANTITIES
//   As in `stable`: log_lik (when compute_log_lik = 1) and theta_filtered /
//   theta_filtered_sd (when compute_theta_filtered = 1), both using mu_eff
//   in place of mu_intercept.
functions {
  // <<< BEGIN GENERATED partial_log_lik (source: inst/stan/include/partial_log_lik.stanfunctions) >>>
  // Do not hand-edit between these markers -- edit the source file
  // above and rerun `Rscript data-raw/sync_stan_functions.R`
  // (checked by tests/testthat/test_stan_includes.R).
  // Shared reduce_sum likelihood for all bilatr Stan model variants.
  //
  // Indifferent to how `theta`, `alpha` and `mu_intercept` were constructed
  // upstream (random walk or OU process, with or without country offsets):
  // it only consumes the already-built `theta`/`alpha`/`mu_intercept`/`phi`
  // values.
  //
  // This is the single canonical source. It is not included via Stan's
  // `#include`, which cmdstanr cannot resolve at sampling time when the
  // include path contains a space
  // (https://github.com/stan-dev/cmdstanr/issues/820). Instead, its contents
  // are copied verbatim into a marker-delimited block in each model's
  // `.stan` file, and the package tests check that the copies match.
  real partial_log_lik(array[] int slice_d,
                        int start, int end,
                        int T, int A,
                        array[,] int is_obs,
                        array[,,] int Y,
                        array[,] real theta,
                        vector mu_intercept,
                        vector phi,
                        vector alpha) {
    real lp = 0;
    for (d in start:end) {
      for (t in 1:T) {
        if (is_obs[d, t] == 1) {
          vector[A] eta = alpha .* rep_vector(theta[d, t], A) - mu_intercept;
          vector[A] p = softmax(eta);
          vector[A] conc = phi[d] * p;
          lp += dirichlet_multinomial_lpmf(Y[d, t] | conc);
        }
      }
    }
    return lp;
  }

  // Country-offset variant of partial_log_lik(), used by the `stable_gamma`
  // model (and present, unused, in the other models' copies of this file).
  //
  // The per-dyad country offset g_d is built once per dyad, outside the `t`
  // loop, and folded into the intercept (mu_eff = mu_intercept + g_d), so
  // the per-cell expression has the same shape and autodiff cost as
  // partial_log_lik()'s. g_d is not stored: a per-dyad matrix in
  // transformed parameters would be written to the output for every draw.
  // w_send[d] == 1 covers directed data (g_d = gamma[, ctry_a[d]] exactly,
  // ctry_b[d] unused); any other value mixes ctry_a's and ctry_b's columns
  // by event share (see the `stable_gamma` model's header).
  real partial_log_lik_offset(array[] int slice_d,
                               int start, int end,
                               int T, int A,
                               array[,] int is_obs,
                               array[,,] int Y,
                               array[,] real theta,
                               vector mu_intercept,
                               vector phi,
                               vector alpha,
                               matrix gamma,
                               array[] int ctry_a,
                               array[] int ctry_b,
                               vector w_send) {
    real lp = 0;
    for (d in start:end) {
      // fold the country offset into the intercept ONCE PER DYAD: the
      // per-cell expression below is then character-for-character the
      // same shape, and the same autodiff cost, as partial_log_lik()'s.
      vector[A] mu_eff = mu_intercept + (w_send[d] == 1.0
                         ? gamma[, ctry_a[d]]
                         : w_send[d] * gamma[, ctry_a[d]] + (1 - w_send[d]) * gamma[, ctry_b[d]]);
      for (t in 1:T) {
        if (is_obs[d, t] == 1) {
          vector[A] eta = alpha .* rep_vector(theta[d, t], A) - mu_eff;
          vector[A] p = softmax(eta);
          vector[A] conc = phi[d] * p;
          lp += dirichlet_multinomial_lpmf(Y[d, t] | conc);
        }
      }
    }
    return lp;
  }
  // <<< END GENERATED partial_log_lik >>>

  // +1/-1 orientation of a draw, from the reference class's sign. Applied
  // to alpha and theta below (NOT to gamma -- see header); see
  // bilatr_alphanorm.stan's header, "Sign: the ORIENTATION FOLD".
  real orientation_sign(vector alpha_raw) {
    return alpha_raw[1] >= 0 ? 1.0 : -1.0;
  }

  // Per-dyad country offset g_d -- the same construction as in
  // partial_log_lik_offset() above, for use in generated quantities.
  vector country_offset(matrix gamma, int ctry_a_d, int ctry_b_d, real w_send_d) {
    return w_send_d == 1.0
           ? gamma[, ctry_a_d]
           : w_send_d * gamma[, ctry_a_d] + (1 - w_send_d) * gamma[, ctry_b_d];
  }

  // Offset-aware per-dyad-period log-likelihood, not summed -- used only
  // by the compute_log_lik generated quantities below. Its `mu_intercept`
  // argument is mu_eff, already folded with g_d by the caller once per
  // dyad, so the body is the same as bilatr_alphanorm.stan's
  // dyad_period_log_lik(); the separate name marks the call sites.
  real dyad_period_log_lik_offset(int obs_dt, array[] int y_dt, real theta_dt,
                                   int A, vector mu_intercept, real phi_d,
                                   vector alpha) {
    if (obs_dt == 0) {
      return 0;
    }
    vector[A] eta = alpha .* rep_vector(theta_dt, A) - mu_intercept;
    vector[A] p = softmax(eta);
    vector[A] conc = phi_d * p;
    return dirichlet_multinomial_lpmf(y_dt | conc);
  }
}
data {
  int<lower=1> T;                            // number of time points
  int<lower=1> D;                            // number of dyads
  int<lower=1> A;                            // number of action types
  int<lower=1> C;                            // reduce_sum grainsize (NOT the
                                              // country count; that is
                                              // n_countries, below)
  array[D, T] int<lower=0, upper=1> is_obs;  // observed indicator
  array[D, T, A] int<lower=0> Y;             // event counts
  int<lower=0, upper=1> compute_log_lik;     // 1 = also compute per-dyad-period
                                              // log_lik in generated quantities
                                              // (D x T x draws; default 0/off)
  int<lower=0, upper=1> prior_only;          // 1 = skip the likelihood entirely
                                              // (fit the prior only; see
                                              // alpha_prior_moments())
  int<lower=0, upper=1> compute_theta_filtered; // 1 = also compute
                                              // theta_filtered/theta_filtered_sd
                                              // in generated quantities, for the
                                              // dyads in filter_dyads
  int<lower=0> n_filter_dyads;               // length of filter_dyads; 0 if
                                              // compute_theta_filtered is 0
  array[n_filter_dyads] int<lower=1, upper=D> filter_dyads; // which dyads to filter

  // Country-offset data (see header). n_countries is the number of
  // DISTINCT countries present in the retained dyads (set by
  // prepare_stan_data()) -- deliberately NOT named `C` (already the
  // reduce_sum grainsize above).
  int<lower=1> n_countries;
  array[D] int<lower=1, upper=n_countries> ctry_a;   // sender (directed) / side A (undirected)
  array[D] int<lower=1, upper=n_countries> ctry_b;   // receiver (directed) / side B
  vector<lower=0, upper=1>[D] w_send;                // 1 if directed; side-A event share if undirected
}
parameters {
  // Latent states per dyad
  array[D, T] real theta_raw;
  sum_to_zero_vector[A] mu_intercept;   // softmax level-shift; sums to 0 exactly
  sum_to_zero_vector[A] alpha_raw;      // pre-normalization discrimination; sums to 0 exactly
  real<lower=0> sigma_theta0;
  vector[D] z_theta0;                   // no mu_theta0: theta0 has mean exactly 0

  // Hierarchical process/dispersion parameters
  vector[D] log_process_noise_raw;  // non-centered: process_noise built in TP
  real mu_log_noise;                // log ratio to sigma_theta0, see header
  real<lower=0> sigma_log_noise;
  vector<lower=0>[D] phi;
  real mu_log_phi;
  real<lower=0> sigma_log_phi;

  // Country-level category offsets (see header)
  matrix[A, n_countries] gamma_z;      // non-centered
  vector<lower=0>[A] sigma_gamma;      // per-category spread across countries
}
transformed parameters {
  // Normalize alpha to RMS (== population SD, sum-to-zero) 1. This is the
  // scale anchor for the whole model; see bilatr_alphanorm.stan's header,
  // "Scale" and "Radial direction". The sign comes from
  // orientation_sign(); see that header, "Sign: the ORIENTATION FOLD".
  vector[A] alpha = orientation_sign(alpha_raw)
                    * alpha_raw * sqrt((1.0 * A) / dot_self(alpha_raw));

  // Country-level category offsets. See header, "IDENTIFICATION: GAMMA
  // CONSTRAINTS" -- the four-step order matters: step 2 (row-center
  // across countries) leaves step 3's per-column centering unable to
  // reintroduce a row mean (the constants step 3 subtracts average to
  // zero across countries, since each row already sums to zero before
  // step 3 runs); step 4's alpha-projection subtracts a multiple of
  // alpha, which itself sums to zero (built from a sum_to_zero_vector),
  // so it cannot perturb the column (per-country) sums step 3 already
  // zeroed, and its own per-country coefficients average to zero across
  // countries so it cannot reintroduce a row mean either. After all four
  // steps, constraints 1-3 hold EXACTLY, simultaneously.
  matrix[A, n_countries] gamma;
  {
    // 1. raw = diag_pre_multiply(sigma_gamma, gamma_z)
    matrix[A, n_countries] raw = diag_pre_multiply(sigma_gamma, gamma_z);

    // 2. subtract each row's (category's) mean across countries
    matrix[A, n_countries] row_centred;
    for (k in 1:A) {
      row_centred[k] = raw[k] - mean(raw[k]);
    }

    // 3. subtract each column's (country's) mean across categories
    matrix[A, n_countries] col_centred;
    for (c in 1:n_countries) {
      col_centred[, c] = row_centred[, c] - mean(row_centred[, c]);
    }

    // 4. project out alpha, per country -- dot_self(alpha), not A: see
    // header (correct in general, not just when alpha's RMS is 1)
    for (c in 1:n_countries) {
      vector[A] v = col_centred[, c];
      gamma[, c] = v - (dot_product(v, alpha) / dot_self(alpha)) * alpha;
    }
  }

  // theta0 folds the same sign as alpha, so alpha .* theta (all the
  // likelihood ever sees) is invariant to which raw-space basin the
  // sampler is in -- see bilatr_alphanorm.stan's header. gamma is NOT
  // folded (see this file's header) -- it needs no orientation_sign()
  // multiplication anywhere.
  vector[D] theta0 = orientation_sign(alpha_raw) * sigma_theta0 * z_theta0;

  // Per-dyad process noise, a RATIO to sigma_theta0 (non-centered
  // lognormal hierarchy on the ratio, not on an absolute theta-unit
  // quantity -- see bilatr_alphanorm.stan's header).
  vector<lower=0>[D] process_noise =
    sigma_theta0 * exp(mu_log_noise + sigma_log_noise * log_process_noise_raw);

  array[D, T] real theta;
  for (d in 1:D) {
    theta[d, 1] = theta0[d] + orientation_sign(alpha_raw) * process_noise[d] * theta_raw[d, 1];
    for (t in 2:T) {
      theta[d, t] = theta[d, t - 1] + orientation_sign(alpha_raw) * process_noise[d] * theta_raw[d, t];
    }
  }
}
model {
  // hyperpriors for pooling (lognormal hierarchy on the process_noise /
  // sigma_theta0 ratio)
  mu_log_noise ~ normal(log(0.2), 0.5);
  sigma_log_noise ~ normal(0, 0.5);
  log_process_noise_raw ~ std_normal();

  mu_log_phi ~ normal(0, 1);
  sigma_log_phi ~ normal(0, 1);
  phi ~ lognormal(mu_log_phi, sigma_log_phi);

  // theta random-walk innovations
  for (d in 1:D) {
    for (t in 1:T) {
      theta_raw[d, t] ~ std_normal();
    }
  }

  sigma_theta0 ~ normal(0, 2);
  z_theta0 ~ std_normal();

  mu_intercept ~ std_normal();

  // load-bearing, not merely regularizing: identifies the radial
  // direction of alpha_raw (see bilatr_alphanorm.stan's header, "Radial
  // direction"). The orientation fold (transformed parameters, above)
  // leaves this prior untouched.
  alpha_raw ~ std_normal();

  // Country-level category offsets (see header, "ADDITIONAL PARAMETERS
  // AND PRIORS"). gamma_z ~ std_normal() (via to_vector(), element-wise):
  // gamma_z's projected-out directions (the alpha direction, the all-ones
  // direction, and the across-country mean) are informed ONLY by this
  // prior, not by the likelihood, which is why the convergence
  // diagnostics leave gamma_z out (not sign ambiguity, but the same
  // treatment).
  sigma_gamma ~ normal(0, 0.3);
  to_vector(gamma_z) ~ std_normal();

  // likelihood, chunked via reduce_sum -- skipped entirely if prior_only,
  // gating this exactly the way compute_log_lik gates its own generated
  // quantity below
  if (!prior_only) {
    array[D] int dyad_seq = linspaced_int_array(D, 1, D);
    target += reduce_sum(partial_log_lik_offset, dyad_seq, C,
                          T, A, is_obs, Y, theta, mu_intercept, phi, alpha,
                          gamma, ctry_a, ctry_b, w_send);
  }
}
generated quantities {
  array[compute_log_lik ? D : 0, compute_log_lik ? T : 0] real log_lik;

  if (compute_log_lik) {
    for (d in 1:D) {
      // fold the offset once per dyad, as in partial_log_lik_offset(), so
      // every site that uses g_d reads the same way.
      vector[A] mu_eff = mu_intercept + country_offset(gamma, ctry_a[d], ctry_b[d], w_send[d]);
      for (t in 1:T) {
        log_lik[d, t] = dyad_period_log_lik_offset(
          is_obs[d, t], Y[d, t], theta[d, t], A, mu_eff, phi[d], alpha
        );
      }
    }
  }

  // Forward-filtered theta: conditional on THIS draw's hyperparameters,
  // a Fisher-scoring (West-Harrison linear-Bayes) filter over the
  // observations -- an approximation to p(theta_t | y_1:t), not the exact
  // marginal. See bilatr_alphanorm.stan's generated quantities for the
  // derivation of the score/information; the offset changes exactly one thing here
  // -- eta uses mu_eff (mu_intercept folded with g_d) in place of
  // mu_intercept -- because d(eta_k)/d(theta) = alpha_k regardless of
  // g_d (g_d does not depend on theta), so the score and information
  // keep their form and are simply evaluated at
  // p = softmax(alpha * m_pred - mu_eff). mu_eff is computed once per
  // filtered dyad, outside the t loop, same reasoning as
  // partial_log_lik_offset()'s per-dyad fold.
  array[compute_theta_filtered ? n_filter_dyads : 0,
        compute_theta_filtered ? T : 0] real theta_filtered;
  array[compute_theta_filtered ? n_filter_dyads : 0,
        compute_theta_filtered ? T : 0] real theta_filtered_sd;

  if (compute_theta_filtered) {
    for (i in 1:n_filter_dyads) {
      int d = filter_dyads[i];
      vector[A] mu_eff = mu_intercept + country_offset(gamma, ctry_a[d], ctry_b[d], w_send[d]);
      real m = 0;
      real p_var = square(sigma_theta0);
      for (t in 1:T) {
        real m_pred = m;
        real p_pred = p_var + square(process_noise[d]);
        if (is_obs[d, t] == 1) {
          vector[A] eta = alpha .* rep_vector(m_pred, A) - mu_eff;
          vector[A] p = softmax(eta);
          real a_bar = dot_product(p, alpha);
          int n = sum(Y[d, t]);
          real info = n * (1 + phi[d]) / (n + phi[d]) * dot_product(p, square(alpha - a_bar));
          real g = 0;
          for (k in 1:A) {
            real conc_k = phi[d] * p[k];
            g += conc_k * (alpha[k] - a_bar) * (digamma(Y[d, t, k] + conc_k) - digamma(conc_k));
          }
          p_var = 1 / (1 / p_pred + info);
          m = m_pred + p_var * g;
        } else {
          m = m_pred;
          p_var = p_pred;
        }
        theta_filtered[i, t] = m;
        theta_filtered_sd[i, t] = sqrt(p_var);
      }
    }
  }
}
