// bilatr: hierarchical dynamic IRT model for dyadic conflict trajectories,
// with a country-level category-offset vector gamma_c added on top of
// `stable` (inst/stan/bilatr_alphanorm.stan).
//
// Registered as `stable_gamma` (R/model_registry.R, status =
// "experimental") since 0.7.0. Fit via
// fit_dyad_ts_dev()/fit_panel_dev(stan_model = "stable_gamma"); not
// reachable via the exported fit_dyad_ts()/fit_panel() (those always fit
// `.BILATR_DEFAULT_MODEL`).
//
// EVERYTHING in bilatr_alphanorm.stan's header applies here unchanged --
// alpha/theta/mu_intercept identification, the radial degeneracy, the
// orientation fold, and their rationale are exactly `stable`'s (this file
// starts as a copy of that one). This header covers only what gamma adds.
//
// MOTIVATION (see dev/claude_code_prompt_0.7.0_country_offsets.md for the
// full derivation this file implements): check_compositional_residuals()
// (0.6.0) found a dyad-level compositional residual (implied RMS ~= 0.217
// on the linear-predictor scale) that a dyad-level beta_d cannot resolve
// -- each dyad already has a free theta, so beta_d is not identified, and
// the check itself showed the per-dyad statistic is nearly powerless. A
// COUNTRY-level offset is identifiable: a country appears in dozens to
// hundreds of dyads.
//
// LINEAR PREDICTOR: eta_dt = alpha * theta_dt - mu_intercept - g_d, where
// g_d is built from the SAME gamma in both directed and undirected data:
//   - directed, d = (i -> j): g_d = gamma_i (the sender generates the
//     events, so its repertoire applies)
//   - undirected, d = {i, j}: g_d = w_send_d * gamma_i + (1 - w_send_d) *
//     gamma_j, w_send_d the observed share of that pair's events with i
//     as actor (data, not a parameter)
// One program serves both: w_send = 1 for directed data reduces the
// undirected formula to the directed one exactly. This is the point of
// the design -- both datasets estimate the same gamma, with the same
// prior and the same interpretation, and the difference in how it enters
// follows mechanically from what a dyad-year IS in each. The models also
// nest: as sigma_gamma -> 0, both this program and `stable` coincide
// (see tests/testthat/test_gamma_offset.R's exact-nesting test at
// n_countries = 1).
//
// IDENTIFICATION: GAMMA CONSTRAINTS. Three constraints on gamma, enforced
// by construction in transformed parameters below, in this order (the
// order is not arbitrary -- see the per-step comments there):
//   1. Orthogonal to alpha. Shifting gamma_c by k*alpha shifts every
//      dyad containing c by exactly what that dyad's own free theta can
//      already absorb -- the confounding is exact, since theta is free
//      per dyad-period. Projected out via dot_self(alpha), NOT via a
//      hardcoded /A: that division is correct only because alpha is
//      normalised to RMS 1 (so dot_self(alpha) == A) -- writing the
//      general form means a future change to alpha's normalisation
//      cannot silently break this projection.
//   2. Orthogonal to 1 (per country). Softmax shift invariance -- a
//      constant added to every category of gamma_c cancels in softmax(),
//      so it is unidentified; centred out per country.
//   3. Centred across countries (per category). The common,
//      country-independent part of any category's offset is absorbed by
//      mu_intercept; centred out per category.
// Verified by construction (not taken on faith) in
// tests/testthat/test_gamma_offset.R that steps 2-4 (row-centre,
// column-centre, alpha-projection) all hold SIMULTANEOUSLY after all four
// steps -- see that test for why each later step provably preserves the
// earlier ones' zero-mean property.
//
// gamma is NOT multiplied by orientation_sign(alpha_raw): the projection
// v - (dot_product(v, alpha) / dot_self(alpha)) * alpha is invariant to
// alpha -> -alpha (the sign cancels: it appears once in the numerator's
// alpha and once in the denominator's alpha, and once more in the
// dot_product's alpha... concretely, replacing alpha with -alpha leaves
// dot_product(v, alpha)/dot_self(alpha) unchanged, since both numerator
// and denominator are even in alpha's sign, and the subtracted term
// therefore uses -alpha, so the WHOLE expression v - (...)*(-alpha) using
// the negated coefficient equals the original -- the projection removes
// exactly the same component of v regardless of alpha's sign). So gamma,
// and hence g_d and eta, are already orientation-free -- unlike alpha and
// theta, which the fold explicitly flips. See
// tests/testthat/test_gamma_offset.R's orientation-invariance test.
//
// PRIOR: sigma_gamma ~ normal(0, 0.3) (half, via vector<lower=0>[A]),
// consistent with the dyad-level implied RMS of 0.217 being an UPPER
// bound on the country-level component (a country-level effect can only
// be part of what a dyad-level check sees, since a dyad-level residual
// pools whatever country-level and genuinely-dyad-level misfit exists).
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
functions {
  // <<< BEGIN GENERATED partial_log_lik (source: inst/stan/include/partial_log_lik.stanfunctions) >>>
  // Do not hand-edit between these markers -- edit the source file
  // above and rerun `Rscript data-raw/sync_stan_functions.R`
  // (checked by tests/testthat/test_stan_includes.R).
  // Shared reduce_sum likelihood for all bilatr Stan model variants.
  //
  // Indifferent to how `alpha` and `mu_intercept` were constructed upstream
  // (fixed-reference vs. sum-to-zero, static vs. OU-derived theta, ...): it
  // only consumes the already-built `theta`/`alpha`/`mu_intercept`/`phi`
  // vectors.
  //
  // This is the single canonical source. It is NOT included at compile/
  // sample time via Stan's `#include` (cmdstanr breaks `#include` resolution
  // at $sample()-time whenever include_paths contains a space -- see
  // https://github.com/stan-dev/cmdstanr/issues/820 -- which bites this
  // project's own devtools::load_all() working tree). Instead,
  // `data-raw/sync_stan_functions.R` splices this file's contents verbatim
  // into a marker-delimited block in each registered model's `.stan` file;
  // see `R/stan_includes.R`. Edit only this file, then rerun the sync
  // script -- do not hand-edit the generated blocks, and do not fork this
  // file per variant.
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

  // Country-offset variant of partial_log_lik(), added alongside it (0.7.0,
  // bilatr_alphanorm_gamma.stan / registered `stable_gamma`) -- not a fork,
  // per this file's own header: partial_log_lik() itself is untouched, and
  // every other registered program (stable, ou, and the legacy soft-anchor
  // pair) gets this spliced in too via data-raw/sync_stan_functions.R, as
  // an unused extra function, so their .stan text changes (one-time
  // recompile) but their data/parameter blocks -- and hence fits already on
  // disk -- are untouched.
  //
  // g_d, the per-dyad country offset, is built ONCE PER DYAD, outside the
  // `t` loop -- not stored anywhere (see bilatr_alphanorm_gamma.stan's
  // header, "Part 0b": a per-dyad-period or per-dyad matrix in transformed
  // parameters would be written to the output CSV at prohibitive cost).
  // w_send[d] == 1 covers directed data (g_d = gamma[, ctry_a[d]] exactly,
  // ctry_b[d] unused) and the directed-as-degenerate-undirected case alike;
  // any other value mixes ctry_a's and ctry_b's columns by event share (see
  // bilatr_alphanorm_gamma.stan's header for the directed/undirected
  // design this implements).
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
      // once per dyad, not per dyad-period
      vector[A] g = w_send[d] == 1.0
                    ? gamma[, ctry_a[d]]
                    : w_send[d] * gamma[, ctry_a[d]] + (1 - w_send[d]) * gamma[, ctry_b[d]];
      for (t in 1:T) {
        if (is_obs[d, t] == 1) {
          vector[A] eta = alpha .* rep_vector(theta[d, t], A) - mu_intercept - g;
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
  // bilatr_alphanorm.stan's header, "IDENTIFICATION: ORIENTATION FOLD".
  // Hand-duplicated here, same reasoning as bilatr_alphanorm_ou.stan's
  // own copy: tiny, not worth extending the sync tooling for.
  real orientation_sign(vector alpha_raw) {
    return alpha_raw[1] >= 0 ? 1.0 : -1.0;
  }

  // Per-dyad country offset g_d -- same construction as
  // partial_log_lik_offset() (inst/stan/include/
  // partial_log_lik.stanfunctions), duplicated here for generated
  // quantities' use (GQ-only, small, hand-verified identical -- same
  // convention as bilatr_alphanorm.stan's dyad_period_log_lik()).
  vector country_offset(matrix gamma, int ctry_a_d, int ctry_b_d, real w_send_d) {
    return w_send_d == 1.0
           ? gamma[, ctry_a_d]
           : w_send_d * gamma[, ctry_a_d] + (1 - w_send_d) * gamma[, ctry_b_d];
  }

  // Offset-aware per-dyad-period log-likelihood, not reduced/summed --
  // used only by the compute_log_lik generated quantities block below.
  // Must stay numerically identical to the per-cell term inside
  // partial_log_lik_offset() above (same eta/softmax/conc/
  // dirichlet_multinomial_lpmf, now with g_d precomputed by the caller
  // once per dyad -- see country_offset() above). Duplicated rather than
  // shared via the GENERATED mechanism for the same reason
  // dyad_period_log_lik() is in bilatr_alphanorm.stan: GQ-only, small.
  real dyad_period_log_lik_offset(int obs_dt, array[] int y_dt, real theta_dt,
                                   int A, vector mu_intercept, real phi_d,
                                   vector alpha, vector g_d) {
    if (obs_dt == 0) {
      return 0;
    }
    vector[A] eta = alpha .* rep_vector(theta_dt, A) - mu_intercept - g_d;
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
                                              // country count -- see Part 0a
                                              // of the build prompt; that is
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

  // Country-offset data (0.7.0; see header). n_countries is the number
  // of DISTINCT countries actually present in the retained dyads (see
  // R/stan_data.R's assemble_stan_data()) -- deliberately NOT named `C`
  // (already the reduce_sum grainsize above).
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
  real mu_log_noise;                // NOTE: now a ratio to sigma_theta0, see header
  real<lower=0> sigma_log_noise;
  vector<lower=0>[D] phi;
  real mu_log_phi;
  real<lower=0> sigma_log_phi;

  // Country-level category offsets (0.7.0; see header)
  matrix[A, n_countries] gamma_z;      // non-centred
  vector<lower=0>[A] sigma_gamma;      // per-category spread across countries
}
transformed parameters {
  // Normalize alpha to RMS (== population SD, sum-to-zero) 1. This is the
  // scale anchor for the whole model; see bilatr_alphanorm.stan's header
  // for the radial degeneracy this creates. Sign is NOT fixed here --
  // see orientation_sign() below and that header, "IDENTIFICATION:
  // ORIENTATION FOLD".
  vector[A] alpha = orientation_sign(alpha_raw)
                    * alpha_raw * sqrt((1.0 * A) / dot_self(alpha_raw));

  // Country-level category offsets. See header, "IDENTIFICATION: GAMMA
  // CONSTRAINTS" -- the four-step order matters: step 2 (row-centre
  // across countries) leaves step 3's per-column centring unable to
  // reintroduce a row mean (the constants step 3 subtracts average to
  // zero across countries, since each row already sums to zero before
  // step 3 runs); step 4's alpha-projection subtracts a multiple of
  // alpha, which itself sums to zero (built from a sum_to_zero_vector),
  // so it cannot perturb the column (per-country) sums step 3 already
  // zeroed, and its own per-country coefficients average to zero across
  // countries so it cannot reintroduce a row mean either. After all four
  // steps, constraints 1-3 hold EXACTLY, simultaneously -- see
  // tests/testthat/test_gamma_offset.R.
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
  // direction of alpha_raw (see bilatr_alphanorm.stan's header, "RADIAL
  // DEGENERACY"). The orientation fold (transformed parameters, above)
  // leaves this prior untouched.
  alpha_raw ~ std_normal();

  // Country-level category offsets (0.7.0; see header). sigma_gamma's
  // half-normal(0, 0.3) is calibrated against the dyad-level implied RMS
  // this model exists to test (0.217), as an upper bound on the
  // country-level component. gamma_z ~ std_normal() (via to_vector(),
  // element-wise): note gamma_z's projected-out directions (the alpha
  // direction, the all-ones direction, and the across-country mean) are
  // pinned ONLY by this prior, not by the likelihood -- see
  // R/diagnose_convergence.R's .bilatr_sign_ambiguous_raw_names()/
  // R/orient.R for where this is kept out of the tiered diagnostics
  // tables (not sign ambiguity, prior-only directions -- a different
  // reason, same tier = NA treatment).
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
      vector[A] g_d = country_offset(gamma, ctry_a[d], ctry_b[d], w_send[d]);
      for (t in 1:T) {
        log_lik[d, t] = dyad_period_log_lik_offset(
          is_obs[d, t], Y[d, t], theta[d, t], A, mu_intercept, phi[d], alpha, g_d
        );
      }
    }
  }

  // Forward-filtered theta: conditional on THIS draw's hyperparameters,
  // a Fisher-scoring (West-Harrison linear-Bayes) filter over the
  // observations -- an approximation to p(theta_t | y_1:t), not the exact
  // marginal. See bilatr_alphanorm.stan's header for the full derivation
  // of the score/information; the offset changes exactly one thing here
  // -- eta gains `- g_d` -- because d(eta_k)/d(theta) = alpha_k
  // regardless of g_d (g_d does not depend on theta), so the score and
  // information keep their form and are simply evaluated at
  // p = softmax(alpha * m_pred - mu_intercept - g_d). g_d is computed
  // once per filtered dyad, outside the t loop, same reasoning as
  // partial_log_lik_offset()'s per-dyad g.
  array[compute_theta_filtered ? n_filter_dyads : 0,
        compute_theta_filtered ? T : 0] real theta_filtered;
  array[compute_theta_filtered ? n_filter_dyads : 0,
        compute_theta_filtered ? T : 0] real theta_filtered_sd;

  if (compute_theta_filtered) {
    for (i in 1:n_filter_dyads) {
      int d = filter_dyads[i];
      vector[A] g_d = country_offset(gamma, ctry_a[d], ctry_b[d], w_send[d]);
      real m = 0;
      real p_var = square(sigma_theta0);
      for (t in 1:T) {
        real m_pred = m;
        real p_pred = p_var + square(process_noise[d]);
        if (is_obs[d, t] == 1) {
          vector[A] eta = alpha .* rep_vector(m_pred, A) - mu_intercept - g_d;
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
