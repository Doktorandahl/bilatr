// bilatr: hierarchical dynamic IRT model for dyadic conflict trajectories.
//
// EXPERIMENTAL variant `ou` (registered in R/model_registry.R,
// status = "experimental"). Not the default model; fit only via
// fit_dyad_ts_dev()/fit_panel_dev(stan_model = "ou").
//
// Motivation: the stable model's theta is a pure random walk, which has
// no restoring force. Between-dyad SD at t=1 is ~0.22 in the current
// stable fit, while expected accumulated drift over a ~16-year panel at
// process_noise ~ 0.2 is ~0.8. Any dyad can therefore traverse the entire
// cross-dyad range within a decade, so initial ordering carries almost no
// information about ordering a decade later -- this is why dyads that
// should be clearly separated in level (e.g. a low-conflict dyad vs a
// long-running rivalry) can invert during a temporary detente under
// stable.
//
// This variant replaces the random walk with an OU/AR(1) process with a
// dyad-specific equilibrium (mu_dyad) and a single GLOBAL persistence
// parameter rho < 1, so cross-dyad ordering has a permanent component
// (the equilibria) and a temporary one (deviations around it), and a
// detente reads as a bounded excursion rather than a new permanent level.
//
// IDENTIFICATION IS DELIBERATELY LEFT UNCHANGED FROM STABLE (alpha[1] = 1,
// mu_intercept[1] = 0, mu_theta_bar/sigma_mu playing exactly the role
// stable's mu_theta0/sigma_theta0 played), so that the comparison against
// stable isolates the effect of the dynamics alone. The affine-ridge fix
// (variant `alphanorm`) is a separate, orthogonal change; variant
// `alphanorm_ou` combines both. Do NOT "fix" the identification here.
//
// Design points that are load-bearing, not stylistic:
//   - The STATIONARY sd (sd_stat) is what gets a hierarchical prior; the
//     innovation sd (process_noise) is DERIVED from it
//     (process_noise = sd_stat * sqrt(1 - rho^2)), not the reverse. If
//     the prior sat on the innovation sd instead, the implied stationary
//     sd would be innovation_sd / sqrt(1 - rho^2), which at rho = 0.9
//     more than doubles it -- so the prior would mean something
//     different at every value of rho. Parameterizing sd_stat instead
//     makes the between/within ratio (sigma_mu vs exp(mu_log_sd_stat))
//     directly priorable, and keeps the model numerically well-behaved
//     as rho -> 1.
//   - theta[d, 1] is initialized from the STATIONARY distribution
//     (sd_stat[d], not process_noise[d]). Without this the first period
//     would be systematically under-dispersed relative to the rest of
//     the series (process_noise's marginal sd is smaller than the
//     stationary sd whenever rho > 0).
//   - rho is GLOBAL, not per-dyad: with T ~ 35 a per-dyad rho would be
//     weakly identified. A per-dyad rho is a possible future extension,
//     not implemented here.
//
// PRIOR UNITS: mu_theta_bar/sigma_mu play stable's mu_theta0/sigma_theta0
// role exactly (same units, same interpretation), so those two priors are
// a genuine port, not a reinterpretation. mu_log_sd_stat/sigma_log_sd_stat
// (the hierarchy on the per-dyad stationary sd) and the rho prior
// (rho_prior_a/rho_prior_b, default Beta(8, 2), i.e. weighted toward
// strong persistence) are NEW and have NOT been checked against prior
// predictive simulation.
//
// GENERATED QUANTITIES: sd_stat and mu_dyad are ordinary top-level
// transformed parameters, so they are already saved on every posterior
// draw without needing to be re-exposed here. generated quantities below
// adds within_between_ratio = sigma_mu / mean(sd_stat) (the single number
// this model's cross-dyad-comparability story reduces to) unconditionally
// (cheap, O(D)), plus log_lik gated behind compute_log_lik (D x T x
// draws, expensive -- default off).
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
                        vector alpha,
                        vector dyad_weight,
                        vector period_weight,
                        vector action_weight) {
    real lp = 0;
    for (d in start:end) {
      for (t in 1:T) {
        if (is_obs[d, t] == 1) {
          vector[A] eta = alpha .* rep_vector(theta[d, t], A) - mu_intercept;
          vector[A] p = softmax(eta);
          // action_weight rescales concentration per action type
          vector[A] conc = phi[d] * (action_weight .* p);
          // dyad_weight/period_weight scale each dyad-period's contribution
          lp += dyad_weight[d] * period_weight[t] *
            dirichlet_multinomial_lpmf(Y[d, t] | conc);
        }
      }
    }
    return lp;
  }
  // <<< END GENERATED partial_log_lik >>>

  // Per-dyad-period log-likelihood, not reduced/summed -- used only by the
  // compute_log_lik generated quantities block below. Must stay
  // numerically identical to the per-cell term inside partial_log_lik
  // above (same eta/softmax/conc/dirichlet_multinomial_lpmf). Duplicated
  // rather than shared via the GENERATED mechanism because it is
  // GQ-only, needed by the three experimental variants only (not
  // stable), and small enough that hand-verified identity across the
  // three files is lower-risk than extending the sync tooling for it.
  real dyad_period_log_lik(int obs_dt, array[] int y_dt, real theta_dt,
                            int A, vector mu_intercept, real phi_d,
                            vector alpha, real dyad_weight_d,
                            real period_weight_t, vector action_weight) {
    if (obs_dt == 0) {
      return 0;
    }
    vector[A] eta = alpha .* rep_vector(theta_dt, A) - mu_intercept;
    vector[A] p = softmax(eta);
    vector[A] conc = phi_d * (action_weight .* p);
    return dyad_weight_d * period_weight_t * dirichlet_multinomial_lpmf(y_dt | conc);
  }
}
data {
  int<lower=1> T;                            // number of time points
  int<lower=1> D;                            // number of dyads
  int<lower=1> A;                            // number of action types
  int<lower=1> C;                            // reduce_sum grainsize
  array[D, T] int<lower=0, upper=1> is_obs;  // observed indicator
  array[D, T, A] int<lower=0> Y;             // event counts
  vector<lower=0>[D] dyad_weight;            // per-dyad reweighting, default 1s
  vector<lower=0>[T] period_weight;          // per-period reweighting, default 1s
  vector<lower=0>[A] action_weight;          // per-action-type reweighting, default 1s
  int<lower=0, upper=1> compute_log_lik;     // 1 = also compute per-dyad-period
                                              // log_lik in generated quantities
                                              // (D x T x draws; default 0/off)
  real<lower=0> rho_prior_a;                 // Beta(rho_prior_a, rho_prior_b) on
  real<lower=0> rho_prior_b;                 // rho; default 8, 2 (weighted toward
                                              // strong persistence)
}
parameters {
  array[D, T] real theta_raw;
  vector[A - 1] mu_intercept_raw;   // global action baseline (index 2:A); identification unchanged from stable
  vector[A - 1] alpha_raw;          // alpha[1] = 1, rest free; identification unchanged from stable

  real mu_theta_bar;                 // takes the role of stable's mu_theta0
  real<lower=0> sigma_mu;            // takes the role of stable's sigma_theta0
  vector[D] mu_dyad_raw;
  real<lower=0, upper=1> rho;        // global OU/AR(1) persistence

  real mu_log_sd_stat;                 // hierarchy on the STATIONARY sd (not the
  real<lower=0> sigma_log_sd_stat;     // innovation sd -- see header)
  vector[D] log_sd_stat_raw;

  vector<lower=0>[D] phi;
  real mu_log_phi;
  real<lower=0> sigma_log_phi;
}
transformed parameters {
  vector[A] alpha;
  alpha[1] = 1;
  for (i in 2:A) {
    alpha[i] = alpha_raw[i - 1];
  }

  vector[A] mu_intercept;
  mu_intercept[1] = 0;
  for (i in 2:A) {
    mu_intercept[i] = mu_intercept_raw[i - 1];
  }

  vector[D] mu_dyad = mu_theta_bar + sigma_mu * mu_dyad_raw;

  // Stationary sd per dyad (hierarchical); innovation sd derived from it
  // and rho, not the reverse -- see header.
  vector<lower=0>[D] sd_stat =
    exp(mu_log_sd_stat + sigma_log_sd_stat * log_sd_stat_raw);
  vector<lower=0>[D] process_noise = sd_stat * sqrt(1 - square(rho));

  array[D, T] real theta;
  for (d in 1:D) {
    // stationary initialization at t = 1 (sd_stat, not process_noise)
    theta[d, 1] = mu_dyad[d] + sd_stat[d] * theta_raw[d, 1];
    for (t in 2:T) {
      theta[d, t] = mu_dyad[d] + rho * (theta[d, t - 1] - mu_dyad[d])
                    + process_noise[d] * theta_raw[d, t];
    }
  }
}
model {
  rho ~ beta(rho_prior_a, rho_prior_b);
  mu_theta_bar ~ normal(0, 1);
  sigma_mu ~ normal(0, 0.5);
  mu_dyad_raw ~ std_normal();

  mu_log_sd_stat ~ normal(log(0.5), 0.5);
  sigma_log_sd_stat ~ normal(0, 0.5);
  log_sd_stat_raw ~ std_normal();

  mu_log_phi ~ normal(0, 1);
  sigma_log_phi ~ normal(0, 1);
  phi ~ lognormal(mu_log_phi, sigma_log_phi);

  // theta OU/AR(1) innovations
  for (d in 1:D) {
    for (t in 1:T) {
      theta_raw[d, t] ~ std_normal();
    }
  }

  mu_intercept_raw ~ std_normal();
  alpha_raw ~ std_normal();

  // likelihood, chunked via reduce_sum
  array[D] int dyad_seq = linspaced_int_array(D, 1, D);
  target += reduce_sum(partial_log_lik, dyad_seq, C,
                        T, A, is_obs, Y, theta, mu_intercept, phi, alpha,
                        dyad_weight, period_weight, action_weight);
}
generated quantities {
  // the single number the cross-dyad-comparability problem reduces to
  real within_between_ratio = sigma_mu / mean(sd_stat);

  array[compute_log_lik ? D : 0, compute_log_lik ? T : 0] real log_lik;

  if (compute_log_lik) {
    for (d in 1:D) {
      for (t in 1:T) {
        log_lik[d, t] = dyad_period_log_lik(
          is_obs[d, t], Y[d, t], theta[d, t], A, mu_intercept, phi[d], alpha,
          dyad_weight[d], period_weight[t], action_weight
        );
      }
    }
  }
}
