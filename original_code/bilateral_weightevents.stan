functions {
  real partial_log_lik(array[] int slice_d,
                       int start, int end,
                       int T, int A,
                       array[,] int is_obs,
                       array[,,] int Y,
                       array[,] real theta,
                       vector mu_intercept,
                       vector phi,
                       vector alpha,
                       vector action_weight) {      // NEW
    real lp = 0;
    for (d in start:end) {
      for (t in 1:T) {
        if (is_obs[d, t] == 1) {
          vector[A] eta = alpha .* rep_vector(theta[d,t], A) - mu_intercept;
          vector[A] p = softmax(eta);
          // action_weight rescales concentration per action type
          // obs_weight downweights high-count dyad-periods
          vector[A] conc = phi[d] * (action_weight .* p);
          lp += dirichlet_multinomial_lpmf(Y[d,t] | conc);
        }
      }
    }
    return lp;
  }
}
data {
  int<lower=1> T;                 // Number of time points (days)
  int<lower=1> D;                 // Number of dyads
  int<lower=1> A;                 // Number of action types
  int<lower=1> C;                 // Chunk size
  array[D, T] int<lower=0, upper=1> is_obs; // Indicator for observed data
  array[D, T, A] int<lower=0> Y;  // Count event indicators
  vector<lower=0>[A] action_weight;        // action type weights
}

parameters {
  // Latent states per dyad
array[D, T] real theta_raw;
vector[A-1] mu_intercept_raw;  // global action baseline
real mu_theta0;
real<lower=0> sigma_theta0;
vector[D] z_theta0; 

  // Hierarchial parameters
  vector<lower=0>[D] process_noise;    // Process noise for dynamic trait
real mu_log_noise;
real<lower=0> sigma_log_noise;
vector<lower=0>[D] phi;
  real mu_log_phi;
  real<lower=0> sigma_log_phi;


  // Global parameters
  vector[A-2] alpha_raw;
  real<lower=0> alpha_hostile;      // known positive, used as negative

  
}

transformed parameters {

  vector[A] alpha;
  alpha[1] = 1;
  alpha[A] = -alpha_hostile;

for(i in 2:(A-1)) {
  alpha[i] = alpha_raw[i-1];
}

vector[A] mu_intercept;
mu_intercept[1] = 0;
for (i in 2:A) {
  mu_intercept[i] = mu_intercept_raw[i-1];
}

 vector[D] theta0;
  for (d in 1:D) {
    theta0[d] = mu_theta0 + sigma_theta0 * z_theta0[d];
  }

array[D, T] real theta;
 for(d in 1:D) {
  theta[d, 1] = theta0[d] + process_noise[d]*theta_raw[d,1];
  for (t in 2:T) {
    theta[d,t] = theta[d,t-1] + process_noise[d] * theta_raw[d,t];
  }
}

}

model {
  // hyperpriors for pooling (lognormal hierarchy)
  mu_log_noise ~ normal(log(0.2), 0.5);
  sigma_log_noise ~ normal(0, 0.5);  
  process_noise ~ lognormal(mu_log_noise, sigma_log_noise);

  mu_log_phi ~ normal(0, 1);
  sigma_log_phi ~ normal(0, 1);
  phi ~ lognormal(mu_log_phi, sigma_log_phi);

  // Priors for theta raw
for (d in 1:D)
  for (t in 1:T)
    theta_raw[d, t] ~ std_normal();


mu_theta0 ~ normal(0, 1);
sigma_theta0 ~ normal(0, 0.5);       // half-normal, constrained positive
  z_theta0 ~ std_normal();      // 

mu_intercept_raw ~ std_normal();


// Priors for global parameters
  alpha_raw ~ std_normal();
  alpha_hostile ~ normal(0, 1);

  // Likelihood
array[D] int dyad_seq = linspaced_int_array(D, 1, D);
  target += reduce_sum(partial_log_lik, dyad_seq, C,
                     T, A, is_obs, Y, theta, mu_intercept, phi, alpha,
                     action_weight);

}
