// generated with brms 2.15.0
// Modified by Chris Waller 2021-07-22
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_loss;  // number of observations
  vector[N_loss] Y_loss;  // response variable
  int<lower=-1,upper=2> cens_loss[N_loss];  // indicates censoring
  real lb_loss[N_loss];  // lower truncation bounds;
  int<lower=1> K_loss_s1;  // number of population-level effects
  matrix[N_loss, K_loss_s1] X_loss_s1;  // population-level design matrix
  int<lower=1> K_sigma_loss;  // number of population-level effects
  matrix[N_loss, K_sigma_loss] X_sigma_loss;  // population-level design matrix
  int<lower=1> N_claimcount;  // number of observations
  int Y_claimcount[N_claimcount];  // response variable
  int<lower=1> K_claimcount_f1;  // number of population-level effects
  matrix[N_claimcount, K_claimcount_f1] X_claimcount_f1;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
  vector[N] ded;
  vector[N] sev;
  vector[N] freq;
}

parameters {
  vector[K_loss_s1] b_loss_s1;  // population-level effects
  vector[K_sigma_loss] b_sigma_loss;  // population-level effects
  vector[K_claimcount_f1] b_claimcount_f1;  // population-level effects
}

model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N_loss] nlp_loss_s1 = X_loss_s1 * b_loss_s1;
    // initialize non-linear predictor term
    vector[N_loss] mu_loss;
    // initialize linear predictor term
    vector[N_loss] sigma_loss = X_sigma_loss * b_sigma_loss;

    vector[N_loss] target_loss;
    // initialize linear predictor term
    vector[N_claimcount] nlp_claimcount_f1 = X_claimcount_f1 * b_claimcount_f1;
    // initialize non-linear predictor term
    vector[N_claimcount] mu_claimcount;
    for (n in 1:N_loss) {
      // apply the inverse link function
      sigma_loss[n] = exp(sigma_loss[n]);
    }
    for (n in 1:N_loss) {
      // compute non-linear predictor values
      mu_loss[n] = nlp_loss_s1[n];
    }
    for (n in 1:N_claimcount) {
      // compute non-linear predictor values
      mu_claimcount[n] = nlp_claimcount_f1[n];
    }
    for (n in 1:N_loss) {
    // special treatment of censored data
      if (cens_loss[n] == 0) {
        target_loss[n] = lognormal_lpdf(Y_loss[n] | mu_loss[n], sigma_loss[n]) -
        lognormal_lccdf(lb_loss[n] | mu_loss[n], sigma_loss[n]);
      } else if (cens_loss[n] == 1) {
        target_loss[n] = lognormal_lccdf(Y_loss[n] | mu_loss[n], sigma_loss[n]) -
        lognormal_lccdf(lb_loss[n] | mu_loss[n], sigma_loss[n]);
      } else if (cens_loss[n] == -1) {
        target_loss[n] = lognormal_lcdf(Y_loss[n] | mu_loss[n], sigma_loss[n]) -
        lognormal_lccdf(lb_loss[n] | mu_loss[n], sigma_loss[n]);
      }
    }
    
    target_loss = target_loss .* sev;
    
    target += sum(target_loss);
    target += freq * poisson_log_lpmf(Y_claimcount | mu_claimcount + 
  log(1 - lognormal_cdf(
      ded, 
      mu_loss,
      sigma_loss
      )
));

  }
  // priors including constants
  target += normal_lpdf(b_loss_s1[1] | 8, 1);
  target += normal_lpdf(b_sigma_loss | 0, 1);
  target += normal_lpdf(b_claimcount_f1[1] | 0, 1);
}

