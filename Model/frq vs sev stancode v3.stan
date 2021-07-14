// generated with brms 2.15.0
functions {
    /* poisson lognormal composite lpdf 
   * Args: 
   *   y: the loss amount
   *   y1: the claim count
   *   mu: mean parameter of the lognormal distribution 
   *   sigma: sd parameter of the lognormal distribution
   *   lambda: the mean of the poisson distribution
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
  real poisson_lognormal_lpdf(real y, int y1, real ded, real mu, real sigma, real lambda) {
    if (y == 0) {
      return poisson_lpmf(y1 | lambda * (1 - lognormal_cdf(ded, mu, sigma)));
    } else {
      return lognormal_lpdf(y | mu, sigma) - lognormal_lccdf(ded | mu, sigma);
    }
  }
  
}

data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_claimcount;  // number of observations
  int Y_claimcount[N_claimcount];   // response variable
  int<lower=1> K_claimcount_f1;  // number of population-level effects
  matrix[N_claimcount, K_claimcount_f1] X_claimcount_f1;  // population-level design matrix
  // covariate vectors for non-linear functions
  vector[N_claimcount] C_claimcount_1;
  int<lower=1> N_loss;  // number of observations
  vector[N_loss] Y_loss;  // response variable
  real lb_loss[N];
  int<lower=1> K_loss;  // number of population-level effects
  matrix[N_loss, K_loss] X_loss;  // population-level design matrix
  // covariate vectors for non-linear functions
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector[K_claimcount_f1] b_claimcount_f1;  // population-level effects
  vector[K_loss] b_loss;  // population-level effects
  real<lower=0> sigma_loss;  // residual SD
}
transformed parameters {
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N_claimcount] nlp_claimcount_f1 = X_claimcount_f1 * b_claimcount_f1;
    // initialize non-linear predictor term
    vector[N_claimcount] mu_claimcount;
    // initialize linear predictor term
    vector[N_loss] nlp_loss = X_loss * b_loss;
    // initialize non-linear predictor term
    vector[N_loss] mu_loss;
    for (n in 1:N) {
      // compute non-linear predictor values
      mu_loss[n] = nlp_loss[n];
      mu_claimcount[n] = nlp_claimcount_f1[n] + log(1 - lognormal_cdf(lb_loss[n], mu_loss[n], sigma_loss));
    }
    for (n in 1:N) {
      if (C_claimcount_1[n] == 1) {
        target += poisson_log_lpmf(Y_claimcount[n] | mu_claimcount[n]);
      } else {
          target += lognormal_lpdf(Y_loss[n] | mu_loss[n], sigma_loss) -
        lognormal_lccdf(lb_loss[n] | mu_loss[n], sigma_loss);
      }
    }
  }
  // priors including constants
  target += normal_lpdf(b_claimcount_f1 | 0, 1);
  target += normal_lpdf(b_loss | 5, 2);
  target += student_t_lpdf(sigma_loss | 3, 0, 3.9)
    - 1 * student_t_lccdf(0 | 3, 0, 3.9);
}
generated quantities {
}