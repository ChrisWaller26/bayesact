// generated with brms 2.15.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_claimcount;  // number of observations
  int Y_claimcount[N_claimcount];  // response variable
  int<lower=1> K_claimcount_f1;  // number of population-level effects
  matrix[N_claimcount, K_claimcount_f1] X_claimcount_f1;  // population-level design matrix
  // covariate vectors for non-linear functions
  vector[N_claimcount] C_claimcount_1;
  int<lower=1> N_loss;  // number of observations
  vector[N_loss] Y_loss;  // response variable
  int<lower=-1,upper=2> cens_loss[N_loss];  // indicates censoring
  real lb_loss[N_loss];  // lower truncation bounds;
  int<lower=1> K_loss;  // number of population-level effects
  matrix[N_loss, K_loss] X_loss;  // population-level design matrix
  int<lower=1> K_sigma_loss;  // number of population-level effects
  matrix[N_loss, K_sigma_loss] X_sigma_loss;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc_loss = K_loss - 1;
  matrix[N_loss, Kc_loss] Xc_loss;  // centered version of X_loss without an intercept
  vector[Kc_loss] means_X_loss;  // column means of X_loss before centering
  int Kc_sigma_loss = K_sigma_loss - 1;
  matrix[N_loss, Kc_sigma_loss] Xc_sigma_loss;  // centered version of X_sigma_loss without an intercept
  vector[Kc_sigma_loss] means_X_sigma_loss;  // column means of X_sigma_loss before centering
  for (i in 2:K_loss) {
    means_X_loss[i - 1] = mean(X_loss[, i]);
    Xc_loss[, i - 1] = X_loss[, i] - means_X_loss[i - 1];
  }
  for (i in 2:K_sigma_loss) {
    means_X_sigma_loss[i - 1] = mean(X_sigma_loss[, i]);
    Xc_sigma_loss[, i - 1] = X_sigma_loss[, i] - means_X_sigma_loss[i - 1];
  }
}
parameters {
  vector[K_claimcount_f1] b_claimcount_f1;  // population-level effects
  vector[Kc_loss] b_loss;  // population-level effects
  real Intercept_loss;  // temporary intercept for centered predictors
  vector[Kc_sigma_loss] b_sigma_loss;  // population-level effects
  real Intercept_sigma_loss;  // temporary intercept for centered predictors
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
    vector[N_loss] mu_loss = Intercept_loss + Xc_loss * b_loss;
    // initialize linear predictor term
    vector[N_loss] sigma_loss = Intercept_sigma_loss + Xc_sigma_loss * b_sigma_loss;
    for (n in 1:N_loss) {
      // apply the inverse link function
      sigma_loss[n] = exp(sigma_loss[n]);
    }
    for (n in 1:N_claimcount) {
      // compute non-linear predictor values
      mu_claimcount[n] = nlp_claimcount_f1[n] + log(1 - lognormal_cdf(lb_loss[n], mu_loss[n], sigma_loss[n]));
    }
    
    for (n in 1:N) {
      if (C_claimcount_1[n] == 1) {
          target += poisson_log_lpmf(Y_claimcount[n] | mu_claimcount[n]);
      } else {
        if (cens_loss[n] == 0) {
          target += lognormal_lpdf(Y_loss[n] | mu_loss[n], sigma_loss[n]) -
              lognormal_lccdf(lb_loss[n] | mu_loss[n], sigma_loss[n]);
        } else if (cens_loss[n] == 1) {
          target += lognormal_lccdf(Y_loss[n] | mu_loss[n], sigma_loss[n]) -
              lognormal_lccdf(lb_loss[n] | mu_loss[n], sigma_loss[n]);
        } else if (cens_loss[n] == -1) {
          target += lognormal_lcdf(Y_loss[n] | mu_loss[n], sigma_loss[n]) -
              lognormal_lccdf(lb_loss[n] | mu_loss[n], sigma_loss[n]);
        }
      }
    }
  }
  // priors including constants
  target += normal_lpdf(b_claimcount_f1[1] | 0, 2);
}

generated quantities {
  // actual population-level intercept
  real b_loss_Intercept = Intercept_loss - dot_product(means_X_loss, b_loss);
  // actual population-level intercept
  real b_sigma_loss_Intercept = Intercept_sigma_loss - dot_product(means_X_sigma_loss, b_sigma_loss);
}