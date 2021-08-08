// generated with brms 2.15.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_loss;  // number of observations
  vector[N_loss] Y_loss;  // response variable
  int<lower=-1,upper=2> cens_loss[N_loss];  // indicates censoring
  real lb_loss[N_loss];  // lower truncation bounds;
  int<lower=1> K_loss_s1;  // number of population-level effects
  matrix[N_loss, K_loss_s1] X_loss_s1;  // population-level design matrix
  // covariate vectors for non-linear functions
  vector[N_loss] C_loss_1;
  int<lower=1> K_sigma_loss;  // number of population-level effects
  matrix[N_loss, K_sigma_loss] X_sigma_loss;  // population-level design matrix
  int<lower=1> N_claimcount;  // number of observations
  int Y_claimcount[N_claimcount];  // response variable
  int<lower=1> K_claimcount_f1;  // number of population-level effects
  matrix[N_claimcount, K_claimcount_f1] X_claimcount_f1;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
  vector[10417] ded;
  vector[10417] sev;
  vector[10417] freq;
}
transformed data {
  int Kc_sigma_loss = K_sigma_loss - 1;
  matrix[N_loss, Kc_sigma_loss] Xc_sigma_loss;  // centered version of X_sigma_loss without an intercept
  vector[Kc_sigma_loss] means_X_sigma_loss;  // column means of X_sigma_loss before centering
  for (i in 2:K_sigma_loss) {
    means_X_sigma_loss[i - 1] = mean(X_sigma_loss[, i]);
    Xc_sigma_loss[, i - 1] = X_sigma_loss[, i] - means_X_sigma_loss[i - 1];
  }
}
parameters {
  vector[K_loss_s1] b_loss_s1;  // population-level effects
  vector[Kc_sigma_loss] b_sigma_loss;  // population-level effects
  real Intercept_sigma_loss;  // temporary intercept for centered predictors
  vector[K_claimcount_f1] b_claimcount_f1;  // population-level effects
}
transformed parameters {
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N_loss] nlp_loss_s1 = X_loss_s1 * b_loss_s1;
    // initialize non-linear predictor term
    vector[N_loss] mu_loss;
    // initialize linear predictor term
    vector[N_loss] sigma_loss = Intercept_sigma_loss + Xc_sigma_loss * b_sigma_loss;
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
      mu_loss[n] = nlp_loss_s1[n] * C_loss_1[n] / C_loss_1[n];
    }
    for (n in 1:N_claimcount) {
      // compute non-linear predictor values
      mu_claimcount[n] = nlp_claimcount_f1[n];
    }
    for (n in 1:N_loss) {
    // special treatment of censored data
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
    target += poisson_log_lpmf(Y_claimcount | mu_claimcount);
  }
  // priors including constants
  target += normal_lpdf(b_loss_s1[1] | 8, 1);
  target += lognormal_lpdf(Intercept_sigma_loss | 0, 1);
  target += normal_lpdf(b_claimcount_f1[1] | 0, 1);
}
generated quantities {
  // actual population-level intercept
  real b_sigma_loss_Intercept = Intercept_sigma_loss - dot_product(means_X_sigma_loss, b_sigma_loss);
  b_sigma_loss_Intercept = b_sigma_loss_Intercept + dot_product(means_X_sigma_loss, b_sigma_loss);
}

