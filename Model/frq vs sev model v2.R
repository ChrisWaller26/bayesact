library(dplyr)
library(brms)


source('model/data simulate.R')
set.seed(12345)
full_data <- data_sim_simple(n_pols = 5000, Intercept_lambda = 0.5, Intercept_mu = 8, 
                             b_lambda = c(0, 0.2), b_mu = c(0, -1, 1), sigma = 2, 
                             deductibles = c(0, 5e3, 10e3, 25e3, 5e4))

fit_freq = 
  bf(claimcount ~ f1 + freq2 * 0 ,
     f1 ~ freq_lvl,
     nl = TRUE) + 
  poisson()

fit_sev = 
  bf(loss | trunc(lb = ded) ~ mu_lvl) + 
  lognormal()

stan_data <- brms::make_standata(
  formula = fit_freq + fit_sev + set_rescor(FALSE),
  data = full_data, 
  prior = c(prior(normal(0, 2), class = 'b', coef = 'Intercept', resp = 'claimcount', nlpar = 'f1'))
)

library(rstan)
options(mc.cores = parallel::detectCores())
fit <- rstan::stan(file = 'Model/frq vs sev stancode v3.stan', 
                   data = stan_data,
                   chains = 4, iter = 2000, warmup = 1000, 
                   cores = 4
)
fit