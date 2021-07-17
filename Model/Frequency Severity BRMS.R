#### Load Packages ####

library(brms)
library(dplyr)
library(purrr)
library(tidyr)
library(renv)
library(stringr)
library(rstan)
library(rstanarm)

options(stringsAsFactors = FALSE,
        mc.cores = parallel::detectCores())

#### Simulate Frequency Data ####

#' Assuming one rating factor, region, with one group.
#' This is just a placeholder until the next stage where we adapt the model
#' to work with multiple rating factors.

set.seed(123456)

regions = c("EMEA", "USC")

freq_n = 5e3
freq_lambda_vec = c(EMEA = 2, USC = 3)

# Defines a non-linear function for lambda to test model still works

lambda_fun = function(expo) 1

freq_data =
  data.frame(
    pol_id =  seq(freq_n),
    freq = TRUE,
    sev = FALSE,
    expo = runif(freq_n, 1, 100),
    ded = runif(freq_n, 1e3, 10e3),
    lim = runif(freq_n, 25e3, 100e3),
    region = sample(regions, freq_n, replace = T)
  ) %>%
  mutate(
    freq_lambda = freq_lambda_vec[region] * lambda_fun(expo),
    claimcount_fgu = rpois(freq_n, freq_lambda),
    loss = 1
  )

#### Simulate severity Data ####

mu_fun = function(expo) 1

sev_mu_vec = c(EMEA = 8, USC = 9)
sev_sigma_vec = c(EMEA = 1, USC = 1.5)

sev_data =
  data.frame(
    ded = rep(freq_data$ded,
              freq_data$claimcount_fgu),
    lim = rep(freq_data$lim,
              freq_data$claimcount_fgu),
    region = rep(freq_data$region,
                 freq_data$claimcount_fgu),
    expo = rep(freq_data$expo,
                 freq_data$claimcount_fgu)
  ) %>%
  mutate(
    loss_uncapped =
      unlist(
        lapply(
          seq(freq_n),
          function(i){
            
            rlnorm(freq_data$claimcount_fgu[i], 
                   sev_mu_vec[freq_data$region[i]] *
                     mu_fun(freq_data$expo[i]), 
                   sev_sigma_vec[freq_data$region[i]])
            
          }
        )
      )
  ) %>%
  mutate(
    pol_id = rep(seq(freq_n), freq_data$claimcount_fgu)
  ) %>% 
  filter(
    loss_uncapped > ded
  ) %>%
  mutate(
    claim_id = row_number(),
    freq = FALSE,
    sev = TRUE,
    lim_exceed = as.integer(loss_uncapped >= lim),
    loss = pmin(loss_uncapped, lim),
    claimcount = 0,
    claimcount_fgu = 0
  )

freq_data_net =
  freq_data %>%
  left_join(
    sev_data %>%
      group_by(
        pol_id
      ) %>%
      summarise(
        claimcount = n()
      ) %>%
      ungroup(),
    by = "pol_id"
  ) %>%
  mutate(
    lim_exceed = 0,
    claimcount = coalesce(claimcount, 0L)
  )

#### Join Data ####

full_data = 
  bind_rows(freq_data_net, 
            sev_data)

#### Multivariate Model ####

fit_freq = 
  bf(claimcount | subset(freq) ~ f1,
     f1 ~ 1 + region,
     nl = TRUE) + 
  poisson(link = "log")

fit_sev = 
  bf(loss | subset(sev) + trunc(lb = ded) + cens(lim_exceed) ~ 
       log(s1) + 7,
     s1 ~ 1 + region,
     sigma ~ 1 + region,
     nl = TRUE
  ) + 
  lognormal()


replace_freq_block <- function(stan_code, 
                               sev_dist = "lnorm", 
                               resp_freq = "claimcount",
                               resp_sev = "loss",
                               freq_par = c("f1"),
                               sev_par = c("s1")) {
  
  cdf_fun =
    case_when(
      sev_dist == "lnorm"  ~ "lognormal_cdf",
      sev_dist == "gamma"  ~ "gamma_cdf",
      sev_dist == "norm"   ~ "normal_cdf",
      sev_dist == "pareto" ~ "pareto_cdf"
    )
  
  sev_par_name =
    case_when(
      sev_dist == "lnorm"  ~ c("mu", "sigma"),
      sev_dist == "gamma"  ~ c("alpha", "beta"),
      sev_dist == "norm"   ~ c("mu", "sigma"),
      sev_dist == "pareto" ~ c("y_min", "alpha")
    )
    
  freq_par_name = paste0("mu_", resp_freq) 
  
  raw_code <- 
    stan_code %>% 
    as.character()
  
  ### Modify the mu_claimcount line
  
  new_code <- 
    gsub(
      str_glue('\\| {freq_par_name}'),
      str_glue(
      '\\| {freq_par_name} + 
        log(1 - {cdf_fun}(
            ded, 
            
            X_{resp_freq}_{freq_par}[, 1:K_{resp_sev}_{sev_par}] * 
            b_{resp_sev}_{sev_par}, 
            
            exp(Intercept_{sev_par_name[2]}_{resp_sev} + 
            X_{resp_freq}_f1[, 2:K_{sev_par_name[2]}_{resp_sev}] * 
            b_{sev_par_name[2]}_{resp_sev}
            )
      )
      )'
      ),
      raw_code
    )
  
  return(new_code)
}

mv_model_formula = fit_sev + fit_freq + set_rescor(FALSE)

stanvars = c(
  stanvar(
    x = freq_data_net$ded,
    name = "ded"
  )
)

priors = c(prior(normal(0, 1),
                 class = b,
                 coef = Intercept,
                 resp = claimcount,
                 nlpar = f1),
           
           prior(normal(8, 1),
                 class = b,
                 coef = Intercept,
                 resp = loss,
                 nlpar = s1),
           
           prior(lognormal(0, 1),
                 class = Intercept,
                 dpar = sigma,
                 resp = loss)
           )

mv_model_code =
  make_stancode(
    mv_model_formula,
    data = full_data,
    prior = priors,
    stanvars = stanvars
  ) %>%
  replace_freq_block()

mv_model_data =
  make_standata(
    mv_model_formula,
    data = full_data,
    prior = priors,
    stanvars = stanvars
  )

mv_model_fit_stan =
  stan(
    model_code = mv_model_code,
    data = mv_model_data,
    chains = 1,
    iter = 1000,
    warmup = 250,
    control = 
      list(adapt_delta = 0.999,
           max_treedepth = 15)
    )

## Convert back to BRMS fit object

mv_model_fit <- 
  brm( formula = mv_model_formula,
       data = full_data, 
       prior = priors, 
       empty = TRUE
       )

mv_model_fit$fit <- mv_model_fit_stan

mv_model_fit <- rename_pars(mv_model_fit )

#### Results ####

model_post_samples =
  posterior_samples(
    mv_model_fit
  )

save(
  full_data,
  mv_model_fit,
  model_post_samples,
  file = "Model/mv_model_fit.RData"
)

base::load(file = "Model/mv_model_fit.RData")

sev_output =
  full_data %>%
  filter(sev) %>%
  mutate(
    mu_pred = 
      posterior_linpred(
        mv_model_fit,
        resp = "loss"
      ) %>%
      colMeans(),
    mu_pred_q025 = 
      posterior_linpred(
        mv_model_fit,
        resp = "loss"
      ) %>%
      apply(2, function(x) quantile(x, 0.025)),
    mu_pred_q975 = 
      posterior_linpred(
        mv_model_fit,
        resp = "loss"
      ) %>%
      apply(2, function(x) quantile(x, 0.975)),
    
    sigma_pred = 
      posterior_epred(
        mv_model_fit,
        resp = "loss",
        dpar = "sigma"
        ) %>%
      colMeans(),
    sigma_pred_q025 = 
      posterior_epred(
        mv_model_fit,
        resp = "loss",
        dpar = "sigma"
      ) %>%
      apply(2, function(x) quantile(x, 0.025)),
    sigma_pred_q975 = 
      posterior_epred(
        mv_model_fit,
        resp = "loss",
        dpar = "sigma"
      ) %>%
      apply(2, function(x) quantile(x, 0.975))
  )

# This is just a high-level check for reasonableness and does not reflect
# the "true" posterior of ground up lambdas predicted by the model

sev_output_mean =
  sev_output %>%
  group_by(
    region
  ) %>%
  summarise(
    mu_pred_mean = mean(mu_pred),
    sigma_pred_mean = mean(sigma_pred)
  ) %>%
  ungroup()

sev_output_mean =
  sev_output %>%
  group_by(
    region
  ) %>%
  summarise(
    mu_pred = mean(mu_pred),
    sigma_pred = mean(sigma_pred)
  ) %>%
  ungroup()

freq_output =
  freq_data_net %>%
  left_join(
    sev_output_mean,
    by = "region"
  ) %>%
  mutate(
    lambda_pred = 
      posterior_linpred(
        mv_model_fit,
        newdata = 
          full_data %>%
          left_join(
            sev_output_mean,
            by = "region"
          ) %>%
          mutate(ded = 0),
        resp = "claimcount"
        ) %>%
      colMeans() %>%
      exp() /
      (1 - plnorm(ded, mu_pred, sigma_pred)),
    
    lambda_pred_q025 = 
      posterior_linpred(
        mv_model_fit,
        newdata = 
          full_data %>%
          left_join(
            sev_output_mean,
            by = "region"
          ) %>%
          mutate(ded = 0),
        resp = "claimcount"
      ) %>%
      exp() %>%
      apply(2, function(x) quantile(x, 0.025)) /
      (1 - plnorm(ded, mu_pred, sigma_pred)),
    
    lambda_pred_q975 = 
      posterior_linpred(
        mv_model_fit,
        newdata = 
          full_data %>%
          left_join(
            sev_output_mean,
            by = "region"
          ) %>%
          mutate(ded = 0),
        resp = "claimcount"
      ) %>%
      exp() %>%
      apply(2, function(x) quantile(x, 0.975)) /
      (1 - plnorm(ded, mu_pred, sigma_pred))
  )

model_compare =
  freq_output %>%
  group_by(
    region
  ) %>%
  summarise_at(
    vars(
      lambda_pred,
      lambda_pred_q025,
      lambda_pred_q975
    ),
    mean
  ) %>%
  ungroup() %>%
  left_join(
    sev_output %>%
    group_by(
      region
    ) %>%
      summarise_at(
        vars(
          mu_pred,
          mu_pred_q025,
          mu_pred_q975,
          
          sigma_pred,
          sigma_pred_q025,
          sigma_pred_q975
        ),
        mean
      ) %>%
      ungroup(),
    by = "region"
  ) %>%
  mutate(
    mu_actual     = sev_mu_vec[region],
    sigma_actual  = sev_sigma_vec[region],
    lambda_actual = freq_lambda_vec[region]
  ) %>%
  select(sort(names(.)))
