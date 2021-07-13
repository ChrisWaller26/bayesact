#### Load Packages ####

library(brms)
library(dplyr)
library(purrr)
library(tidyr)
library(renv)

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

freq_data =
  data.frame(
    pol_id =  seq(freq_n),
    freq = TRUE,
    sev = FALSE,
    ded = runif(freq_n, 1e3, 10e3),
    lim = runif(freq_n, 25e3, 50e3),
    region = sample(regions, freq_n, replace = T)
  ) %>%
  mutate(
    freq_lambda = freq_lambda_vec[region],
    claimcount_fgu = rpois(freq_n, freq_lambda),
    loss = 1
  )

#### Simulate severity Data ####

sev_mu_vec = c(EMEA = 8, USC = 9)
sev_sigma_vec = c(EMEA = 1, USC = 1.5)

sev_data =
  data.frame(
    ded = rep(freq_data$ded,
              freq_data$claimcount_fgu),
    lim = rep(freq_data$lim,
              freq_data$claimcount_fgu),
    region = rep(freq_data$region,
              freq_data$claimcount_fgu)
  ) %>%
  mutate(
    loss_uncapped =
      unlist(
        lapply(
          seq(freq_n),
          function(i){
            
            rlnorm(freq_data$claimcount_fgu[i], 
                   sev_mu_vec[freq_data$region[i]], 
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
    claimcount = coalesce(claimcount, 0L)
  )

#### Join Data ####

full_data = 
  bind_rows(freq_data_net, 
            sev_data)

#### Multivariate Model ####

stanvars_lik =
  "
  nlp_claimcount_f1  = 
    nlp_claimcount_f1 + 
      log(1 - lognormal_cdf(
              ded, 
              Intercept_loss + X_claimcount_f1[, 2:K_loss] * b_loss, 
              exp(Intercept_sigma_loss + X_claimcount_f1[, 2:K_loss] * b_sigma_loss)
              )
              );
  "

fit_freq = 
  bf(claimcount | subset(freq) ~ f1 ,
     f1 ~ 1 + region,
     nl = TRUE) + 
  poisson()

fit_sev = 
  bf(loss | subset(sev) + trunc(lb = ded) + cens(lim_exceed) ~ 1 + region,
     sigma ~ 1 + region
  ) + 
  lognormal()


mv_model_fit =
  brm(fit_freq + fit_sev + set_rescor(FALSE),
      data = full_data,
      stanvars =
        c(stanvar(
          scode = stanvars_lik,
          block = "likelihood"
        ),
        stanvar(
          x = freq_data_net$ded,
          name = "ded"
        )
        ),
      prior = 
        c(prior(normal(1, 1),
                class = b,
                coef = Intercept,
                resp = claimcount,
                nlpar = f1),
          
          prior(normal(8, 2),
                class = Intercept,
                resp = loss),
          
          prior(lognormal(0, 1),
                class = Intercept,
                dpar = sigma,
                resp = loss)
        ),
      chains = 1,
      iter = 1000,
      warmup = 250,
      control = list(adapt_delta = 0.999,
                     max_treedepth = 10)
  )

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
    Intercept_loss = mean(mu_pred),
    Intercept_sigma_loss = mean(sigma_pred)
  ) %>%
  ungroup()

freq_output =
  freq_data_net %>%
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
      colMeans(),
    
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
      apply(2, function(x) quantile(x, 0.025)),
    
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
      apply(2, function(x) quantile(x, 0.975))
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
