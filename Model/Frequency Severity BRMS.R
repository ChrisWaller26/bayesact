#### Load Packages ####

library(brms)
library(dplyr)
library(purrr)
library(tidyr)
library(renv)

options(stringsAsFactors = FALSE)

#### Simulate Frequency Data ####

#' Assuming one rating factor, region, with one group.
#' This is just a placeholder until the next stage where we adapt the model
#' to work with multiple rating factors.

regions = c("EMEA")

freq_n = 25e3
freq_lambda_vec = c(EMEA = 2)

freq_data =
  data.frame(
    pol_id =  seq(freq_n),
    freq = TRUE,
    sev = FALSE,
    ded = runif(freq_n, 1e3, 10e3),
    region = sample(regions, freq_n, replace = T)
  ) %>%
  mutate(
    freq_lambda = freq_lambda_vec[region],
    claimcount_fgu = rpois(freq_n, freq_lambda),
    loss = 1
  )

#### Simulate severity Data ####

sev_mu_vec = c(EMEA = 8)
sev_sigma_vec = c(EMEA = 1.5)

sev_data =
  data.frame(
    ded = rep(freq_data$ded,
              freq_data$claimcount_fgu),
    region = rep(freq_data$region,
              freq_data$claimcount_fgu)
  ) %>%
  mutate(
    loss =
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
    loss > ded
  ) %>%
  mutate(
    claim_id = row_number(),
    freq = FALSE,
    sev = TRUE,
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

stanvars =
  "
  nlp_claimcount_f1  = 
    nlp_claimcount_f1 + 
      log(1 - lognormal_cdf(
              ded, 
              Intercept_loss, 
              Intercept_sigma_loss)
              );
  "


fit_freq = 
  bf(claimcount | subset(freq) ~ f1,
     f1 ~ 1 + region,
     nl = TRUE) + 
  poisson()

fit_sev = 
  bf(loss | subset(sev) + trunc(lb = ded) ~ 1,
     sigma ~ 1
  ) + 
  lognormal()


mv_model_fit =
  brm(fit_freq + fit_sev + set_rescor(FALSE),
      data = full_data,
      stanvars =
        c(stanvar(
          scode = stanvars,
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
      control = list(adapt_delta = 0.99,
                     max_treedepth = 10)
  )

#### Results ####


model_pred_sev_mu =
  posterior_linpred(
    mv_model_fit,
    resp = "loss"
  )

model_pred_sev_sigma =
  posterior_epred(
    mv_model_fit,
    resp = "loss",
    dpar = "sigma"
  )

model_pred_freq =
  posterior_epred(
    mv_model_fit,
    resp = "claimcount"
  )

model_compare =
  freq_data_net %>%
  mutate(
    lambda = freq_lambda,
    pred_lambda   = colMeans(model_pred_freq) /
      (1 - plnorm(ded,
                  mean(model_pred_sev_mu),
                  mean(model_pred_sev_sigma))
      )
  ) %>%
  summarise(
    pred_lambda = mean(pred_lambda)
  ) %>%
  ungroup() %>%
  bind_cols(
    sev_data %>%
      mutate(
        mu    = colMeans(model_pred_sev_mu),
        sigma = colMeans(model_pred_sev_sigma)
      ) %>%
      summarise(
        pred_mu    = mean(mu),
        pred_sigma = mean(sigma)
      ) %>%
      ungroup()
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "par",
    values_to = "pred"
  ) %>%
  mutate(
    actual = c(freq_lambda_vec, sev_mu_vec, sev_sigma_vec),
    diff = paste0(round(100 * (1 - pred / actual), 1), "%"),
    pred = round(pred, 3)
  )
