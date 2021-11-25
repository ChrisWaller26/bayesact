#### Load Packages ####

library(brms)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(rstan)
library(rstanarm)
library(readr)

# Package created in this repo

library(bayesact)

options(stringsAsFactors = FALSE,
        mc.cores = parallel::detectCores())

#### Simulate Frequency Data ####

#' Assuming one rating factor, region, with one group.
#' This is just a placeholder until the next stage where we adapt the model
#' to work with multiple rating factors.

set.seed(123456)

regions = c("EMEA", "USC")

freq_n = 5e3

# Defines a non-linear function for lambda to test model still works

freq_mu_fun = function(expo, region){
  exp(c(EMEA = 1, USC = 1.4)[region])
}

freq_data =
  data.frame(
    pol_id =  seq(freq_n),
    expo = runif(freq_n, 1, 100),
    ded = runif(freq_n, 1e3, 5e3),
    lim = runif(freq_n, 50e3, 100e3),
    region = sample(regions, freq_n, replace = T)
  ) %>%
  mutate(
    freq_mu = freq_mu_fun(expo, region),
    claimcount_fgu =
      rpois(freq_n, freq_mu)
  )

#### Simulate severity Data ####

mu_fun = function(expo, region){
  c(EMEA = 8, USC = 9)[region]
}

sigma_vec = exp(c(EMEA = 0, USC = 0))

alpha_vec = expp1(c(EMEA = -1, USC = -1))

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

            rlnormpower(
              freq_data$claimcount_fgu[i],
                   mu_fun(freq_data$expo[i],
                            freq_data$region[i]),
              sigma_vec[freq_data$region[i]],
              alpha_vec[freq_data$region[i]]
              )

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
    lim_exceed = as.integer(loss_uncapped >= lim),
    loss = pmin(loss_uncapped, lim)
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
    claimcount = coalesce(claimcount, 0)
  )

#### Run Model ####

mv_model_fit =
  brms_freq_sev(

    freq_formula =
      bf(claimcount ~ 1 + region),

    sev_formula =
      bf(loss | trunc(lb = ded) + cens(lim_exceed) ~
           1 + region,
         sigma ~ 1,
         alpha ~ 1
      ),

    freq_family = poisson(),
    sev_family = bayesact::lnormpower(),

    freq_data = freq_data_net,
    sev_data = sev_data,

    prior = c(prior(normal(0, 1),
                     class = Intercept,
                     resp = claimcount),
               prior(normal(0, 1),
                     class = b,
                     resp = claimcount),

              prior(normal(8, 1),
                    class = Intercept,
                    resp = loss),
              prior(normal(0, 1),
                    class = b,
                    resp = loss),

              prior(normal(0, 1),
                    class = Intercept,
                    resp = loss,
                    dpar = sigma),

              prior(normal(-1, 1),
                    class = Intercept,
                    resp = loss,
                    dpar = alpha)
    ),

    ded_name = "ded",
    ded_adj_min = 0.0001,
    use_cmdstan = F,

    chains = 1,
    iter = 1000,
    warmup = 500,

    refresh = 100,

    mle = FALSE,
    sample_prior = "no",
    freq_adj_fun = NULL,
    stanvars     = lnormpower_stanvars
  )

#### Results ####

model_post_samples =
  posterior_samples(
    mv_model_fit
  ) %>%
  transmute(
    s1_emea = b_loss_s1_Intercept,
    s1_usc  = b_loss_s1_Intercept +
      b_loss_s1_regionUSC,

    sigma = exp(b_sigma_loss_Intercept),

    alpha = exp(b_alpha_loss_Intercept) - 1,

    f1_emea = exp(b_claimcount_f1_Intercept),
    f1_usc  = exp(b_claimcount_f1_Intercept +
                    b_claimcount_f1_regionUSC)
  )

model_output =
  model_post_samples %>%
  sapply(
    function(x) c(lower = quantile(x, 0.025),
                  mean  = mean(x),
                  upper = quantile(x, 0.975))
  ) %>%
  as.data.frame() %>%
  bind_rows(
    data.frame(
      s1_emea = mu_fun(1, "EMEA"),
      s1_usc  = mu_fun(1, "USC"),

      sigma = sigma_vec["EMEA"],

      alpha = alpha_vec["EMEA"],

      f1_emea = freq_mu_fun(1, "EMEA"),
      f1_usc  = freq_mu_fun(1, "USC")
    )
  )

rownames(model_output) =
  c("Lower 2.5%",
    "Mean",
    "Upper 97.5%",
    "Actual")

# Check test has passed

(model_output[4,] > model_output[1,]) &
  (model_output[4,] < model_output[3,])
