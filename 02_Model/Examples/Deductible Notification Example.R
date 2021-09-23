#### Load Packages ####

library(brms)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(rstan)
library(rstanarm)
library(readr)

library(renv)

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
    ded = 2e3, # runif(freq_n, 1e3, 5e3),
    lim = runif(freq_n, 50e3, 100e3),
    region = sample(regions, freq_n, replace = T)
  ) %>%
  mutate(
    notified_size = ded * 0.5,
    freq_mu = freq_mu_fun(expo, region),
    claimcount_fgu = 
      rpois(freq_n, freq_mu)
  )

#### Simulate severity Data ####

mu_fun = function(expo, region){
  c(EMEA = 8, USC = 9)[region]
}

sev_par2_vec = exp(c(EMEA = 0, USC = 0.5))

sev_data =
  data.frame(
    ded = rep(freq_data$ded,
              freq_data$claimcount_fgu),
    notified_size = 
      rep(freq_data$notified_size,
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
                   mu_fun(freq_data$expo[i],
                            freq_data$region[i]),
                   sev_par2_vec[freq_data$region[i]]
                   )
            
          }
        )
      )
  ) %>%
  mutate(
    pol_id = rep(seq(freq_n), freq_data$claimcount_fgu)
  ) %>% 
  filter(
    loss_uncapped > notified_size
  ) %>%
  mutate(
    claim_id = row_number(),
    
    # Losses are capped at the limit - amounts above this are unknown
    
    loss = pmin(loss_uncapped, lim),
    
    #' If loss > ded, then ground-up value of loss is known
    #' If loss < notified_size, it will not be reported at all
    #' If notified_size <= loss <= ded, then loss will be notified
    #' but we won't know the exact size (it will just be reported as a
    #' zero loss)
    
    censor = 
      case_when(
        loss_uncapped >= lim ~ 1L,  # Right-censored
        loss_uncapped <= ded ~ 2L,  # Interval censored
        TRUE ~ 0L
      ),
    
    loss = 
      case_when(
        loss_uncapped < ded ~ notified_size,
        TRUE ~ loss
      ),
    
    weight = 1e3 / nrow(.)
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
    claimcount = coalesce(claimcount, 0),
    weight = 1e3 / nrow(.)
  )

#### Run Model ####

mv_model_fit =
  brms_freq_sev(
    
    freq_formula = 
      bf(claimcount | weights(weight) ~ 1 + region),
    
    sev_formula = 
      bf(loss | 
           trunc(lb = notified_size) + 
           cens(censor, y2 = ded) + 
           weights(weight) ~ 
           1 + region,
         sigma ~ 1 + region
      ),
    
    freq_family = poisson(),
    sev_family = lognormal(),
    
    freq_data = freq_data_net,
    sev_data = sev_data,
    
    priors = c(prior(normal(0, 1),
                     class = Intercept,
                     resp = claimcount),
               
               prior(normal(0, 1),
                     class = b,
                     resp = claimcount),
               
               prior(normal(8, 1),
                     class = Intercept,
                     resp = loss),
               
               prior(normal(0, 0.5),
                     class = Intercept,
                     dpar = sigma,
                     resp = loss),
               
               prior(normal(0, 0.5),
                     class = b,
                     dpar = sigma,
                     resp = loss)
    ),
    
    ded_name = "notified_size",
    ded_adj_min = 0.0001,
    use_cmdstan = T,
    
    chains = 1,
    iter = 300,
    warmup = 150,

    refresh = 100,
    adapt_delta = 0.99,
    max_treedepth = 15,
    
    mle = FALSE,
    sample_prior = "no",
    freq_adj_fun = NULL,
    stanvars     = NULL
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
    
    sigma_emea = exp(b_sigma_loss_Intercept), 
    sigma_usc  = exp(b_sigma_loss_Intercept 
                     + b_sigma_loss_regionUSC
                     ),
    
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
      
      sigma_emea = sev_par2_vec["EMEA"], 
      sigma_usc  = sev_par2_vec["USC"],
      
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
