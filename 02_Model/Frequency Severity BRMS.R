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

# Defines a non-linear function for lambda to test model still works

lambda_fun = function(expo, region){
  c(EMEA = 1, USC = 2)[region]
} 

freq_data =
  data.frame(
    pol_id =  seq(freq_n),
    freq = TRUE,
    sev = FALSE,
    expo = runif(freq_n, 1, 100),
    ded = runif(freq_n, 0, 1000),
    lim = runif(freq_n, 25e3, 100e3),
    region = sample(regions, freq_n, replace = T)
  ) %>%
  mutate(
    freq_lambda = lambda_fun(expo, region),
    claimcount_fgu = rpois(freq_n, freq_lambda),
    loss = ded + 10000 # arbitrary loss which will be later ignored by model
  )

#### Simulate severity Data ####

mu_fun = function(expo, region){
  c(EMEA = 7, USC = 8)[region]
}

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
                   mu_fun(freq_data$expo[i], freq_data$region[i]), 
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
    loss_uncapped = loss,
    lim_exceed = 0,
    claimcount = coalesce(claimcount, 0L)
  )

#### Join Data ####

full_data = 
  bind_rows(freq_data_net, 
            sev_data)

#### Multivariate Model ####

freq_link = c(mu = "log")
sev_link = c(mu = "identity", sigma = "log")

fit_freq = 
  bf(claimcount ~ f1,
     f1 ~ 1 + region,
     nl = TRUE) + 
  poisson(link = freq_link[["mu"]])

fit_sev = 
  bf(loss | trunc(lb = ded) + cens(lim_exceed) ~ 
       s1 * expo / expo,
     s1 ~ 1 + region,
     sigma ~ 1 + region,
     nl = TRUE
  ) + 
  lognormal(
    link = sev_link[["mu"]],
    link_sigma = sev_link[["sigma"]]
  )

mv_model_formula = fit_sev + fit_freq + set_rescor(FALSE)

stanvars = c(
  stanvar(
    x = full_data$ded,
    name = "ded"
  ),
  stanvar(
    x = full_data$sev,
    name = "sev"
  ),
  stanvar(
    x = full_data$freq,
    name = "freq"
  )
  # ,
  # stanvar(
  #   scode = 
  #     "b_sigma_loss_Intercept = b_sigma_loss_Intercept + dot_product(means_X_sigma_loss, b_sigma_loss);",
  #   block = "genquant",
  #   position = "end"
  # )
)

priors = c(prior(normal(0, 1),
                 class = b,
                 coef = Intercept,
                 resp = claimcount,
                 nlpar = f1),
           
           # prior(normal(0, 1),
           #       class = b,
           #       coef = Intercept,
           #       resp = claimcount,
           #       nlpar = f2),
           
           prior(normal(8, 1),
                 class = b,
                 coef = Intercept,
                 resp = loss,
                 nlpar = s1),
           
           # prior(normal(8, 1),
           #       class = b,
           #       coef = Intercept,
           #       resp = loss,
           #       nlpar = s2),
           
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
  )

mv_model_data =
  make_standata(
    mv_model_formula,
    data = full_data,
    prior = priors,
    stanvars = stanvars
  )

mv_model_fit <- 
  brm( formula = mv_model_formula,
       data = full_data, 
       prior = priors, 
       empty = TRUE
  )

sev_dist = fit_sev$family$family
freq_dist = fit_freq$family$family

sev_resp = fit_sev$resp
freq_resp = fit_freq$resp

sev_arg =
  case_when(
    sev_dist == "lognormal" ~ c("mu", "sigma"),
    sev_dist == "gamma"     ~ c("alpha", "beta"),
    sev_dist == "normal"    ~ c("mu", "sigma"),
    sev_dist == "pareto"    ~ c("y_min", "alpha")
  )
freq_arg = c("mu")

sev_par = setdiff(names(fit_sev$pforms), sev_arg)
freq_par = setdiff(names(fit_freq$pforms), freq_arg)

inv_logit = function(x) exp(x) / (1 + exp(x))

sev_inv_link = 
  case_when(
    sev_link == "identity" ~ "",
    sev_link == "log"      ~ "exp",
    sev_link == "logit"    ~ "inv_logit"
  )
freq_inv_link = 
  case_when(
    freq_link == "identity" ~ "",
    freq_link == "log"      ~ "exp",
    freq_link == "logit"    ~ "inv_logit"
  )

sev_formula = as.character(fit_sev$formula[3])
freq_formula = as.character(fit_freq$formula[3])

sev_feature_stan = grep(str_glue("C_{sev_resp}_"), 
                      names(mv_model_data), 
                      value = TRUE,
                      fixed = TRUE)
freq_feature_stan = grep("C_{freq_resp}_", 
                      names(mv_model_data), 
                      value = TRUE, 
                      fixed = TRUE)

## Severity Features

sev_feature =
  gsub(
    '[^0-9|A-Z|a-z|_| ]',
    " ",
    str_squish(as.character(mv_model_fit$formula$forms[[sev_resp]]$formula[3]))
  ) %>%
  str_squish()

for(par in sev_par){
  
  sev_feature = trimws(gsub(par, "", sev_feature))
  
}

sev_feature = unique(str_split(sev_feature, " ")[[1]])

## Frequency Features

freq_feature =
  gsub(
    '[^0-9|A-Z|a-z|_| ]',
    " ",
    str_squish(as.character(mv_model_fit$formula$forms[[freq_resp]]$formula[3]))
  ) %>%
  str_squish()

for(par in freq_par){
  
  freq_feature = trimws(gsub(par, "", freq_feature))
  
}

freq_feature = unique(str_split(freq_feature, " ")[[1]])

## Severity Stan Formula

sev_formula_stan = 
  paste0(
    sev_inv_link[1],
    "(",
    sev_formula,
    ")"
    )

if(length(sev_feature_stan) >= 1){
  
  for(i in seq(length(sev_feature_stan))){
    
    sev_formula_stan =
      gsub(
        sev_feature[i],
        paste0(sev_feature_stan[i], "[n]"),
        sev_formula_stan
      )
    
  }
}

## Frequency Stan Formula

freq_formula_stan = 
  paste0(
    freq_inv_link[1],
    "(",
    freq_formula,
    ")"
  )

if(length(freq_feature_stan) >= 1){
  
  for(i in seq(length(freq_feature_stan))){
    
    freq_formula_stan =
      gsub(
        freq_feature[i],
        paste0(freq_feature_stan[i], "[n]"),
        freq_formula_stan
      )
    
  }
  
}



## Convert to Stan function format

for(par in sev_par){
  
  sev_formula_stan =
    gsub(
      par,
      paste0("nlp_", sev_resp, "_", par, "[n]"),
      sev_formula_stan
    )
  
}

for(par in freq_par){
  
  freq_formula_stan =
    gsub(
      par,
      paste0("nlp_", freq_resp, "_", par, "[n]"),
      freq_formula_stan
    )
  
}

## Modify Template Code

code_model_template =
  str_glue(
    paste(readLines("02_Model/template_model.stan"), collapse = "\n"),
    .open = "!!{",
    .close = "}!!"
    )

adjusted_code = 
  paste(
    substr(
      mv_model_code,
      1,
      str_locate(mv_model_code, "model \\{")[, 1] - 1
    ),

    substr(
      mv_model_code,
      str_locate(mv_model_code, "model \\{")[, 1],
      str_locate(mv_model_code, "for \\(n in 1\\:N_loss\\)")[, 1] - 1
    ),
    
    code_model_template,
    
    substr(
      mv_model_code,
      str_locate(mv_model_code, "// priors")[, 1],
      nchar(mv_model_code)
    ),
    
    sep = "\n"
    )


mv_model_fit_stan =
  stan(
    model_code = adjusted_code,
    data = mv_model_data,
    chains = 1,
    iter = 1000,
    warmup = 250,
    control = 
      list(adapt_delta = 0.999,
           max_treedepth = 15)
  )

## Convert back to BRMS fit object

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
  file = "02_Model/mv_model_fit.RData"
)

base::load(file = "02_Model/mv_model_fit.RData")

sev_output =
  full_data %>%
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
  ) %>%
  filter(sev)

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
          filter(freq) %>%
          left_join(
            sev_output_mean,
            by = "region"
          ) %>%
          mutate(ded = 0),
        resp = "claimcount"
      ) %>%
      colMeans() %>%
      exp(),
    
    lambda_pred_q025 = 
      posterior_linpred(
        mv_model_fit,
        newdata = 
          full_data %>%
          filter(freq) %>%
          left_join(
            sev_output_mean,
            by = "region"
          ) %>%
          mutate(ded = 0),
        resp = "claimcount"
      ) %>%
      exp() %>%
      apply(2, function(x) quantile(x, 0.025)),
    
    lambda_pred_q975 = 
      posterior_linpred(
        mv_model_fit,
        newdata = 
          full_data %>%
          filter(freq) %>%
          left_join(
            sev_output_mean,
            by = "region"
          ) %>%
          mutate(ded = 0),
        resp = "claimcount"
      ) %>%
      exp() %>%
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
    mu_actual     = mu_fun(region),
    sigma_actual  = sev_sigma_vec[region],
    lambda_actual = lambda_fun(region)
  ) %>%
  select(sort(names(.)))




