library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(tidysimloss)
library(bayesact)
library(brms)
###install tidysimloss by: devtools::install_github('Atan1988/tidysimloss')

# policy_alist <- alist(
#   Exposures ~ rnorm(mean = 0, sd = 1)
# )
# 
# policy_parameters_alist <- alist(
#   ded ~ rdiscrete(ded_list),
#   limit ~ rdiscrete(limit_list),
#   region ~ rdiscrete(region_options),
#   region_options = c("EMEA", "USC"),
#   ded_list = c(0), #c(0, 2500, 5000, 25000, 50000, 100000),
#   limit_list = c(1e6, 2e6, 3e6)
# )
# 
# 
# frequency_alist <- alist(
#   total_claims ~ rpois(lambda),
#   lambda = a_lambda + b_lambda_region,
#   #a_lambda = 0.05,
#   b_lambda_region = case_when(region == 'USC' ~ b2
#                                 , TRUE ~ b1) #c(0.1, 0.05, 0.025)
# )
# 
# frequency_params_components_alist <- alist(
#   a_lambda ~ rnorm(mean = 2, sd = 0.00),
#   b1 ~ rnorm(mean = 0, sd = 0.00),
#   b2 ~ rnorm(mean = 1, sd = 0.00)
# )
# 
# severity_alist <- alist(
#   loss ~ rgamma(shape = 0.5, scale = mu),
#   mu = exp(a_mu + b_mu),
#   b_mu = case_when(region == 'USC' ~ b2
#                    , TRUE ~ b1)
# )
# 
# severity_params_alist <- alist(
#   a_mu ~ rnorm(mean = 9, sd = 0.00),
#   b1 ~ rnorm(mean = 0, sd = 0.00),
#   b2 ~ rnorm(mean = 1, sd = 0.00)
# )
# 
# 
# N_Policies <- 1e4
# full_data <- base_simulator(N_Policies,
#                             policy_exprs = policy_alist, policy_parameters = policy_parameters_alist,
#                             frequency_exprs = frequency_alist, frequnecy_parameters = frequency_params_components_alist,
#                             severity_exprs = severity_alist, severity_paramters = severity_params_alist)
# 
# freq_data_net <- full_data %>% dplyr::filter(freq2 == 1) %>% dplyr::rename(freq = freq2)
# sev_data <- full_data %>% dplyr::filter(freq2 == 0) %>% dplyr::rename(freq = freq2) 

sev_data <- tibble(
  claim_id = seq_len(1e4)
) %>% 
  dplyr::mutate(
    region = sample(c("EMA", "USC"), size = dplyr::n(), replace = T), 
    scale = dplyr::case_when(
      region == "EMA" ~ 9, 
      T ~ 10
    ), 
    ded = 0, 
    lim_exceed = 0,
    loss = rgamma(dplyr::n(), shape = 0.5, scale = exp(scale))
  )
  
sev_formula <-   bf(loss | trunc(lb = ded) + cens(lim_exceed) ~ s1/1000,
                    s1 ~ 1 + region, 
                    nl = T)

sev_formula2 <- bf(loss | trunc(lb = ded) + cens(lim_exceed) ~ s1 + 5,
                  s1 ~ 1 + region, 
                  nl = T)

brmsfit <- brm(
  sev_formula, 
  data = sev_data, # %>% dplyr::mutate(loss = loss / 1000), 
  family = Gamma(),
  prior = c(
             # prior(constant(9), 
             #       nlpar = a1), 
             # 
             # prior(normal(0, 1),
             #       class = b,
             #       nlpar = s1)
        prior(normal(0, 1), class = b, nlpar = s1)
  ),
  chains = 1,
  iter = 2000,
  warmup = 1000,
  refresh = 50,
  control =
    list(adapt_delta = 0.8,
         max_treedepth = 10)
)
brmsfit

brmsfit2 <- brm(
  sev_formula2, 
  data = sev_data, # %>% dplyr::mutate(loss = loss / 1000), 
  family = Gamma(link = 'log'),
  prior = c(
    # prior(constant(9), 
    #       nlpar = a1), 
    # 
    # prior(normal(0, 1),
    #       class = b,
    #       nlpar = s1)
    prior(normal(0, 1), class = b, nlpar = s1)
  ),
  chains = 1,
  iter = 2000,
  warmup = 1000,
  refresh = 50,
  control =
    list(adapt_delta = 0.8,
         max_treedepth = 10)
)
brmsfit2

mv_model_fit =
  brms_freq_sev(

    freq_formula =
      bf(claimcount ~ 1 + region),

    sev_formula =
      bf(loss | trunc(lb = ded) + cens(lim_exceed) ~ a1 + s1,
          a1 ~ 1, 
          s1 ~ 1 + region, 
          nl = T
      ),

    freq_family = poisson(),
    sev_family = Gamma(link = 'log'),

    freq_data = freq_data_net,
    sev_data = sev_data,

    priors = c(prior(normal(0, 1),
                     class = b,
                     coef = Intercept, 
                     nlpar = f1, 
                     resp = claimcount),

               prior(normal(0, 1),
                     class = b,
                     nlpar = f1, 
                     resp = claimcount),

               prior(constant(8), 
                     nlpar = a1, 
                     resp = loss), 
               
               prior(normal(0, 1),
                     class = b,
                     coef = Intercept, 
                     nlpar = s1,
                     resp = loss)
    ),

    ded_name = "ded",

    chains = 1,
    iter = 1000,
    warmup = 250,
    refresh = 5,
    control =
      list(adapt_delta = 0.999,
           max_treedepth = 15)
  )

