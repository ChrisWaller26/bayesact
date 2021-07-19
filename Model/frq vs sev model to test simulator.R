library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(tidysimloss)
###install tidysimloss by: devtools::install_github('Atan1988/tidysimloss')

policy_alist <- alist(
  Exposures ~ rtrunc(FUN = 'norm', Att = 5e5, mean, sd),
  mean = a + b_Industries[Industries]
)

policy_parameters_alist <- alist(
  ded ~ rdiscrete(ded_list),
  limit ~ rdiscrete(limit_list),
  Industries ~ rdiscrete(Industry_options),
  Industry_options = c('Energy', 'Construction', 'Healthcare'),
  Eff_yrs = seq(2010, 2016, 1),
  sd = 10e6, b_Industries = c(10e6, 25e6, 0), a = 50e6,
  ded_list = c(0, 2500, 5000, 25000, 50000, 100000),
  limit_list = c(1e6, 2e6, 3e6)
)


frequency_alist <- alist(
  total_claims ~ rpois(lambda),
  lambda = a_lambda * Exposures^b_lambda_industry,
  #a_lambda = 0.05,
  b_lambda_industry = case_when(Industries == 'Healthcare' ~ b3
                                , Industries == 'Construction' ~ b2
                                , TRUE ~ b1) #c(0.1, 0.05, 0.025)
)

frequency_params_components_alist <- alist(
  a_lambda ~ rnorm(mean = 0.3, sd = 0.00),
  b1 ~ rnorm(mean = 0.1, sd = 0.00),
  b2 ~ rnorm(mean = 0.075, sd = 0.00),
  b3 ~ rnorm(mean = 0.05, sd = 0.00)
)

severity_alist <- alist(
  loss ~ rlnorm(meanlog = mu, sdlog = 2),
  mu = a_mu + b_mu,
  b_mu = case_when(Industries == 'Healthcare' ~ b3
                   , Industries == 'Construction' ~ b2
                   , TRUE ~ b1)
)

severity_params_alist <- alist(
  a_mu ~ rnorm(mean = 9, sd = 0.00),
  b1 ~ rnorm(mean = 0, sd = 0.00),
  b2 ~ rnorm(mean = 1, sd = 0.00),
  b3 ~ rnorm(mean = -1, sd = 0.00)
)


N_Policies <- 1e4
full_data <- base_simulator(N_Policies,
                            policy_exprs = policy_alist, policy_parameters = policy_parameters_alist,
                            frequency_exprs = frequency_alist, frequnecy_parameters = frequency_params_components_alist,
                            severity_exprs = severity_alist, severity_paramters = severity_params_alist)

library(brms)
library(dplyr)
fit_freq =
  bf(claimcount ~ f1 + freq2 * 0 ,
     f1 ~ 1 + Industries:log(Exposures),
     nl = TRUE) +
  poisson()

fit_sev =
  bf(loss | trunc(lb = ded) + cens(lim_exceed) ~ s1,
     s1 ~ 1 + Industries,
     sigma ~ 1,
     nl = TRUE
  ) +
  lognormal()

mvfit_formula <- fit_sev + fit_freq +  set_rescor(FALSE)

stan_data <- brms::make_standata(
  formula = mvfit_formula,
  data = full_data
)

stan_code <- brms::make_stancode(
  formula = mvfit_formula,
  data = full_data,
  prior = c(prior(normal(-1, 0.1), class = 'b', resp = 'claimcount', coef = 'Intercept', nlpar = 'f1'),
            prior(normal(5, 2), class = 'b', resp = 'loss', nlpar = 's1')
  )
)


replace_freq_block <- function(stan_code) {
  raw_code <- stan_code %>% as.character() %>% strsplit('\\\n') %>% .[[1]]
  ###modify the mu_claimcount line
  raw_code <- stringr::str_replace(raw_code, 'mu_claimcount\\[n\\] = nlp_claimcount_f1\\[n\\]',
                                   'mu_claimcount[n] = nlp_claimcount_f1[n] + log(1 - lognormal_cdf(lb_loss[n], mu_loss[n], sigma_loss[n]))')
  
  ###identify likelihood line for loss
  loss_likelihood_start <- which(grepl('special treatment of censored', raw_code))
  raw_code1 <- c(raw_code[1:(loss_likelihood_start - 1)],
                 c(paste0(stringr::str_pad(" ", 5, 'left', ' '), 'if (C_claimcount_1[n] == 1) {'),
                   paste0(stringr::str_pad(" ", 8, 'left', ' '), 'target += poisson_log_lpmf(Y_claimcount[n] | mu_claimcount[n]);'),
                   paste0(stringr::str_pad(" ", 5, 'left', ' '), '} else {' )),
                 raw_code[loss_likelihood_start:length(raw_code)])
  ###replace poisson likelihood with closing bracket
  poisson_likelihood_position <- which(grepl('poisson_log_lpmf', raw_code1))
  raw_code1[max(poisson_likelihood_position)] <- '   }'
  return(raw_code1 %>% paste(collapse = " \n"))
}

modified_stan_code <- replace_freq_block(stan_code)
library(rstan)
options(mc.cores = parallel::detectCores())
stanfit <- rstan::stan(#file = 'Model/modified stan code.stan',
  model_code = modified_stan_code,
  data = stan_data,
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4
)
stanfit

# feed the Stan model back into brms
brmsfit <- brm( formula = mvfit_formula,
                data = full_data,
                prior = c(
                  prior(normal(-1, 0.5), class = 'b', resp = 'claimcount', coef = 'Intercept', nlpar = 'f1'),
                  prior(normal(5, 2), class = 'b', resp = 'loss', nlpar = 's1')
                ), empty = TRUE)
brmsfit$fit <- stanfit
brmsfit <- rename_pars(brmsfit )

brmsfit
