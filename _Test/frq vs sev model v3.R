library(dplyr)
library(brms)


source('model/data simulate.R')
set.seed(12345)
# full_data <- data_sim_simple(n_pols = 5000, Intercept_lambda = 0.5, Intercept_mu = 8, 
#                              b_lambda = c(0, 0.2), b_mu = c(0, -1, 1), sigma = 2, 
#                              deductibles = c(0, 5e3, 10e3, 25e3, 5e4))

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


full_data1 <- full_data %>% 
  dplyr::mutate(
    freq2 = ifelse(freq, 1, 0), 
    loss = ifelse(freq, ded + 0.1, loss),
    lim_exceed = dplyr::coalesce(lim_exceed, 0)
  )

fit_freq =
  bf(claimcount ~ f1 + freq2 * 0 ,
     f1 ~ 1 + region,
     nl = TRUE) +
  poisson()

fit_sev =
  bf(loss | trunc(lb = ded) + cens(lim_exceed) ~ s1,
     s1 ~ 1 + region,
     sigma ~ 1 + region,
     nl = TRUE
     ) +
  lognormal()

# fit_freq = 
#   bf(claimcount | subset(freq) ~ f1,
#      f1 ~ 1 + region,
#      nl = TRUE) + 
#   poisson()
# 
# fit_sev = 
#   bf(loss | subset(sev) + trunc(lb = ded) + cens(lim_exceed) ~ s1,
#      s1 ~ 1 + region,
#      sigma ~ 1 + region,
#      nl = TRUE
#   ) + 
#   lognormal()



stan_data <- brms::make_standata(
  formula = fit_sev + fit_freq +  set_rescor(FALSE),
  data = full_data1, 
  prior = c(prior(normal(0, 2), class = 'b', coef = 'Intercept', resp = 'claimcount', nlpar = 'f1'))
)

stan_code <- brms::make_stancode(
  formula = fit_sev + fit_freq + set_rescor(FALSE),
  data = full_data1, 
  prior = c(prior(normal(0, 2), class = 'b', coef = 'Intercept', resp = 'claimcount', nlpar = 'f1'), 
            prior(normal(5, 2), class = 'b', coef = 'Intercept', resp = 'loss', nlpar = 's1')
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
writeLines(replace_freq_block(stan_code), con = 'Model/modified stan code.stan')

library(rstan)
options(mc.cores = parallel::detectCores())
stanfit <- rstan::stan(#file = 'Model/modified stan code.stan', 
                   model_code = modified_stan_code,
                   data = stan_data,
                   chains = 1, iter = 200, warmup = 100, 
                   cores = 1
)
stanfit

# feed the Stan model back into brms
brmsfit <- brm( formula = fit_sev + fit_freq +  set_rescor(FALSE),
                data = full_data1, 
                prior = c(prior(normal(0, 2), class = 'b', coef = 'Intercept', resp = 'claimcount', nlpar = 'f1'), 
                          prior(normal(5, 2), class = 'b', coef = 'Intercept', resp = 'loss', nlpar = 's1')
                ), empty = TRUE)
brmsfit$fit <- stanfit
brmsfit <- rename_pars(brmsfit )

brmsfit