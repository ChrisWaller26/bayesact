
library(cmdstanr)
library(dplyr)


N = 1000
par_lambda = 0.5
par_mu = 10
par_sigma = 1

rlnormpois =
  function(n, lambda, mu, sigma){
    sapply(rpois(n, lambda),
           function(x) sum(rlnorm(x, mu, sigma)))
  }

agg_data =
  data.frame(
    pol_id = seq(N),
    loss_agg = rlnormpois(N, par_lambda, par_mu, par_sigma)
  )

stan_data =
  list(
    N = N,
    loss_agg = loss_agg
  )



