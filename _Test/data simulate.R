library(dplyr)
#' @title simulate data simple
#' @param n_pols number of policies
#' @param Intercept_lambda poisson rate
#' @param Intercept_mu  lognormal mu
#' @param b_lambda betas for lamda
#' @param b_mu   betas for mu
#' @param sigma  lognormal sigma
#' @param deductibles   deductibles amount
#' @export
data_sim_simple <- function(n_pols, Intercept_lambda, Intercept_mu, b_lambda, b_mu, sigma, deductibles) {
  
  pol_df <- tibble::tibble(
    pol_id = seq_len(n_pols), 
    freq_lvl = sample(seq_len(length(b_lambda)), size = n_pols, replace = T), 
    mu_lvl = sample(seq_len(length(b_mu)), size = n_pols, replace = T), 
    ded = sample(deductibles, size = n_pols, replace = T)
  ) %>% 
    dplyr::mutate(
      freq_factor = b_lambda[freq_lvl] + Intercept_lambda, 
      mu_factor = b_mu[mu_lvl] + Intercept_mu, 
      claim_fgu = rpois(dplyr::n(), lambda = freq_factor)
    )
  
  sev_df <- pol_df %>% 
    dplyr::filter(claim_fgu > 0) %>% 
    dplyr::group_nest(pol_id) %>% 
    dplyr::mutate(
      claim_id = purrr::map(data, ~seq_len(.$claim_fgu)),
      loss = purrr::map(data, ~rlnorm(.$claim_fgu, .$mu_factor, sigma))
    ) %>% 
    tidyr::unnest(cols = c(data, claim_id, loss)) %>% 
    dplyr::filter(loss > ded) %>% 
    dplyr::mutate(
      freq2 = 0, 
      claim_id = paste0(pol_id, "_", claim_id), 
      claimcount = 0
    )
  
   frq_df <- pol_df %>% 
     dplyr::select(pol_id, freq_lvl, mu_lvl, ded) %>% 
     dplyr::left_join(
       sev_df %>% dplyr::group_by(pol_id) %>% dplyr::summarise(claimcount = dplyr::n())
     ) %>% 
     dplyr::mutate(
       claimcount = dplyr::coalesce(claimcount, 0), 
       freq2 = 1, 
       loss = ded + 0.1, 
       claim_id = '0'
     )
   
   return(
     dplyr::bind_rows(
       frq_df, sev_df %>% dplyr::select(-freq_factor, -mu_factor, -claim_fgu)
     ) %>% 
       dplyr::mutate_at(
         dplyr::vars(dplyr::ends_with('lvl')), as.character
       )
   )
}
