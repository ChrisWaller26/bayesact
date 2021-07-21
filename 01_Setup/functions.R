
freq_form = fit_freq
sev_form = fit_sev

formula_to_stan =
  function(freq_form, sev_form){
    
    # Frequency parameters
    
    freq_resp = freq_form$resp
    
    freq_pforms = freq_form$pforms
    freq_pnames = names(freq_pforms)
    
    freq_family = as.character(freq_form$family$family)
    
    freq_dpars = freq_form$family$dpars
    
    freq_link_name = freq_form$family$link
    
    freq_inv_link_name = 
      case_when(
        freq_link_name == "log"      ~ "exp",
        freq_link_name == "identity" ~ "",
        freq_link_name == "logit"    ~ "inv_logit"
      )
    
    freq_link_name =
      case_when(
        freq_link_name == "identity" ~ "",
        TRUE ~ freq_link_name
      )
    
    # Severity parameters
    
    sev_formula = as.character(sev_form$formula[3])
    
    sev_resp = sev_form$resp
    
    sev_pforms = sev_form$pforms
    sev_dpars = sev_form$family$dpars
    sev_pnames = setdiff(names(sev_pforms), sev_dpars)
    
    sev_family = as.character(sev_form$family$family)
    
    sev_link_name = 
      c(sev_form$family$link,
        sev_form$family[[paste0("link_", sev_dpars[2])]]
      )
    
    sev_inv_link_name = 
      case_when(
        sev_link_name == "log"      ~ "exp",
        sev_link_name == "identity" ~ "",
        sev_link_name == "logit"    ~ "inv_logit"
      )
    
    sev_link_name =
      case_when(
        sev_link_name == "identity" ~ "",
        TRUE ~ sev_link_name
      )
    
    #' TODO: Change this so it converts the severity formula, regardless
    #' of what formula is used.
    
    sev_factors = gsub("[^0-9|a-z|A-Z]", "", sev_formula)
    sev_formula_stan = sev_formula
    
    for(n in seq(length(sev_pnames))){
      
      sev_factors =
        gsub(sev_pnames[n], " ", sev_factors)
      
    }
    
    sev_factors = unlist(str_split(str_squish(sev_factors), " "))
    
    for(n in seq(length(sev_pnames))){
      
      sev_formula_stan =
        gsub(
          sev_factors[n],
          paste0("C_", "_", freq_resp, "_" , n, "[n]"),
          gsub(sev_pnames[n], 
               paste0("nlp_", freq_resp, "_", freq_pnames[n], "[n]"), 
               sev_formula_stan)
        )
      
      
    }

    stan_formula =
      str_glue(
        "{freq_link_name}(
  {freq_inv_link_name}({freq_dpars[1]}_{freq_resp}) * 
  (1 - lognormal_cdf(
    ded, 
    
    {sev_formula_stan}, 
    
    {sev_inv_link_name[2]}(Intercept_{sev_dpars[2]}_{sev_resp} + 
          X_{freq_resp}_{freq_dpars[1]}[, 2:K_{sev_dpars[2]}_{sev_resp}] * 
          b_{sev_dpars[2]}_{sev_resp}
    )
  )
  )
)"
      )
    
    return(stan_formula)
    
  }


formula_to_stan(fit_freq, fit_sev)


