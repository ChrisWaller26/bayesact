freq_formula = 
  bf(claimcount ~ f1,
     f1 ~ 1 + region,
     nl = TRUE)
sev_formula = 
  bf(loss | trunc(lb = ded) + cens(lim_exceed) ~ 
       s1,
     s1 ~ 1 + region,
     sigma ~ 1 + region,
     nl = TRUE
  )

freq_family = poisson
sev_family = lognormal

freq_data = freq_data_net
sev_data = sev_data

freq_link = c(mu = "log")
sev_link = c(mu = "identity", sigma = "log")

priors = c(prior(normal(0, 1),
                 class = b,
                 coef = Intercept,
                 resp = claimcount,
                 nlpar = f1),
           
           prior(normal(8, 1),
                 class = b,
                 coef = Intercept,
                 resp = loss,
                 nlpar = s1),

           prior(lognormal(0, 1),
                 class = Intercept,
                 dpar = sigma,
                 resp = loss)
           )

ded_name = "ded"


brms_freq_sev =
  function(
    freq_formula,
    sev_formula,
    freq_family = poisson,
    sev_family = lognormal,
    freq_data,
    sev_data,
    freq_link = c("log"),
    sev_link = c("identity", "log"),
    priors,
    ded_name = "ded",
    ...
    ){
    
    library(brms)
    library(dplyr)
    library(purrr)
    library(tidyr)
    library(stringr)
    library(rstan)
    library(rstanarm)
    library(readr)
    
    #### Multivariate Model ####
    
    fit_freq = 
      freq_formula + 
      freq_family(link = freq_link[[1]])
    
    fit_sev = 
      sev_formula + 
      sev_family(
        link = sev_link[[1]],
        link_sigma = sev_link[[2]]
      )
    
    mv_model_formula = fit_sev + fit_freq + set_rescor(FALSE)
    
    #### Join Data ####
    
    full_data = 
      bind_rows(
        freq_data %>%
          mutate(
            freq = 1,
            !!fit_sev$resp := get(ded_name) + 1 
          ),
        sev_data %>%
          mutate(
            freq = 0,
            !!fit_freq$resp := 0 
          )
      )
    
    # Replace NAs
    
    full_data =
      full_data %>%
      mutate_all(
        function(x){
          coalesce(x, get(paste0("as.", class(x)))(1))
        }
        )
      
    
    stanvars = c(
      stanvar(
        x = full_data[[ded_name]],
        name = "ded"
      ),
      stanvar(
        x = full_data$freq,
        name = "freq"
      )
    )
    
    ## Setup Stan Code and Data
    
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
    
    ## Extract Info from Formula Functions
    
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
    freq_feature_stan = grep(str_glue("C_{freq_resp}_"), 
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
        read_file("02_Model/template_model.stan"),
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
          str_locate(mv_model_code, 
                     str_glue("for \\(n in 1\\:N_{sev_resp}\\)"))[, 1] - 1
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
        ...
      )
    
    ## Convert back to BRMS fit object
    
    mv_model_fit$fit <- mv_model_fit_stan
    
    mv_model_fit <- rename_pars(mv_model_fit)
    
    return(mv_model_fit)
    
  }
