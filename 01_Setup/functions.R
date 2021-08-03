
brms_freq_sev =
  function(
    freq_formula,
    sev_formula,
    freq_family,
    sev_family,
    freq_data,
    sev_data,
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
    
    #### Error Checks ####
    
    if(is.null(freq_formula)){
      
      stop("Frequency Formula Required")
      
    }
    
    if(is.null(sev_formula)){
      
      stop("Severity Formula Required")
      
    }
    
    if(is.null(freq_family)){
      
      stop("Frequency Family Required")
      
    }
    
    if(is.null(sev_family)){
      
      stop("Severity Family Required")
      
    }
    
    if(is.null(freq_data)){
      
      stop("Frequency Data Required")
      
    }
    
    if(is.null(sev_data)){
      
      stop("Severity Data Required")
      
    }
    
    if(is.null(priors)){
      
      stop("Priors Required")
      
    }
    
    if(is.null(ded_name)){
      
      stop("Name of Deductible Column Required")
      
    }
    
    #### Multivariate Model ####
    
    sev_dist = 
      case_when(
        sev_family$family == "gaussian" ~ "normal",
        sev_family$family == "Gamma"    ~ "gamma",
        TRUE ~ sev_family$family
      )
    
    freq_dist = 
      case_when(
        freq_family$family == "negbinomial" ~ "neg_binomial",
        TRUE ~ freq_family$family
      )
    
    sev_resp = sev_formula$resp
    freq_resp = freq_formula$resp
    
    ## Convert Linear Model to Non-Linear 
    
    if(!attr(freq_formula$formula, "nl")){
      
      freq_formula =
        bf(
          as.formula(
            paste0(as.character(freq_formula$formula[2]), " ~ f1")
          ),
          as.formula(
            paste0("f1 ~ ", as.character(freq_formula$formula[3]))
          ),
          freq_formula$pforms,
          nl = TRUE
        )
      
      priors = 
        priors %>%
        mutate(
          nlpar = 
            case_when(
              (resp == freq_resp) & (dpar == "") ~ "f1",
              TRUE ~ nlpar
            ),
          coef =
            case_when(
              (resp == freq_resp) & 
                (dpar == "") & 
                class == "Intercept" ~ "Intercept",
              TRUE ~ coef
            ),
          class = 
            case_when(
              (resp == freq_resp) & (dpar == "") ~ "b",
              TRUE ~ class
            )
        )
      
    }
    
    if(!attr(sev_formula$formula, "nl")){
      
      sev_formula =
        bf(
          as.formula(
            paste0(as.character(sev_formula$formula[2]), " ~ s1")
          ),
          as.formula(
            paste0("s1 ~ ", as.character(sev_formula$formula[3]))
          ),
          sev_formula$pforms,
          nl = TRUE
        )
      
      priors = 
        priors %>%
        mutate(
          nlpar = 
            case_when(
              (resp == sev_resp) & (dpar == "") ~ "s1",
              TRUE ~ nlpar
            ),
          coef =
            case_when(
              (resp == sev_resp) & 
                (dpar == "") & 
                class == "Intercept" ~ "Intercept",
              TRUE ~ coef
            ),
          class = 
            case_when(
              (resp == sev_resp) & (dpar == "") ~ "b",
              TRUE ~ class
            )
        )
      
    }
    
    
    freq_link = 
      freq_family[
        grep("link", 
             setdiff(names(freq_family), c("linkinv", "linkfun")),
             value = TRUE
        )
      ] %>%
      unlist()
    
    sev_link = 
      sev_family[
        grep("link", 
             setdiff(names(sev_family), c("linkinv", "linkfun")),
             value = TRUE
        )
      ] %>%
      unlist()
    
    #### Convert to censored if uncensored
    
    ## Frequency
    
    if(!grepl("|", freq_formula$formula[2], fixed = TRUE)){
      
      freq_formula =
        bf(
          as.formula(
            str_c(
              as.character(freq_formula$formula[2]),
              " | cens(no_censoring) + trunc(lb = 0, ub = 999999999L)~ ", 
              as.character(freq_formula$formula[3])
            )
          ),
          freq_formula$pforms,
          nl = TRUE
        )
      
    }
    
    if(!grepl("cens", freq_formula$formula[2])){
      
      freq_formula =
        bf(
          as.formula(
            str_c(
              gsub(
                "|",
                " | cens(no_censoring) + ",
                as.character(freq_formula$formula[2]),
                fixed = TRUE
              ),
              " ~ ", 
              as.character(freq_formula$formula[3])
            )
          ),
          freq_formula$pforms,
          nl = TRUE
        )
      
    }
    
    ## Severity
    
    if(!grepl("|", sev_formula$formula[2], fixed = TRUE)){
      
      sev_formula =
        bf(
          as.formula(
            str_c(
              as.character(sev_formula$formula[2]),
              " | cens(no_censoring) + trunc(lb = 0, ub = 1e99) ~ ", 
              as.character(sev_formula$formula[3])
            )
          ),
          sev_formula$pforms,
          nl = TRUE
        )
      
    }
    
    if(!grepl("cens", sev_formula$formula[2])){
      
      sev_formula =
        bf(
          as.formula(
            str_c(
              gsub(
                "|",
                " | cens(no_censoring) + ",
                as.character(sev_formula$formula[2]),
                fixed = TRUE
              ),
              " ~ ", 
              as.character(sev_formula$formula[3])
            )
          ),
          sev_formula$pforms,
          nl = TRUE
        )
      
    }
    
    #### Convert unbounded model to truncated
    
    ## Frequency
    
    if(grepl("lb =|lb=|trunc", freq_formula$formula[2]) &
       !grepl("ub =|ub=", freq_formula$formula[2])){
      
      freq_formula =
        bf(
          as.formula(
            str_c(
              gsub(
                "trunc(",
                "trunc(ub = 999999999L, ", 
                as.character(freq_formula$formula[2]),
                fixed = TRUE 
              ),
              " ~ ", 
              as.character(freq_formula$formula[3])
            )
          ),
          freq_formula$pforms,
          nl = TRUE
        )
      
    }else if(!grepl("lb =|lb=", freq_formula$formula[2]) &
             grepl("ub =|ub=", freq_formula$formula[2])){
      
      freq_formula =
        bf(
          as.formula(
            str_c(
              gsub(
                "trunc(",
                "trunc(lb = 0L, ", 
                as.character(freq_formula$formula[2]),
                fixed = TRUE
              ),
              " ~ ", 
              as.character(freq_formula$formula[3])
            )
            
          ),
          freq_formula$pforms,
          nl = TRUE
        )
      
    }else if(!grepl("trunc", freq_formula$formula[2])){
      
      freq_formula =
        bf(
          as.formula(
            str_c(
              as.character(freq_formula$formula[2]),
              " + trunc(lb = 0L, ub = 999999999L) ~ ",
              as.character(freq_formula$formula[3])
            )
          ),
          freq_formula$pforms,
          nl = TRUE
        )
      
    }
    
    ## Severity
    
    if(grepl("lb =|lb=|trunc", sev_formula$formula[2]) &
       !grepl("ub =|ub=", sev_formula$formula[2])){
      
      sev_formula =
        bf(
          as.formula(
            str_c(
              gsub(
                "trunc(",
                "trunc(ub = 1e99, ", 
                as.character(sev_formula$formula[2]),
                fixed = TRUE
              ),
              " ~ ",
              as.character(sev_formula$formula[3])
            )
            
          ),
          sev_formula$pforms,
          nl = TRUE
        )
      
    }else if(!grepl("lb =|lb=", sev_formula$formula[2]) &
             grepl("ub =|ub=", sev_formula$formula[2])){
      
      sev_formula =
        bf(
          as.formula(
            str_c(
              gsub(
                "trunc(",
                "trunc(lb = 0, ", 
                as.character(sev_formula$formula[2]),
                fixed = TRUE 
              ),
              " ~ ",
              as.character(sev_formula$formula[3])
            )
          ),
          sev_formula$pforms,
          nl = TRUE
        )
      
    }else if(!grepl("trunc", sev_formula$formula[2])){
      
      sev_formula =
        bf(
          as.formula(
            str_c(
              as.character(sev_formula$formula[2]),
              " + trunc(lb = 0, ub = 1e99) ~ ",
              as.character(sev_formula$formula[3])
            )
          ),
          sev_formula$pforms,
          nl = TRUE
        )
      
    }
    
    ## Create multivariate model
    
    fit_freq = 
      freq_formula + 
      freq_family
    
    fit_sev = 
      sev_formula + 
      sev_family
    
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
      ) %>%
      mutate(
        no_censoring = 0
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
      ),
      stanvar(
        scode = 
          str_glue(
            "    vector[N_{sev_resp}] sev_target; \n
               vector[N_{freq_resp}] freq_target;"
          ),
        block = "model"
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
    
    sev_arg = fit_sev$family$dpars
    freq_arg = fit_freq$family$dpars
    
    if(sev_dist == "gamma"){
      
      sev_arg = rev(sev_arg)
      
    }
    
    sev_arg_stan = 
      str_flatten(
        str_c(
          sev_arg, 
          "_", 
          sev_resp, 
          "[n]"
        ), 
        ", ")
    
    freq_arg_stan = 
      str_flatten(
        str_c(
          freq_arg, 
          "_", 
          freq_resp, 
          "[n]"
        ), 
        ", ")
    
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
    
    sev_formula_r = as.character(fit_sev$formula[3])
    freq_formula_r = as.character(fit_freq$formula[3])
    
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
        sev_formula_r,
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
        freq_formula_r,
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
                     str_c(
                       "for \\(n in 1\\:N_", sev_resp, "\\) \\{\n    // special treatment"))[, 1] - 1
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
