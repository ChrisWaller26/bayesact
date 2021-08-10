#' Bayesian Frequency Severity Model with BRMS
#'
#' @description
#' This function allows you to create a multivariate Bayesian
#' frequency-severity model in which the frequency parameter is adjusted
#' by the survival function of the severity model at the deductible/excess/attachment point.
#'
#' This allows for both linear and non-linear functions of the variates for
#' both the frequency and severity model components.
#'
#' @param freq_formula BRMS Formula; Linear/Non-linear formula for frequency model
#' @param sev_formula  BRMS Formula; Linear/Non-linear formula for severity model
#' @param freq_family  Family; Family for frequency model
#' @param sev_family   Family; Family for severity model
#' @param freq_data    Data Frame; The data required for the frequency model
#' @param sev_data     Data Frame; The data required for the severity model
#' @param priors       BRMS Prior; The set of priors for both the frequency and severity models
#' @param ded_name     Character; The column name for the deductible/excess/attachment point in the frequency data
#' @param ...          Additional accepted BRMS fit parameters
#' @return             BRMS Fit
#'
#' @examples
#'
#' #### Simulate Frequency Data ####
#'
#' options(stringsAsFactors = FALSE,
#'         mc.cores = parallel::detectCores())
#'
#' #' Assuming one rating factor: region.
#'
#' set.seed(123456)
#'
#' # Region Names
#'
#' regions = c("EMEA", "USC")
#'
#' # Number of frequency samples
#'
#' freq_n = 5e3
#'
#' # Defines a function for lambda
#'
#' freq_lambda = exp(c(EMEA = 0.5, USC = 1))
#'
#' # Generate samples for ground-up frequency data
#'
#' freq_data =
#'   data.frame(
#'     pol_id =  seq(freq_n),
#'     ded = runif(freq_n, 1e3, 5e3),
#'     lim = runif(freq_n, 50e3, 100e3),
#'     region = sample(regions, freq_n, replace = T)
#'   ) %>%
#'   mutate(
#'     freq_lambda = freq_lambda[region],
#'     claimcount_fgu =
#'       rpois(freq_n, freq_lambda)
#'   )
#'
#' #### Simulate severity Data ####
#'
#' mu_vec = c(EMEA = 8, USC = 9)
#' sigma_vec = exp(c(EMEA = 0, USC = 0.4))
#'
#' sev_data =
#'   data.frame(
#'     ded = rep(freq_data$ded,
#'               freq_data$claimcount_fgu),
#'     lim = rep(freq_data$lim,
#'               freq_data$claimcount_fgu),
#'     region = rep(freq_data$region,
#'                  freq_data$claimcount_fgu)
#'   ) %>%
#'   mutate(
#'     loss_uncapped =
#'       unlist(
#'         lapply(
#'           seq(freq_n),
#'           function(i){
#'
#'             rlnorm(freq_data$claimcount_fgu[i],
#'                    mu_vec[freq_data$region[i]],
#'                    sigma_vec[freq_data$region[i]]
#'             )
#'
#'           }
#'         )
#'       )
#'   ) %>%
#'   mutate(
#'     pol_id = rep(seq(freq_n), freq_data$claimcount_fgu)
#'   ) %>%
#'   filter(
#'     loss_uncapped > ded
#'   ) %>%
#'   mutate(
#'     claim_id = row_number(),
#'     lim_exceed = as.integer(loss_uncapped >= lim),
#'     loss = pmin(loss_uncapped, lim)
#'   )
#'
#' # Frequency data filtered for losses below the deductible
#'
#' freq_data_net =
#'   freq_data %>%
#'   left_join(
#'     sev_data %>%
#'       group_by(
#'         pol_id
#'       ) %>%
#'       summarise(
#'         claimcount = n()
#'       ) %>%
#'       ungroup(),
#'     by = "pol_id"
#'   ) %>%
#'   mutate(
#'     claimcount = coalesce(claimcount, 0)
#'   )
#'
#' #### Run Model ####
#'
#' mv_model_fit =
#'   brms_freq_sev(
#'
#'     freq_formula =
#'       bf(claimcount ~ 1 + region),
#'
#'     sev_formula =
#'       bf(loss | trunc(lb = ded) + cens(lim_exceed) ~
#'            1 + region,
#'          sigma ~ 1 + region
#'       ),
#'
#'     freq_family = poisson(),
#'     sev_family = lognormal(),
#'
#'     freq_data = freq_data_net,
#'     sev_data = sev_data,
#'
#'     priors = c(prior(normal(0, 1),
#'                      class = Intercept,
#'                      resp = claimcount),
#'
#'                prior(normal(0, 1),
#'                      class = b,
#'                      resp = claimcount),
#'
#'                prior(normal(8, 1),
#'                      class = Intercept,
#'                      resp = loss),
#'
#'                prior(lognormal(0, 1),
#'                      class = Intercept,
#'                      dpar = sigma,
#'                      resp = loss),
#'
#'                prior(normal(0, 1),
#'                      class = b,
#'                      dpar = sigma,
#'                      resp = loss)
#'     ),
#'
#'     ded_name = "ded",
#'
#'     chains = 1,
#'     iter = 1000,
#'     warmup = 250,
#'     refresh = 50,
#'     control =
#'       list(adapt_delta = 0.999,
#'            max_treedepth = 15)
#'   )
#'
#' #### Results ####
#'
#' model_post_samples =
#'   posterior_samples(
#'     mv_model_fit
#'   ) %>%
#'   transmute(
#'     s1_emea = b_loss_s1_Intercept,
#'     s1_usc  = b_loss_s1_Intercept +
#'       b_loss_s1_regionUSC,
#'
#'     sigma_emea = exp(b_sigma_loss_Intercept),
#'     sigma_usc  = exp(b_sigma_loss_Intercept +
#'                        b_sigma_loss_regionUSC),
#'
#'     f1_emea = b_claimcount_f1_Intercept,
#'     f1_usc  = b_claimcount_f1_Intercept +
#'       b_claimcount_f1_regionUSC
#'   )
#'
#' model_output =
#'   model_post_samples %>%
#'   sapply(
#'     function(x) c(lower = quantile(x, 0.025),
#'                   mean  = mean(x),
#'                   upper = quantile(x, 0.975))
#'   ) %>%
#'   as.data.frame() %>%
#'   bind_rows(
#'     data.frame(
#'       s1_emea = 8,
#'       s1_usc  = 9,
#'
#'       sigma_emea = exp(0),
#'       sigma_usc  = exp(0.4),
#'
#'       f1_emea = 0.5,
#'       f1_usc  = 1
#'     )
#'   )
#'
#' @export
#'
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

    template_stan_code =

    '

  for(n in 1:N_!!{freq_resp}!!){

    mu_!!{freq_resp}!![n] =
      mu_!!{freq_resp}!![n] *
       (1 - !!{sev_dist}!!_cdf(ded[n], !!{sev_arg_stan}!!));

  }

  for (n in 1:N_!!{sev_resp}!!) {
    // special treatment of censored data
      if (cens_!!{sev_resp}!![n] == 0) {
        sev_target[n] = !!{sev_dist}!!_lpdf(Y_!!{sev_resp}!![n] | !!{sev_arg_stan}!!) -
        log_diff_exp(!!{sev_dist}!!_lcdf(ub_!!{sev_resp}!![n] |!!{sev_arg_stan}!!), !!{sev_dist}!!_lcdf(lb_!!{sev_resp}!![n] | !!{sev_arg_stan}!!));
      } else if (cens_!!{sev_resp}!![n] == 1) {
        sev_target[n] = !!{sev_dist}!!_lccdf(Y_!!{sev_resp}!![n] | !!{sev_arg_stan}!!) -
        log_diff_exp(!!{sev_dist}!!_lcdf(ub_!!{sev_resp}!![n] | !!{sev_arg_stan}!!), !!{sev_dist}!!_lcdf(lb_!!{sev_resp}!![n] | !!{sev_arg_stan}!!));
      } else if (cens_!!{sev_resp}!![n] == -1) {
        sev_target[n] = !!{sev_dist}!!_lcdf(Y_!!{sev_resp}!![n] | !!{sev_arg_stan}!!) -
        log_diff_exp(!!{sev_dist}!!_lcdf(ub_!!{sev_resp}!![n] | !!{sev_arg_stan}!!), !!{sev_dist}!!_lcdf(lb_!!{sev_resp}!![n] | !!{sev_arg_stan}!!));
      }

      // if(is_nan(sev_target[n]) || is_inf(sev_target[n])){
      //
      //   print(
      //     "Y_loss: ", Y_loss[n],
      //     "     n: ", n,
      //     "     mu_loss: ", mu_loss[n],
      //     "     shape_loss: ", shape_loss[n],
      //     "     Intercept_shape_loss: ", Intercept_shape_loss,
      //     "     b_shape_loss: ", b_shape_loss,
      //     "     b_loss: ", b_loss_s1,
      //   "     lb lcdf: ", gamma_lcdf(lb_loss[n] | shape_loss[n], mu_loss[n]), "________",
      //   "     ub lcdf: ", gamma_lcdf(ub_loss[n] | shape_loss[n], mu_loss[n]),
      //   "    log diff 2:", log_diff_exp(!!{sev_dist}!!_lcdf(ub_!!{sev_resp}!![n] | !!{sev_arg_stan}!!), !!{sev_dist}!!_lcdf(lb_!!{sev_resp}!![n] | !!{sev_arg_stan}!!))
      //   );
      //
      // }

    }

    target += sum(sev_target .* (1 - freq));

    for (n in 1:N_!!{freq_resp}!!) {
    // special treatment of censored data
      if (cens_!!{freq_resp}!![n] == 0) {
        freq_target[n] = !!{freq_dist}!!_lpmf(Y_!!{freq_resp}!![n] | !!{freq_arg_stan}!!) -
        log_diff_exp(!!{freq_dist}!!_lcdf(ub_!!{freq_resp}!![n] | !!{freq_arg_stan}!!), !!{freq_dist}!!_lcdf(lb_!!{freq_resp}!![n] - 1 | !!{freq_arg_stan}!!));
      } else if (cens_!!{freq_resp}!![n] == 1) {
        freq_target[n] = !!{freq_dist}!!_lccdf(Y_!!{freq_resp}!![n] | !!{freq_arg_stan}!!) -
        log_diff_exp(!!{freq_dist}!!_lcdf(ub_!!{freq_resp}!![n] | !!{freq_arg_stan}!!), !!{freq_dist}!!_lcdf(lb_!!{freq_resp}!![n] - 1 | !!{freq_arg_stan}!!));
      } else if (cens_!!{freq_resp}!![n] == -1) {
        freq_target[n] = !!{freq_dist}!!_lcdf(Y_!!{freq_resp}!![n] | !!{freq_arg_stan}!!) -
        log_diff_exp(!!{freq_dist}!!_lcdf(ub_!!{freq_resp}!![n] | !!{freq_arg_stan}!!), !!{freq_dist}!!_lcdf(lb_!!{freq_resp}!![n] - 1 | !!{freq_arg_stan}!!));
      }

    }

    target += sum(freq_target .* freq);

  }

    '

    code_model_template =
      str_glue(
        template_stan_code,
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
