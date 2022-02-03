#' Bayes Factors for Bayesact Class
#'
#' @export
#'
bayes_factor = function(x1, x2, resp = NULL, newdata = NULL, sev_samples = NULL, sample_max = 1e6, custom_pfun = NULL,...){

  if(is.null(resp)){

    stop("Response Variable Required")

  }

  if(is.bayesact(x1) != is.bayesact(x2)){

    stop("Both models should be of bayesact or non-bayesact class")

  }

  if(x1$bayesact$sev_formula$resp  != x2$bayesact$sev_formula$resp |
     x1$bayesact$freq_formula$resp != x2$bayesact$freq_formula$resp){

    stop("Both models should have the same response variable")

  }

  if(!identical(x1$bayesact$freq_data, x2$bayesact$freq_data) |
     !identical(x1$bayesact$sev_data, x2$bayesact$sev_data)){

    stop("Both models should use the same data")

  }

  if(is.bayesact(x1)){

    freq_link1 = get(x1$bayesact$freq_family$link)
    freq_link2 = get(x2$bayesact$freq_family$link)

    sev_resp = x1$bayesact$sev_formula$resp
    freq_resp = x1$bayesact$freq_formula$resp

    sev_pars1 = x1$family[[sev_resp]]$dpars
    sev_pars_n1 = length(sev_pars1)
    sev_pars2 = x2$family[[sev_resp]]$dpars
    sev_pars_n2 = length(sev_pars2)

    iter_tot = x1$bayesact$chains * (x1$bayesact$iter - x1$bayesact$warmup)

    freq_sev =
      ifelse(
        sev_resp == resp,
        "sev",
        "freq"
      )

    if(freq_sev == "sev"){

      sev_prior1 =
        x1$prior %>%
        filter(resp == sev_resp) %>%
        mutate(resp = "")
      sev_prior2 =
        x2$prior %>%
        filter(resp == sev_resp) %>%
        mutate(resp = "")

      x1_sev =
        brm(
          formula       = x1$formula$forms[[sev_resp]],
          family        = x1$bayesact$sev_family,
          data          = x1$data %>% filter(freq == 0),
          prior         = sev_prior1,
          chains        = x1$bayesact$chains,
          iter          = x1$bayesact$iter,
          warmup        = x1$bayesact$warmup,
          sample_prior  = x1$bayesact$sample_prior,
          stanvars      = x1$bayesact$stanvars,
          control       = x1$bayesact$control,
          backend       = "rstan",
          save_pars     = x1$save_pars
        )

      x2_sev =
        brm(
          formula       = x2$formula$forms[[sev_resp]],
          family        = x2$bayesact$sev_family,
          data          = x2$data %>% filter(freq == 0),
          prior         = sev_prior2,
          chains        = x2$bayesact$chains,
          iter          = x2$bayesact$iter,
          warmup        = x2$bayesact$warmup,
          sample_prior  = x2$bayesact$sample_prior,
          stanvars      = x2$bayesact$stanvars,
          control       = x2$bayesact$control,
          backend       = "rstan",
          save_pars     = x2$save_pars
        )

      return(
        brms::bayes_factor(
          x1 = x1_sev,
          x2 = x2_sev,
          ...
        )
      )

    }else{

      new_freq_data =
        x1$data %>%
        filter(freq == 1) %>%
        mutate(
          data_row_id = row_number()
        )

      get_surv = function(x, sev_pars_n, ded, par1, par2, par3, par4, par5){

        if(x$bayesact$sev_family$family == "custom"){

          if(is.null(custom_pfun)){
            pfun = pfun_map[[x$bayesact$sev_family$name]]
          }else{
            pfun = custom_pfun
          }

          if(is.null(pfun)){
            stop("Must specify custom_pfun for custom families")
          }

        }else{
          pfun = pfun_map[[x$bayesact$sev_family$family]]
        }

        if(sev_pars_n == 1){
          1 - pfun(ded, par1)
        }else if(sev_pars_n == 2){
          1 - pfun(ded, par1, par2)
        }else if(sev_pars_n == 3){
          1 - pfun(ded, par1, par2, par3)
        }else if(sev_pars_n == 4){
          1 - pfun(ded, par1, par2, par3, par4)
        }else{
          1 - pfun(ded, par1, par2, par3, par4, par5)
        }

        }

      if(is.null(sev_samples)){

        sev_samples =
          sample(seq(iter_tot),
                 min(ceiling(sample_max / nrow(new_freq_data)), iter_tot)
                 )

      }

      sev_pars = c(sev_pars1, sev_pars2)
      x_i = rep(1:2, c(sev_pars_n1, sev_pars_n2))

      new_freq_data =
        lapply(
          seq(length(x_i)),
          function(i){

            par = sev_pars[[i]]
            x = list(x1, x2)[[x_i[i]]]

            output =
              posterior_epred(
                x,
                resp = sev_resp,
                dpar = par,
                newdata = new_freq_data,
                draw_ids = sev_samples
              ) %>%
                as.data.frame() %>%
                mutate(
                  iter_number = row_number()
                ) %>%
                pivot_longer(
                  cols = -iter_number,
                  names_to = "data_row_id",
                  values_to = paste0(par, x_i[i])
                ) %>%
              mutate(
                data_row_id = as.integer(substr(data_row_id, 2, 1000))
              )

            if(i > 1){

              output =
                output %>%
                  select(
                    -iter_number,
                    -data_row_id
                  )

            }

            output

          }
        ) %>%
        bind_cols() %>%
        left_join(
          new_freq_data,
          by = "data_row_id"
        ) %>%
        mutate(
          ded_offset1 =
            pmax(x1$bayesact$ded_adj_min,
                 get_surv(x1,
                          sev_pars_n1,
                          get(x1$bayesact$ded_name),
                          get(paste0(sev_pars1[1], 1)),
                          get(paste0(sev_pars1[2], 1)),
                          get(paste0(sev_pars1[3], 1)),
                          get(paste0(sev_pars1[4], 1)),
                          get(paste0(sev_pars1[5], 1)))),

          ded_offset2 =
            pmax(x2$bayesact$ded_adj_min,
                 get_surv(x2,
                          sev_pars_n2,
                          get(x2$bayesact$ded_name),
                          get(paste0(sev_pars2[1], 2)),
                          get(paste0(sev_pars2[2], 2)),
                          get(paste0(sev_pars2[3], 2)),
                          get(paste0(sev_pars2[4], 2)),
                          get(paste0(sev_pars2[5], 2))))
        ) %>%
        group_by_at(
          names(new_freq_data)
        ) %>%
        summarise(
          ded_offset1 = mean(freq_link1(ded_offset1)),
          ded_offset2 = mean(freq_link2(ded_offset2)),
          .groups = "keep"
        ) %>%
        ungroup()

      new_freq_formula1 =
        x1$formula$forms[[freq_resp]]
      new_freq_formula2 =
        x2$formula$forms[[freq_resp]]

      new_freq_formula1$pforms[[1]][3] =
        str2expression(
          paste(as.character(new_freq_formula1$pforms[[1]][3]),
                "+ offset(ded_offset1)")
        )
      new_freq_formula2$pforms[[1]][3] =
        str2expression(
          paste(as.character(new_freq_formula2$pforms[[1]][3]),
                "+ offset(ded_offset2)")
        )

      freq_prior1 =
        x1$prior %>%
        filter(resp == freq_resp) %>%
        mutate(resp = "")
      freq_prior2 =
        x2$prior %>%
        filter(resp == freq_resp) %>%
        mutate(resp = "")

      cat("Recompiling frequency model 1\n\n")

      x1_freq =
        brm(
          formula       = new_freq_formula1,
          family        = x1$bayesact$freq_family,
          data          = new_freq_data,
          prior         = freq_prior1,
          chains        = x1$bayesact$chains,
          iter          = x1$bayesact$iter,
          warmup        = x1$bayesact$warmup,
          sample_prior  = x1$bayesact$sample_prior,
          stanvars      = x1$bayesact$stanvars,
          control       = x1$bayesact$control,
          backend       = "rstan",
          save_pars     = x1$save_pars
        )

      cat("Recompiling frequency model 2\n\n")

      x2_freq =
        brm(
          formula       = new_freq_formula2,
          family        = x2$bayesact$freq_family,
          data          = new_freq_data,
          prior         = freq_prior2,
          chains        = x2$bayesact$chains,
          iter          = x2$bayesact$iter,
          warmup        = x2$bayesact$warmup,
          sample_prior  = x2$bayesact$sample_prior,
          stanvars      = x2$bayesact$stanvars,
          control       = x2$bayesact$control,
          backend       = "rstan",
          save_pars     = x2$save_pars
        )

      cat("Running Bayes Factor Function...\n\n")

      return(
        brms::bayes_factor(
          x1 = x1_freq,
          x2 = x2_freq,
          ...
        )
      )

    }

  }else{

    return(
      brms::bayes_factor(
        x1 = x1,
        x2 = x2,
        resp = resp,
        ...)
    )

  }

}
