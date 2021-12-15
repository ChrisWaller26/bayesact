#' Leave One Out Criterion for Bayesact Class
#'
#' @export
#'
loo = function(x, resp, newdata = NULL, sev_samples = NULL, ...){

  if(is.bayesact(x)){

    freq_link = get(x$bayesact$freq_family$link)

    sev_resp = x$bayesact$sev_formula$resp
    freq_resp = x$bayesact$freq_formula$resp

    sev_pars = x$family[[sev_resp]]$dpars
    sev_pars_n = length(sev_pars)

    iter_tot = x$bayesact$chains * (x$bayesact$iter - x$bayesact$warmup)

    freq_sev =
      ifelse(
        sev_resp == resp,
        "sev",
        "freq"
      )

    if(freq_sev == "sev"){

      brms::loo(
        x = x,
        resp = resp,
        newdata = x$data %>% filter(freq == 0),
        ...
      )

    }else{

      new_freq_data =
        x$data %>%
        filter(freq == 1) %>%
        mutate(
          data_row_id = row_number()
        )

      get_surv = function(ded, par1, par2, par3){

        pfun = pfun_map[[x$bayesact$sev_family$family]]

        if(sev_pars_n == 1){

          1 - pfun(ded, par1)

        }else if(sev_pars_n == 2){

          1 - pfun(ded, par1, par2)

        }else{

          1 - pfun(ded, par1, par2, par3)

        }

        }

      if(is.null(sev_samples)){

        sev_samples =
          sample(seq(iter_tot),
                 min(ceiling(1e6 / nrow(new_freq_data)), iter_tot)
                 )

      }

      new_freq_data =
        lapply(
          sev_pars,
          function(par){

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
                  values_to = par
                ) %>%
              mutate(
                data_row_id = as.integer(substr(data_row_id, 2, 1000))
              )

            if(match(par, sev_pars) > 1){

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
          ded_offset =
            pmax(x$bayesact$ded_adj_min,
                 get_surv(get(x$bayesact$ded_name),
                          get(sev_pars[1]),
                          get(sev_pars[2]),
                          get(sev_pars[3])))
        ) %>%
        group_by_at(
          names(new_freq_data)
        ) %>%
        summarise(
          ded_offset = mean(freq_link(ded_offset)),
          .groups = "keep"
        ) %>%
        ungroup()

      new_freq_formula =
        x$formula$forms[[freq_resp]]

      new_freq_formula$pforms$f1[3] =
        str2expression(
            paste(as.character(new_freq_formula$pforms$f1[3]),
                  "+ offset(ded_offset)")
        )

      freq_prior =
        x$prior %>%
        filter(resp == freq_resp) %>%
        mutate(resp = "")

      cat("Recompiling frequency model\n\n")

      if(x$backend == "rstan"){

        x_freq =
          brm(
            formula       = new_freq_formula,
            family        = x$bayesact$freq_family,
            data          = new_freq_data,
            prior         = freq_prior,
            chains        = x$bayesact$chains,
            iter          = x$bayesact$iter,
            warmup        = x$bayesact$warmup,
            sample_prior  = x$bayesact$sample_prior,
            stanvars      = x$bayesact$stanvars,
            control       = x$bayesact$control,
            backend       = x$backend
          )

      }else{

        x_freq =
          brm(
            formula       = new_freq_formula,
            family        = x$bayesact$freq_family,
            data          = new_freq_data,
            prior         = freq_prior,
            chains        = x$bayesact$chains,
            iter          = x$bayesact$iter,
            warmup        = x$bayesact$warmup,
            adapt_delta   = x$bayesact$adapt_delta,
            max_treedepth = x$bayesact$max_treedepth,
            sample_prior  = x$bayesact$sample_prior,
            stanvars      = x$bayesact$stanvars,
            backend       = x$backend
          )

      }

      cat("Running LOO Function...\n\n")

      brms::loo(
        x = x_freq,
        ...
      )

    }

  }else{

    brms::loo(x = x, resp = resp, ...)

  }

}
