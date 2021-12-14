#' @export
#'
update = function(object,
                  freq_formula  = NULL,
                  sev_formula   = NULL,
                  freq_family   = NULL,
                  sev_family    = NULL,
                  freq_data     = NULL,
                  sev_data      = NULL,
                  ded_name      = NULL,
                  freq_adj_fun  = NULL,
                  mle           = NULL,
                  ded_adj_min   = NULL,
                  prior         = NULL,
                  chains        = NULL,
                  iter          = NULL,
                  warmup        = NULL,
                  refresh       = NULL,
                  adapt_delta   = NULL,
                  max_treedepth = NULL,
                  sample_prior  = NULL,
                  stanvars      = NULL,
                  control       = NULL,
                  backend       = NULL,
                  save_pars      = NULL,
                  ...){

  if(is.bayesact(object)){

    brms_freq_sev(

      freq_formula =
        if(is.null(freq_formula)){
               object$bayesact$freq_formula
          }else{
               freq_formula},

      sev_formula =
        if(is.null(sev_formula)){
               object$bayesact$sev_formula
          }else{
               sev_formula},

      freq_family =
        if(is.null(freq_family)){
               object$bayesact$freq_family
          }else{
               freq_family},
      sev_family =
        if(is.null(sev_family)){
               object$bayesact$sev_family
          }else{
               sev_family},

      freq_data =
        if(is.null(freq_data)){
               object$bayesact$freq_data
          }else{
               freq_data},
      sev_data =
        if(is.null(sev_data)){
               object$bayesact$sev_data
          }else{
               sev_data},

      prior =
        if(is.null(prior)){
               object$bayesact$prior
          }else{
               prior},

      ded_name =
        if(is.null(ded_name)){
               object$bayesact$ded_name
          }else{
               ded_name},
      ded_adj_min =
        if(is.null(ded_adj_min)){
               object$bayesact$ded_adj_min
          }else{
               ded_adj_min},
      backend =
        if(is.null(backend)){
               object$backend
          }else{
               backend},

      chains =
        if(is.null(chains)){
               object$bayesact$chains
          }else{
               chains},
      iter =
        if(is.null(iter)){
               object$bayesact$iter
          }else{
               iter},
      warmup =
        if(is.null(warmup)){
               object$bayesact$warmup
          }else{
               warmup},

      refresh =
        if(is.null(refresh)){
               object$bayesact$refresh
          }else{
               refresh},
      adapt_delta =
        if(is.null(adapt_delta)){
               object$bayesact$adapt_delta
          }else{
               adapt_delta},
      max_treedepth =
        if(is.null(max_treedepth)){
               object$bayesact$max_treedepth
          }else{
               max_treedepth},

      mle =
        if(is.null(mle)){
               object$bayesact$mle
          }else{
               mle},
      sample_prior =
        if(is.null(sample_prior)){
               object$bayesact$sample_prior
          }else{
               sample_prior},
      freq_adj_fun =
        if(is.null(freq_adj_fun)){
               object$bayesact$freq_adj_fun
          }else{
               freq_adj_fun},
      stanvars     =
        if(is.null(stanvars)){
               object$bayesact$stanvars
          }else{
               stanvars},

      save_pars =
        if(is.null(save_pars)){
               object$save_pars}
      else{
               save_pars},
      ...
    )

  }else{

    stats::update(object, ...)

  }

}
