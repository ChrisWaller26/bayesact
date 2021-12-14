
library(brms)

#' Distribution Function Mapping
#'
#' @export
dfun_map =
  list(
    lognormal          = dlnorm,
    Gamma              = function(x, scale, shape) dgamma(x, shape = shape, scale = scale),
    frechet            = dfrechet,
    weibull            = dweibull,
    gaussian           = dnorm,
    shifted_lognormal  = dshifted_lnorm,
    skew_normal        = dskew_normal,
    exponential        = dexp,
    gen_extreme_value  = dgen_extreme_value,
    exgaussian         = dexgaussian,
    Beta               = dbeta,
    hurdle_gamma       = dhurdle_gamma,
    hurdle_lognormal   = dhurdle_lognormal,
    zero_inflated_beta = dzero_inflated_beta,
    student            = dstudent_t
  )
#' @export
pfun_map =
  list(
    lognormal          = plnorm,
    Gamma              = function(x, scale, shape) pgamma(x, shape = shape, scale = scale),
    frechet            = pfrechet,
    weibull            = pweibull,
    gaussian           = pnorm,
    shifted_lognormal  = pshifted_lnorm,
    skew_normal        = pskew_normal,
    exponential        = pexp,
    gen_extreme_value  = pgen_extreme_value,
    exgaussian         = pexgaussian,
    Beta               = pbeta,
    hurdle_gamma       = phurdle_gamma,
    hurdle_lognormal   = phurdle_lognormal,
    zero_inflated_beta = pzero_inflated_beta,
    student            = pstudent_t
  )
