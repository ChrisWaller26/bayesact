#' Modified Log-Normal Power Law Distribution
#'
#' @description Density, distribution function and random
#' generation for the log normal power law distribution.
#'
#' This is based on applications to modelling masses of stars
#' (https://iopscience.iop.org/article/10.3847/1538-4357/ab88a6) but is also
#' applicable to insurance losses to model both attrirtional and large losses
#' in one distribution.
#' @export
dlnormpower = function(x, mu = 0, sigma = 1, alpha = 1){

  alpha * exp(alpha * mu + 0.5 * alpha ^ 2 * sigma ^ 2) * x ^ (-(1 + alpha)) *
    pnorm((log(x) - mu) / sigma - alpha * sigma)

}
#'
#' @export
plnormpower = function(q, mu = 0, sigma = 1, alpha = 1){

  plnorm(q, mu, sigma) -
    exp(alpha * mu + 0.5 * alpha ^ 2 * sigma ^ 2) * q ^ (-alpha) *
    pnorm((log(q) - mu) / sigma - alpha * sigma)

}
#' @export
rlnormpower = function(n, mu = 0, sigma = 1, alpha = 1){

  exp(mu + sigma * rnorm(n) - log(runif(n)) / alpha)

}
#' @export
lnormpower_exp = function(mu = 0, sigma = 1, alpha = 1){

  alpha / (alpha - 1) * exp(0.5 * sigma ^ 2 + mu)

}
#' @export
lnormpower_mk = function(k = 1, mu = 0, sigma = 1, alpha = 1){

  alpha / (alpha - k) * exp(0.5 * k ^ 2 * sigma ^ 2 + mu * k)

}
#'
#' LogNormal-Power BRMS Custom Family
#'
#' @export
lnormpower <-
  function(link = "identity",
           link_sigma = "log",
           link_alpha = "log1p"){

    custom_family(
      "lnormpower",
      dpars = c("mu", "sigma", "alpha"),
      links = c(link, link_sigma, link_alpha),
      lb = c(NA, 0, 1),
      type = "real"
    )

  }
#'
#' LogNormal-Power BRMS Custom Family Stan Functions
#'
#' @export
lnormpower_stan_funs <- "
  real lnormpower_lpdf(real y, real mu, real sigma, real alpha) {

    real output;

    output = log(alpha) + (alpha * mu + 0.5 * alpha ^ 2 * sigma ^ 2) + log(y) * (-(1 + alpha)) +
    log(Phi((log(y) - mu) / sigma - alpha * sigma));

    return output;
  }

  real lnormpower_lcdf(real y, real mu, real sigma, real alpha) {

    real output;

    output = log(lognormal_cdf(y, mu, sigma) -
    exp(alpha * mu + 0.5 * alpha ^ 2 * sigma ^ 2) * y ^ (-alpha) *
    Phi((log(y) - mu) / sigma - alpha * sigma));

    return output;
  }

  real lnormpower_cdf(real y, real mu, real sigma, real alpha) {

    real output;

    output = exp(lnormpower_lcdf(y | mu, sigma, alpha));

    return output;

  }

  real lnormpower_lccdf(real y, real mu, real sigma, real alpha) {

    real output;

    output = log(1 - lnormpower_cdf(y, mu, sigma, alpha));

    return output;

  }

  real lnormpower_rng(real mu, real sigma, real alpha) {

    real output;

    output = exp(mu + sigma * normal_rng(0, 1) - log(uniform_rng(0, 1)) / alpha);

    return output;

  }

"
#'
#' LogNormal-Power BRMS Custom Family Stan Variables
#'
#' @export
lnormpower_stanvars <- brms::stanvar(scode = lnormpower_stan_funs, block = "functions")

#'
#' LogNormal-Power Expectation
#'
#' @export
posterior_epred_lnormpower <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  alpha <- brms::get_dpar(prep, "alpha")
  dims = dim(mu)
  lb <-
    case_when(
      is.null(prep$data$lb) ~ 0,
      TRUE ~ prep$data$lb
    )
  ub <-
    case_when(
      is.null(prep$data$ub) ~ Inf,
      TRUE ~ prep$data$ub
    )

  matrix(lnormpower_exp(c(mu), c(sigma), c(alpha)),
         nrow = dims[1],
         ncol = dims[2])
}

#'
#' LogNormal-Power Predictions
#'
#' @export
posterior_predict_lnormpower <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  alpha <- brms::get_dpar(prep, "alpha", i = i)
  lb <-
    case_when(
      is.null(prep$data$lb[i]) ~ 0,
      TRUE ~ prep$data$lb[i]
    )
  ub <-
    case_when(
      is.null(prep$data$ub[i]) ~ Inf,
      TRUE ~ prep$data$ub[i]
    )
  rlnormpower(length(mu), mu, sigma, alpha)
}

#'
#' LogNormal-Power Log-Likelihood
#'
#' @export
log_lik_lnormpower <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  alpha <- brms::get_dpar(prep, "alpha", i = i)
  lb <-
    case_when(
      is.null(prep$data$lb[i]) ~ 0,
      TRUE ~ prep$data$lb[i]
    )
  ub <-
    case_when(
      is.null(prep$data$ub[i]) ~ Inf,
      TRUE ~ prep$data$ub[i]
    )
  y <- prep$data$Y[i]
  log(dlnormpower(y, mu, sigma, alpha)) -
    log(plnormpower(ub, mu, sigma, alpha) -
          plnormpower(lb, mu, sigma, alpha))
}
