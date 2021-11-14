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
