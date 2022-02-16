functions{

  real lognormal_poisson_lpdf(real y, real lambda, real mu, real sigma){

    int N_COUNT = 20;
    real output[N_COUNT + 1] = rep_vector(0, N_COUNT + 1);

    if(y == 0){

      output[1] = poisson_lpmf(0 | lambda);

    }else{

      for(i in 1:N_COUNT){

        if(i == 1){

          output[2] =
            poisson_lpmf(1 | lambda) +
            lognormal_lpdf(y | mu, sigma);

        }

      }

    }

    return sum(output);

  }

}

data {
  int<lower=0> N;xc
  vector[N] loss_agg;
}

parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=0> lambda;
}

model {

  for(i in 1:N){
    target += lognormal_poisson_lpdf(loss_agg[i], lambda, mu, sigma);
  }

}

