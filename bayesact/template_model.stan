
  for(n in 1:N_!!{freq_resp}!!){

    mu_!!{freq_resp}!![n] =
      mu_!!{freq_resp}!![n] *
       (1 - !!{sev_dist}!!_cdf(ded[n], !!{sev_arg_stan}!!));

  }
