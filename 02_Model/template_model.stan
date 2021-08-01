  for(n in 1:N_!!{sev_resp}!!){
    
    mu_!!{sev_resp}!![n] = !!{sev_formula_stan}!!;
    
  }
  
  for(n in 1:N_!!{sev_resp}!!){
    
    mu_!!{freq_resp}!![n] = !!{freq_formula_stan}!!;
    
  }
  
  for (n in 1:N_!!{sev_resp}!!) {
      // apply the inverse link function
      !!{sev_arg[2]}!!_!!{sev_resp}!![n] = !!{sev_inv_link[2]}!!(!!{sev_arg[2]}!!_!!{sev_resp}!![n]);
    }
    
  // Severity
  
  for(n in 1:N_!!{sev_resp}!!){
    
      if (cens_!!{sev_resp}!![n] == 0) {
        target += (!!{sev_dist}!!_lpdf(Y_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n]) -
        !!{sev_dist}!!_lccdf(lb_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n])) * (1 - freq[n]);
      } else if (cens_!!{sev_resp}!![n] == 1) {
        target += (!!{sev_dist}!!_lccdf(Y_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n]) -
        !!{sev_dist}!!_lccdf(lb_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n])) * (1 - freq[n]);
      } else if (cens_!!{sev_resp}!![n] == -1) {
        target += (!!{sev_dist}!!_lcdf(Y_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n]) -
        !!{sev_dist}!!_lccdf(lb_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n])) * (1 - freq[n]);
      }
    }
    
  for(n in 1:N_!!{freq_resp}!!){
    
    target += 
      freq[n] * 
        !!{freq_dist}!!_lpmf(Y_!!{freq_resp}!![n]| 
                      mu_!!{freq_resp}!![n] * (1 - !!{sev_dist}!!_cdf(ded[n], mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n])));
    
  }
  
  }