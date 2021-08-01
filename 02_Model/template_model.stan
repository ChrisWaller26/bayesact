  for(n in 1:N_!!{sev_resp}!!){
    
    mu_!!{sev_resp}!![n] = !!{sev_formula_stan}!!;
    
  }
  
  for(n in 1:N_!!{freq_resp}!!){
    
    mu_!!{freq_resp}!![n] = 
      !!{freq_formula_stan}!!  * 
      (1 - !!{sev_dist}!!_cdf(ded[n], mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n]));
    
  }
  
  for (n in 1:N_!!{sev_resp}!!) {
      // apply the inverse link function
      !!{sev_arg[2]}!!_!!{sev_resp}!![n] = !!{sev_inv_link[2]}!!(!!{sev_arg[2]}!!_!!{sev_resp}!![n]);
    }

  for (n in 1:N_!!{sev_resp}!!) {
    // special treatment of censored data
      if (cens_!!{sev_resp}!![n] == 0) {
        sev_target[n] = !!{sev_dist}!!_lpdf(Y_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n]) -
        log_diff_exp(!!{sev_dist}!!_lcdf(ub_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n]), !!{sev_dist}!!_lcdf(lb_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n]));
      } else if (cens_!!{sev_resp}!![n] == 1) {
        sev_target[n] = !!{sev_dist}!!_lccdf(Y_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n]) -
        log_diff_exp(!!{sev_dist}!!_lcdf(ub_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n]), !!{sev_dist}!!_lcdf(lb_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n]));
      } else if (cens_!!{sev_resp}!![n] == -1) {
        sev_target[n] = !!{sev_dist}!!_lcdf(Y_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n]) -
        log_diff_exp(!!{sev_dist}!!_lcdf(ub_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n]), !!{sev_dist}!!_lcdf(lb_!!{sev_resp}!![n] | mu_!!{sev_resp}!![n], !!{sev_arg[2]}!!_!!{sev_resp}!![n]));
      }
    }
    
    target += sum(sev_target .* (1 - freq));
    
    for (n in 1:N_!!{freq_resp}!!) {
    // special treatment of censored data
      if (cens_!!{freq_resp}!![n] == 0) {
        freq_target[n] = !!{freq_dist}!!_lpmf(Y_!!{freq_resp}!![n] | mu_!!{freq_resp}!![n]) -
        log_diff_exp(!!{freq_dist}!!_lcdf(ub_!!{freq_resp}!![n] | mu_!!{freq_resp}!![n]), !!{freq_dist}!!_lcdf(lb_!!{freq_resp}!![n] - 1 | mu_!!{freq_resp}!![n]));
      } else if (cens_!!{freq_resp}!![n] == 1) {
        freq_target[n] = !!{freq_dist}!!_lccdf(Y_!!{freq_resp}!![n] | mu_!!{freq_resp}!![n]) -
        log_diff_exp(!!{freq_dist}!!_lcdf(ub_!!{freq_resp}!![n] | mu_!!{freq_resp}!![n]), !!{freq_dist}!!_lcdf(lb_!!{freq_resp}!![n] - 1 | mu_!!{freq_resp}!![n]));
      } else if (cens_!!{freq_resp}!![n] == -1) {
        freq_target[n] = !!{freq_dist}!!_lcdf(Y_!!{freq_resp}!![n] | mu_!!{freq_resp}!![n]) -
        log_diff_exp(!!{freq_dist}!!_lcdf(ub_!!{freq_resp}!![n] | mu_!!{freq_resp}!![n]), !!{freq_dist}!!_lcdf(lb_!!{freq_resp}!![n] - 1 | mu_!!{freq_resp}!![n]));
      }
    }
    
    target += sum(freq_target .* freq);
  
  
  
  
  
  
  }