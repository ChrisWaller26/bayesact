
  for(n in 1:N_!!{freq_resp}!!){
    
    mu_!!{freq_resp}!![n] = 
      mu_!!{freq_resp}!![n] * 
       (1 - !!{sev_dist}!!_cdf(ded[n], !!{sev_arg_stan}!!));
    
  }
  
  for (n in 1:N_!!{sev_resp}!!) {
    // special treatment of censored data
      if (cens_!!{sev_resp}!![n] == 0) {
        sev_target[n] = !!{sev_dist}!!_lpdf(Y_!!{sev_resp}!![n] | !!{sev_arg_stan}!!) -
        log_diff_exp(!!{sev_dist}!!_lcdf(ub_!!{sev_resp}!![n] |!!{sev_arg_stan}!!), !!{sev_dist}!!_lcdf(lb_!!{sev_resp}!![n] | !!{sev_arg_stan}!!));
      } else if (cens_!!{sev_resp}!![n] == 1) {
        sev_target[n] = !!{sev_dist}!!_lccdf(Y_!!{sev_resp}!![n] | !!{sev_arg_stan}!!) -
        log_diff_exp(!!{sev_dist}!!_lcdf(ub_!!{sev_resp}!![n] | !!{sev_arg_stan}!!), !!{sev_dist}!!_lcdf(lb_!!{sev_resp}!![n] | !!{sev_arg_stan}!!));
      } else if (cens_!!{sev_resp}!![n] == -1) {
        sev_target[n] = !!{sev_dist}!!_lcdf(Y_!!{sev_resp}!![n] | !!{sev_arg_stan}!!) -
        log_diff_exp(!!{sev_dist}!!_lcdf(ub_!!{sev_resp}!![n] | !!{sev_arg_stan}!!), !!{sev_dist}!!_lcdf(lb_!!{sev_resp}!![n] | !!{sev_arg_stan}!!));
      }
      
      // if(is_nan(sev_target[n]) || is_inf(sev_target[n])){
      // 
      //   print(
      //     "Y_loss: ", Y_loss[n],
      //     "     n: ", n,
      //     "     mu_loss: ", mu_loss[n],
      //     "     shape_loss: ", shape_loss[n],
      //     "     Intercept_shape_loss: ", Intercept_shape_loss,
      //     "     b_shape_loss: ", b_shape_loss,
      //     "     b_loss: ", b_loss_s1,
      //   "     lb lcdf: ", gamma_lcdf(lb_loss[n] | shape_loss[n], mu_loss[n]), "________",
      //   "     ub lcdf: ", gamma_lcdf(ub_loss[n] | shape_loss[n], mu_loss[n]),
      //   "    log diff 2:", log_diff_exp(!!{sev_dist}!!_lcdf(ub_!!{sev_resp}!![n] | !!{sev_arg_stan}!!), !!{sev_dist}!!_lcdf(lb_!!{sev_resp}!![n] | !!{sev_arg_stan}!!))
      //   );
      //   
      // }
      
    }
    
    target += sum(sev_target .* (1 - freq));
    
    for (n in 1:N_!!{freq_resp}!!) {
    // special treatment of censored data
      if (cens_!!{freq_resp}!![n] == 0) {
        freq_target[n] = !!{freq_dist}!!_lpmf(Y_!!{freq_resp}!![n] | !!{freq_arg_stan}!!) -
        log_diff_exp(!!{freq_dist}!!_lcdf(ub_!!{freq_resp}!![n] | !!{freq_arg_stan}!!), !!{freq_dist}!!_lcdf(lb_!!{freq_resp}!![n] - 1 | !!{freq_arg_stan}!!));
      } else if (cens_!!{freq_resp}!![n] == 1) {
        freq_target[n] = !!{freq_dist}!!_lccdf(Y_!!{freq_resp}!![n] | !!{freq_arg_stan}!!) -
        log_diff_exp(!!{freq_dist}!!_lcdf(ub_!!{freq_resp}!![n] | !!{freq_arg_stan}!!), !!{freq_dist}!!_lcdf(lb_!!{freq_resp}!![n] - 1 | !!{freq_arg_stan}!!));
      } else if (cens_!!{freq_resp}!![n] == -1) {
        freq_target[n] = !!{freq_dist}!!_lcdf(Y_!!{freq_resp}!![n] | !!{freq_arg_stan}!!) -
        log_diff_exp(!!{freq_dist}!!_lcdf(ub_!!{freq_resp}!![n] | !!{freq_arg_stan}!!), !!{freq_dist}!!_lcdf(lb_!!{freq_resp}!![n] - 1 | !!{freq_arg_stan}!!));
      }
   
    }
    
    target += sum(freq_target .* freq);
    
  }