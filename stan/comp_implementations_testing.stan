functions{
  // Our implementation of things
  real log_COM_Poisson(int k, real log_mu, real nu){
    return k * log_mu - nu * lgamma(k + 1);
  }
  real log_Z_com_poisson_approx_new(real log_mu, real nu) {
    // Based on equations (4) and (31) of doi:10.1007/s10463-017-0629-6
    real nu2 = nu^2;
    real log_resids[4];
    real ans;
    real lcte = (nu * exp(log_mu/nu)) - ( (nu-1)/(2*nu)* log_mu + (nu-1)/2*log(2*pi()) + 0.5 *log(nu));
    real c_1 = (nu2-1)/24;
    real c_2 = (nu2-1)/1152*(nu2 + 23);
    real c_3 = (nu2-1)/414720* (5*square(nu2) - 298*nu2 + 11237);
    if(nu < 1){
      print("Approximation doesn't work great when nu < 1, returning without residuals");
      return(lcte);
    } 
    log_resids[1] = 0;
    log_resids[2] = log(c_1) - 1 * (log(nu) + log_mu/nu);
    log_resids[3] = log(c_2) - 2 * (log(nu) + log_mu/nu);
    log_resids[4] = log(c_3) - 3 * (log(nu) + log_mu/nu);
    ans = lcte + log_sum_exp(log_resids);
    return ans;
  }
  real log_COM_Poisson_constant(real log_mu, real nu, real Eps, int n0_, int maxIter) {
    vector[maxIter+1] storeVal;
    real leps = log(Eps) + log(2);
    int n = 1;
    int n0 = n0_;
    real old_term;
    real new_term;
    
    // Setting up first iterations
    old_term = log_COM_Poisson(n0, log_mu, nu);
    n0 += 1;
    new_term = log_COM_Poisson(n0, log_mu, nu);
    n0 += 1;
    storeVal[1] = old_term;
    
    // Bulk of the distribution
    while (new_term - old_term > 0) {
      // print("Doing inner loop, n:", n);
      old_term = new_term;
      new_term = log_COM_Poisson(n0, log_mu, nu);
      n0 += 1;
      n += 1;
      storeVal[n] = old_term;
      if (n >= maxIter) return(log_sum_exp(storeVal[:n]));
    }
    
    while ( (old_term - log_diff_exp(0, new_term - old_term)) > leps  )  {
      // print("n:", n,  " Delta:", old_term - log_diff_exp(0, new_term - old_term));
      old_term = new_term;
      new_term = log_COM_Poisson(n0, log_mu, nu);
      n0 += 1;
      n += 1;
      storeVal[n] = old_term;
      if (n >= maxIter) return(log_sum_exp(storeVal[:n]));
    }
    storeVal[n+1] = old_term - log_diff_exp(0, new_term - old_term) - log(2);
    
    return log_sum_exp(storeVal[:(n+1)]);
  }
  // Taken from https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_com_poisson.stan
  // log normalizing constant of the COM Poisson distribution
  // implementation inspired by code of Ben Goodrich
  // improved following suggestions of Sebastian Weber (#892)
  // Args:
  //   log_mu: log location parameter
  //   shape: positive shape parameter
  real log_Z_com_poisson(real log_mu, real nu) {
    real log_Z; 
    int k = 2;
    int M = 10000;
    int converged = 0;
    int num_terms = 50;
    if (nu == 1) {
      return exp(log_mu);
    }
    // nu == 0 or Inf will fail in this parameterization
    if (nu <= 0) {
      reject("nu must be positive");
    }
    if (nu == positive_infinity()) {
      reject("nu must be finite");
    }
    // if (log_mu * nu >= log(1.5) && nu >= 1) {
      //   return log_Z_com_poisson_approx(log_mu, nu);
      // }
      // direct computation of the truncated series
      // check if the Mth term of the series is small enough
      // if (M * log_mu - nu*lgamma(M + 1) > -36.0) {
        //   reject("nu is too close to zero.");
        // }
        // first 2 terms of the series
        log_Z = log1p_exp(log_mu);   
        while (converged == 0) {
          // adding terms in batches simplifies the AD tape
          vector[num_terms + 1] log_Z_terms;
          int i = 1;
          log_Z_terms[1] = log_Z;
          while (i < num_terms) {
            log_Z_terms[i + 1] = k * log_mu - nu*lgamma(k + 1);
            k += 1;
            if (log_Z_terms[i + 1] <= -36.0) {
              converged = 1;
              break;
            }
            i += 1;
          }
          log_Z = log_sum_exp(log_Z_terms[1:i]);
        }
        return log_Z;
  }
  
  /* Comparison stuff */
  real  signum(real x) {
    real ans;
    if(x < 0){
      ans = -1;
    }else{
      if(x == 0){
        ans = 0;
      }else{
        ans = 1;
      }
    }
    return ans;
  }
  real robust_difference(real x, real y){
    real sgn = signum(x-y);
    real m = min({x, y});
    real M = max({x, y});
    return(sgn * exp(log_diff_exp(M, m)));
  }
}
data{
  real log_mu;
  real<lower=0> nu;
  real<lower=0> eps;
  int<lower=0> M;
  real true_value;
}
generated quantities {
  real lZ_approx_new = log_Z_com_poisson_approx_new(log_mu, nu);
  real lZ_brute_force_new = log_COM_Poisson_constant(log_mu, nu, eps, 0, M);
  real lZ_brute_force_brms = log_Z_com_poisson(log_mu, nu);
  real diff_approx_new = robust_difference(true_value, lZ_approx_new);
  real diff_brute_force_new = robust_difference(true_value, lZ_brute_force_new);
  real diff_brute_force_brms = robust_difference(true_value, lZ_brute_force_brms);
}
