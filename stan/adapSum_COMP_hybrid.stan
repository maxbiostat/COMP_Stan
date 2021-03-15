functions{
  // Taken from https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_com_poisson.stan
  // log approximate normalizing constant of the COM poisson distribuion
  // approximation based on doi:10.1007/s10463-017-0629-6
  // Args: see log_Z_com_poisson()
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
      reject("Approximation doesn't work great when nu < 1, returning without residuals");
    } 
    log_resids[1] = 0;
    log_resids[2] = log(c_1) - 1 * (log(nu) + log_mu/nu);
    log_resids[3] = log(c_2) - 2 * (log(nu) + log_mu/nu);
    log_resids[4] = log(c_3) - 3 * (log(nu) + log_mu/nu);
    ans = lcte + log_sum_exp(log_resids);
    return ans;
  }
  real log_COM_Poisson(int k, real log_mu, real nu){
    return k * log_mu - nu * lgamma(k + 1);
  }
  real log_COMP_constant_adaptive(real log_mu, real nu, real Eps, int n0_, int maxIter) {
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
  // log normalizing constant of the COM Poisson distribution
  // implementation inspired by code of Ben Goodrich
  // improved following suggestions of Sebastian Weber (#892)
  // Args:
  //   log_mu: log location parameter
  //   shape: positive shape parameter
  real log_Z_com_poisson(real log_mu, real nu, int M, int num_terms, real eps) {
    real log_Z;
    real leps = log(eps);
    int k = 2;
    int converged = 0;
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
    if (log_mu * nu >= log(1.5) && nu > 1) {
      return log_Z_com_poisson_approx_new(log_mu, nu);
    }
    // direct computation of the truncated series
    // check if the Mth term of the series is small enough
    if ( ( M * log_mu - nu*lgamma(M + 1) ) > leps) {
      reject("nu is too close to zero.");
    }
    log_Z = log_COMP_constant_adaptive(log_mu, nu,  eps, 0, M);
    return log_Z;
  }
  // COM Poisson log-PMF for a single response (log parameterization)
  // Args: 
  //   y: the response value 
  //   log_mu: log location parameter
  //   shape: positive shape parameter
  real com_poisson_log_lpmf(int y, real log_mu, real nu, real logZ) {
    if (nu == 1) return poisson_log_lpmf(y | log_mu);
    return nu * (y * log_mu - lgamma(y + 1)) - logZ;
  }
  // COM Poisson log-PMF for a single response
  real com_poisson_lpmf(int y, real mu, real nu, real logZ) {
    if (nu == 1) return poisson_lpmf(y | mu);
    return com_poisson_log_lpmf(y | log(mu), nu, logZ);
  }
}
data{
  int<lower=0> K;
  int<lower=0> n[K];
  int<lower=0> y[K];
  int<lower=0> N;
  real<lower=0> s_mu;
  real<lower=0> r_mu;
  real<lower=0> nu_sd;
  int<lower=0> M;
  real<lower=0> eps;
  int<lower=0> batch_size;
}
parameters{
  real mu;
  real<lower=0> nu;
}
transformed parameters{
  real log_norm_const = log_Z_com_poisson(log(mu), nu, M, batch_size, eps);
}
model{
  mu ~ gamma(s_mu, r_mu);
  nu ~ normal(0, nu_sd);
  // Likelihood
  for(k in 1:K) target += n[k] * com_poisson_lpmf(y[k] | mu, nu, log_norm_const);
}
