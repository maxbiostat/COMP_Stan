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
    if (log_mu * nu >= log(1.5) && nu >= 1.0) {
      return log_Z_com_poisson_approx_new(log_mu, nu);
    }
    // direct computation of the truncated series
    // check if the Mth term of the series is small enough
    if ( M * log_mu - nu*lgamma(M + 1) > leps) {
      reject("nu is too close to zero.");
    }
    // first 2 terms of the series
    log_Z = log1p_exp(nu * log_mu);
    while (converged == 0) {
      if(k >= M) break;
      // adding terms in batches simplifies the AD tape
      vector[num_terms + 1] log_Z_terms;
      int i = 1;
      log_Z_terms[1] = log_Z;
      while (i <= num_terms) {
        log_Z_terms[i + 1] = k * log_mu - nu*lgamma(k + 1);
        k += 1;
        if (log_Z_terms[i + 1] <= leps) {
          converged = 1;
          break;
        }
        i += 1;
      }
      log_Z = log_sum_exp(log_Z_terms[1:i]);
    }
    return log_Z;
  }
  // COM Poisson log-PMF for a single response (log parameterization)
  // Args: 
  //   y: the response value 
  //   log_mu: log location parameter
  //   shape: positive shape parameter
  real com_poisson_log_lpmf(int y, real log_mu, real nu, real logZ) {
    if (nu == 1) return poisson_log_lpmf(y | log_mu);
    return y * log_mu - nu*lgamma(y + 1) - logZ;
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
