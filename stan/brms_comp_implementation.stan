functions{
  // Taken from https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_com_poisson.stan
  // log approximate normalizing constant of the COM poisson distribuion
  // approximation based on doi:10.1007/s10463-017-0629-6
  // Args: see log_Z_com_poisson()
  real log_Z_com_poisson_approx(real log_mu, real nu) {
    real nu_mu = nu * exp(log_mu); 
    real nu2 = nu^2;
    // first 4 terms of the residual series
    real log_sum_resid = log1p(
      nu_mu^(-1) * (nu2 - 1) / 24 + 
      nu_mu^(-2) * (nu2 - 1) / 1152 * (nu2 + 23) +
      nu_mu^(-3) * (nu2 - 1) / 414720 * (5 * nu2^2 - 298 * nu2 + 11237)
    );
    return nu_mu + log_sum_resid  - 
      ((log(2 * pi()) + log_mu) * (nu - 1) / 2 + log(nu) / 2);
  }
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
    if (log_mu * nu >= log(1.5) && log_mu >= log(1.5)) {
      return log_Z_com_poisson_approx(log_mu, nu);
    }
    // direct computation of the truncated series
    // check if the Mth term of the series is small enough
    if (nu * (M * log_mu - lgamma(M + 1)) > -36.0) {
      reject("nu is too close to zero.");
    }
    // first 2 terms of the series
    log_Z = log1p_exp(nu * log_mu);
    while (converged == 0) {
      // adding terms in batches simplifies the AD tape
      vector[num_terms + 1] log_Z_terms;
      int i = 1;
      log_Z_terms[1] = log_Z;
      while (i <= num_terms) {
        log_Z_terms[i + 1] = nu * (k * log_mu - lgamma(k + 1));
        k += 1;
        if (log_Z_terms[i + 1] <= -36.0) {
          converged = 1;
          break;
        }
        i += 1;
      }
      log_Z = log_sum_exp(log_Z_terms[1:(i + 1)]);
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
}
parameters{
  real mu;
  real<lower=0> nu;
}
transformed parameters{
  real log_norm_const = log_Z_com_poisson(log(mu), nu);
}
model{
  mu ~ gamma(s_mu, r_mu);
  nu ~ normal(0, nu_sd);
  // Likelihood
  for(k in 1:K) target += n[k] * com_poisson_lpmf(y[k] | mu, nu, log_norm_const);
}
