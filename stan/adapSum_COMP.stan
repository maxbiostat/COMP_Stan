functions{
  // Taken from https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_com_poisson.stan
  real log_COM_Poisson(int k, real log_mu, real nu){
    return k * log_mu - nu * lgamma(k + 1);
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
  real<lower=0> eps;
  int<lower=0> M;
}
parameters{
  real mu;
  real<lower=0> nu;
}
transformed parameters{
  real log_norm_const = log_COM_Poisson_constant(log(mu), nu, eps, 0, M);
}
model{
  mu ~ gamma(s_mu, r_mu);
  nu ~ normal(0, nu_sd);
  // Likelihood
  for(k in 1:K) target += n[k] * com_poisson_lpmf(y[k] | mu, nu, log_norm_const);
}
