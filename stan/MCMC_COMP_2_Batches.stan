functions{
  #include comp_2_pmf.stan
  #include InfiniteBatches.stan
}
data{
  int<lower=0> K;
  int<lower=0> n[K];
  int<lower=0> y[K];
  int<lower=0> N;
  real<lower=0> s_mu;
  real<lower=0> r_mu;
  real<lower=0> s_nu;
  real<lower=0> r_nu;
  real<lower=0> eps;
  int<lower=0> batchSize;
  int<lower=0> M;
}
parameters{
  real<lower=0> mu;
  real<lower=0> nu;
}
transformed parameters{
  real log_mu = log(mu);
  real log_norm_const[2] = infiniteBatches({log_mu, nu}, batchSize, eps, M, 0);
}
model{
  mu ~ gamma(s_mu, r_mu);
  nu ~ gamma(s_nu, r_nu); // Benson & Friel (2021)
  // Likelihood
  for(k in 1:K){
  target += n[k] * COM_Poisson_2_lpmf(y[k] | log_mu, nu, log_norm_const[1]);
  } 
}
generated quantities{
  real n_iter = log_norm_const[2];
}