functions{
  #include comp_pmf.stan
  #include infiniteSumFixedNmax.stan
  #include Fixed.stan
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
  int<lower=0> N_cap;
}
parameters{
  real mu;
  real<lower=0> nu;
}
transformed parameters{
  real log_mu = log(mu);
  real log_norm_const[2] = log_Z_COMP_fixed(log_mu, nu, N_cap, M);
}
model{
  mu ~ gamma(s_mu, r_mu);
  nu ~ normal(0, nu_sd);
  // Likelihood
  for(k in 1:K){
    target += n[k] * COM_Poisson_lpmf(y[k] | log_mu, nu, log_norm_const[1]);
  } 
}
generated quantities{
  real n_iter = log_norm_const[2];
}