functions{
  #include comp_2_pmf.stan
  #include infiniteSumFixedNmax.stan
  #include fixed_2.stan
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
  int<lower=0> N_max;
}
parameters{
  real mu;
  real<lower=0> nu;
}
transformed parameters{
  real log_theta = log(mu)/nu;
  real log_norm_const[2] = log_Z_COMP_2_fixed(log_theta, nu, N_max, M);
}
model{
  mu ~ gamma(s_mu, r_mu);
  nu ~ normal(0, nu_sd);
  // Likelihood
  for(k in 1:K){
   target += n[k] * COM_Poisson_2_lpmf(y[k] | log_theta, nu, log_norm_const[1]);
  } 
}
generated quantities{
  real theta = exp(log_theta);
  real n_iter = log_norm_const[2];
}