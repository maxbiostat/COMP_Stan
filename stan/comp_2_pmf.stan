/* Conway-Maxwell Poisson probability mass function */ 
real log_COM_Poisson_2(int k, real log_mu, real nu){
    return nu*(k * log_mu - lgamma(k + 1));
}
real COM_Poisson_2_lpmf(int k, real log_mu, real nu, real logZ){
    return nu*(k * log_mu - lgamma(k + 1)) - logZ;
}
real logFunction(int n, real[] p){
  return(log_COM_Poisson_2(n, p[1], p[2]));
}