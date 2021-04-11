/* Implements a direct "naive" strategy that stores values until term < eps and */
/* then sums them in bulk */
real [] log_Z_COMP_naive(real log_mu, real nu, real Eps, int maxIter) {
  vector[maxIter+2] storeVal;
  real leps = log(Eps);
  int n = 3;
  
  storeVal[1] = 0;
  storeVal[2] = log_mu;
  
  // Setting up first iterations
  while (n < exp(log_mu / nu)){
    storeVal[n] = log_COM_Poisson(n - 1, log_mu, nu);
    n += 1;
    if(n >= maxIter) break;
  }
  storeVal[n] = log_COM_Poisson(n - 1, log_mu, nu);
  
  while (storeVal[n] > leps){
    n += 1;
    storeVal[n] = log_COM_Poisson(n - 1, log_mu, nu);
    if(n >= maxIter) break;
  }
  return {log_sum_exp(storeVal[:n]), n};
}

