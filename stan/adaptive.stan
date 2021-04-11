/*Adaptive algorithm -- direct implementation */
real [] log_Z_COMP_adaptive(real log_mu, real nu, real Eps, int maxIter, int n0_) {
  vector[maxIter+2] storeVal;
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
  while (new_term > old_term) {
    // print("Doing inner loop, n:", n);
    old_term = new_term;
    new_term = log_COM_Poisson(n0, log_mu, nu);
    n0 += 1;
    n += 1;
    storeVal[n] = old_term;
    if (n >= maxIter) return({log_sum_exp(storeVal[:n]), n});
  }
  // Tail of the distribution
  while ( (new_term - log_diff_exp(0, new_term - old_term)) > leps  )  {
    // print("n:", n,  " Delta:", old_term - log_diff_exp(0, new_term - old_term));
    old_term = new_term;
    new_term = log_COM_Poisson(n0, log_mu, nu);
    n0 += 1;
    n += 1;
    storeVal[n] = old_term;
    if (n >= maxIter) return({log_sum_exp(storeVal[:n]), n});
  }
  //storeVal[n+1] = new_term;
  storeVal[n+1] = new_term - log_diff_exp(0, new_term - old_term) - log(2);
  storeVal[n+2] = new_term - log(2);
  return {log_sum_exp(storeVal[:(n+2)]), n};
}

