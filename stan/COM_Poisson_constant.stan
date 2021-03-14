functions{
  real log_COM_Poisson(int k, real log_mu, real nu){
    return k * log_mu - nu * lgamma(k + 1);
  }

  real log_COM_Poisson_constant(real mu, real nu, real Eps, int n0_, int maxIter) {
    vector[maxIter+1] storeVal;
    real lmu = log(mu);
    real leps = log(Eps) + log(2);
    int n = 1;
    int n0 = n0_;
    real old_term;
    real new_term;

    // Setting up first iterations
    old_term = log_COM_Poisson(n0, lmu, nu);
    n0 += 1;
    new_term = log_COM_Poisson(n0, lmu, nu);
    n0 += 1;
    storeVal[1] = old_term;

    // Bulk of the distribution
    while (new_term - old_term > 0 || n < maxIter) {
      old_term = new_term;
      new_term = log_COM_Poisson(n0, lmu, nu);
      n0 += 1;
      n += 1;
      storeVal[n] = old_term;
    }

    if (n == maxIter) return(log_sum_exp(storeVal[:n]));

    while ( ((old_term - log_diff_exp(0, new_term - old_term)) > leps) || (n < maxIter))  {
      old_term = new_term;
      new_term = log_COM_Poisson(n0, lmu, nu);
      n0 += 1;
      n += 1;
      storeVal[n] = old_term;
    }

    if (n == maxIter) return(log_sum_exp(storeVal[:n]));

    storeVal[n+1] = old_term - log_diff_exp(0, new_term - old_term) - log(2);

    return log_sum_exp(storeVal[:(n+1)]);
  }
}
