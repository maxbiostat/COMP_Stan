/* Modified brms implementation that sums all terms at once instead of doing it in batches */
// Taken from https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_com_poisson.stan
real[] log_Z_COMP_brms_bulk(real log_mu, real nu, real eps, int M) {
  real log_Z;
  int k = 2;
  real leps = log(eps);
  int converged = 0;
  vector[M] log_Z_terms;
  int i = 1;
  if (nu == 1) {
    return {exp(log_mu), 0};
  }
  // nu == 0 or Inf will fail in this parameterization
  if (nu <= 0) {
    reject("nu must be positive");
  }
  if (nu == positive_infinity()) {
    reject("nu must be finite");
  }
  // first 2 terms of the series
  log_Z_terms[1] = log1p_exp(log_mu);
  while (converged == 0) {
    if(k >= M) break;
    // adding terms in batches simplifies the AD tape
    log_Z_terms[i + 1] = k * log_mu - nu*lgamma(k + 1);
    k += 1;
    if (log_Z_terms[i + 1] <= leps) {
      converged = 1;
      break;
    }
    i += 1;
  }
  log_Z = log_sum_exp(log_Z_terms[1:i]);
  return {log_Z, k};
}


