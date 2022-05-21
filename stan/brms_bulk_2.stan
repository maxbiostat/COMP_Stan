/* Modified brms implementation that sums all terms at once instead of doing it in batches */
// Taken from https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_com_poisson.stan
real[] log_Z_COMP_brms_bulk2(real log_mu, real nu, real eps, int M) {
  int k = 2;
  real leps = log(eps);
  vector[M] log_Z_terms;
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
  log_Z_terms[1] = 0;
  log_Z_terms[2] = nu*log_mu;
  while (log_Z_terms[k] >= leps && k < M) {
    k += 1;
    log_Z_terms[k] = nu*((k - 1) * log_mu - lgamma(k));
  }
  return {log_sum_exp(log_Z_terms[:k]), k};
}
