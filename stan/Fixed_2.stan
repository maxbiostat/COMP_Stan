/* Modified brms implementation that sums all terms at once instead of doing it in batches */
// Taken from https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_com_poisson.stan
real[] log_Z_COMP_2_fixed(real log_mu, real nu, int Nmax, int M) {
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
  return(
    infiniteSumToCap({log_mu, nu}, Nmax, M, 0)
  );
}
