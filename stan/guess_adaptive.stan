/* Adaptive algorithm guessing the number of iterations needed using the algebra solver */
vector rootF(vector x, vector theta, real[] xr, int[] xi){
  vector[1] out;
  out[1] = - theta[2] * lgamma(x[1] + 1) + (x[1] + 1) * theta[1] - log_diff_exp(theta[2] * log(x[1] + 1), theta[1]) - theta[3];
  return out;
}
real[] log_Z_COMP_GuessAdaptive(real logMu, real Nu, real epsilon,
                    int M, data real[] xr, data int[] xi){
  real leps = log(epsilon) + log(2);
  real anp1;
  vector[1] guess = [exp(logMu / Nu)]';
  vector[1] r;
  int n = 0;
  real N;
  vector[M] values;
  vector[3] pars = [logMu, Nu, leps]';
  r = algebra_solver(rootF, guess, pars, xr, xi);
  N = min({r[1], M-3});
  while (n < N) {
    values[n + 1] = log_COM_Poisson(n, logMu, Nu);
    n += 1;
  }
  values[n + 1] = log_COM_Poisson(n, logMu, Nu);
  anp1 = log_COM_Poisson(n + 1, logMu, Nu);
  values[n + 2] = anp1 - log_diff_exp(0, anp1 - values[n + 1]) - log(2);
  values[n + 3] = anp1 - log(2);
  return {log_sum_exp(values[:(n + 3)]), n + 1};
}

