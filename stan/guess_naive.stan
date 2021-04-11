/* "Naive" algorithm guessing the number of iterations needed using the algebra solver */
vector miniRootF(vector x, vector theta, real[] xr, int[] xi){
  vector[1] out;
  out[1] = x[1] * theta[1] - theta[2] * lgamma(x[1] + 1) - theta[3];
  return out;
}

real[] log_Z_COMP_GuessNaive(real logMu, real Nu, real epsilon,
      int M, data real[] xr, data int[] xi){
  real leps = log(epsilon);
  vector[1] guess = [exp(logMu / Nu)]';
  vector[1] r;
  int n = 0;
  real N;
  vector[M] values;
  vector[3] pars = [logMu, Nu, leps]';
  r = algebra_solver(miniRootF, guess, pars, xr, xi);
  N = min({r[1], M-1});
  while (n < N) {
    values[n + 1] = log_COM_Poisson(n, logMu, Nu);
    n += 1;
  }
  values[n + 1] = log_COM_Poisson(n, logMu, Nu);
  return {log_sum_exp(values[:(n+1)]), n + 1};
}

