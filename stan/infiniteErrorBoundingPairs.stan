// Internal use
int adaptiveConvergenceCheck(real oldT, real newT, real lepsilon, real log1mL) {
  // if L = 0, there is a simpler check
  if (!log1mL) return oldT - newT < log1p_exp(newT - lepsilon);
  
  real logZ = newT + oldT - log_diff_exp(oldT, newT);
  real ls = newT - log1mL;
  
  if (logZ > ls){
    return log_diff_exp(logZ, ls) >= lepsilon;
  }else{
    return log_diff_exp(ls, logZ) >= lepsilon;
  }
}

// Adaptive inifinite sum algorithm
// Requires definition of logFunction with two arguments:
// int k and real[] parameters
array[] real infiniteErrorBoundingPairs(array[] real p, real epsilon, int maxIter, real logL, int n0) {
  vector[maxIter + 1] storeVal;
  real leps = log(epsilon) + log2();
  int n = 2;
  int n0_ = n0;
  real temp;
  real log1mL = log_diff_exp(0, logL);
  
  // Setting up first iterations
  storeVal[1] = logFunction(n0_, p);
  n0_ += 1;
  storeVal[2] = logFunction(n0_, p);
  
  // Find the maximum
  while (storeVal[n] > storeVal[n - 1]) {
    n0_ += 1;
    n += 1;
    storeVal[n] = logFunction(n0_, p);
    if (n >= maxIter) return({log_sum_exp(storeVal[:n]), 1. * n}); // Return if maxIter is reached
  }
  // Start testing convergence after the maximum
  while ( adaptiveConvergenceCheck(storeVal[n - 1], storeVal[n], leps, log1mL) )  {
    n0_ += 1;
    n += 1;
    storeVal[n] = logFunction(n0_, p);
    if (n >= maxIter) break;
  }
  temp = storeVal[n];
  storeVal[n] = temp - log_diff_exp(0, temp - storeVal[n - 1]) - log2();
  storeVal[n + 1] = temp - log1mL - log2();
  return {log_sum_exp(sort_asc(storeVal[:(n + 1)])), 1. * n};
}
