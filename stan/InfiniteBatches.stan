// c-folding inifinite sum algorithm
// Requires definition of logFunction with two arguments:
// int k and real[] parameters
array[] real infiniteBatches(array[] real p, int batch_size, real epsilon, int maxIter, int n0) {
  array[batch_size + 1] real storeVal;
  real leps = log(epsilon);
  int n0_ = n0;
  int nextCheckPoint = batch_size - 1;
  int lastCheckPoint = nextCheckPoint + 1;
  real summedIncrement = negative_infinity();
  int n = batch_size;
  real check;
  storeVal[1] = negative_infinity();
  // Setting up first batches
  for (i in 2:(batch_size + 1)) {
    storeVal[i] = logFunction(n0_, p);
    n0_ += 1;
  }
  storeVal[1] = log_sum_exp(sort_asc(storeVal));

  // Find the maximum
  while (n < maxIter && (
    n <= batch_size ||
    summedIncrement > leps ||
    storeVal[batch_size + 1] - storeVal[batch_size] >
      -log1p_exp(storeVal[batch_size + 1] - summedIncrement) ||
    storeVal[batch_size + 1] > storeVal[batch_size] ||
    is_inf(storeVal[batch_size])
  )) {
    storeVal[1] = log_sum_exp(storeVal[1], summedIncrement);
    for (i in 2:(batch_size + 1)) {
      storeVal[i] = logFunction(n0_, p);
      n0_ += 1;
    }
    summedIncrement = log_sum_exp(sort_asc(storeVal[2:]));
    n += batch_size;
    if (n >= maxIter) break; // Return if maxIter is reached
  }
  return {log_sum_exp(storeVal[1], summedIncrement), 1. * n};
}
