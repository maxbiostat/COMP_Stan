/* Functions for comparing approximate results to the true answer */
real  signum(real x) {
  real ans;
  if(x < 0){
    ans = -1;
  }else{
    if(x == 0){
      ans = 0;
    }else{
      ans = 1;
    }
  }
  return ans;
}
real robust_difference(real x, real y){
  real sgn = signum(x-y);
  real m = min({x, y});
  real M = max({x, y});
  return(sgn * exp(log_diff_exp(M, m)));
  //return((log_diff_exp(M, m)));
}

