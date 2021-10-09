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
real robust_difference(real x, real y, int relative){
  real sgn = signum(x-y);
  real m = min({x, y});
  real M = max({x, y});
  real ans;
  if(relative==1){
    ans = sgn * exp(log_diff_exp(M, m)-x);
  }else{
    ans = sgn * exp(log_diff_exp(M, m));
  }
  return(ans);
}

