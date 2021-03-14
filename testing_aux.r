log_sum_exp <- function(x){
  # log(sum(exp(x - max(x)))) + max(x) ## unstable
  ans <- matrixStats::logSumExp(lx = x)
  return(ans)
}
##
log1m <- function(x) {
  return(log1p(-x))
}
##
log1m_exp <- function(a) {
  if (a > -0.693147) {
    return(log(-expm1(a))) 
  } else {
    return(log1m(exp(a)));
  }
}
##
log_diff_exp <- function (x, y){
  return(x + log1m_exp(y - x))
}  
##
robust_difference <- function(x, y, log = FALSE){
  sgn <- sign(x-y)
  if(log){
    ans <- log_diff_exp(max(x, y), min(x, y))
  }else{
    ans <- sgn * exp(log_diff_exp(max(x, y), min(x, y)))
  }
  return(ans)
}
##
relative_difference <- function(x, y){
  obj <- all.equal(as.numeric(x), as.numeric(y), tolerance = 0)
  ans <- suppressWarnings(as.numeric(gsub("Mean relative difference: ", "", obj)))
  suppressWarnings(if(is.na(ans)) ans <- 0 )
  return(ans)
}
