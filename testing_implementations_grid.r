library(cmdstanr)
library(rstan)
################
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())
COMP_lpdf <- function(k, theta){
  lambda <- theta[1]
  nu <- theta[2]
  return(
    k * log(lambda) - nu*lfactorial(k)  
  )
}
#################
fpath <- "stan/compare_COMP_implementations.stan"
implementations <- cmdstanr::cmdstan_model(fpath, include_paths = "./stan/")

compare_implementations <- function(Mu, Nu, Eps, M){
  Theta <- c(Mu, Nu)
  if(Nu == 1){
    TrueValue <- Mu  
  }else{
    if(Nu == 2){
      TrueValue <- log(besselI(2*sqrt(Mu), nu = 0))
    }else{
      lps <- COMP_lpdf(k = 0:(5*M), theta = Theta)
      TrueValue <- matrixStats::logSumExp(lps)
    }
  }
  test.data <- list(
    log_mu = log(Mu),
    nu = Nu,
    eps = Eps,
    M = M,
    true_value = TrueValue
  )
  raw <- implementations$sample(data = test.data, chains = 1, 
                                iter_warmup = 0, iter_sampling = 1,
                                fixed_param = TRUE, show_messages = FALSE)
  results <- stanfit(raw)
  out <- extract(results)
  
  res.asymp <- data.frame(logZ = out$lZ_asymp,
                          niter = NA,
                          error = out$diff_asymp,
                          relative_error = out$rel_diff_asymp,
                          below_tolerance = out$getItRight_asymp,
                          method = "Gaunt")
  
  res.threshold <- data.frame(logZ = out$lZ_threshold[1],
                          niter = out$lZ_threshold[2],
                          error = out$diff_threshold,
                          relative_error = out$rel_diff_threshold,
                          below_tolerance = out$getItRight_threshold,
                          method = "ST")
  
  res.errorBounding <- data.frame(logZ = out$lZ_errorBounding[1],
                             niter = out$lZ_errorBounding[2],
                             error = out$diff_errorBounding,
                             relative_error = out$rel_diff_errorBounding,
                             below_tolerance = out$getItRight_errorBounding,
                             method = "errorBounding")
  
  ans <- do.call(rbind, list(res.asymp,
                             res.threshold,
                             res.errorBounding))
  ### 
  ans$mu <- Mu
  ans$nu <- Nu
  ans$max_iter <- M
  ans$eps <- Eps
  ans$true_value <- TrueValue
  return(ans)
}

epsilons <- 10^(-(1:20))
mus <- c(0.1, 0.2, 0.5, 1, 1.5, 2, 5, 10, 15, 50, 100, 500)
nus <- c(0.1, 0.2, 0.5, 1, 1.5, 2, 5)

grid <- expand.grid(epsilon = epsilons, mu = mus, nu = nus)

comp.time <- system.time(
  all.res <- lapply(1:nrow(grid), function(i){
    compare_implementations(Mu = grid[i,]$mu,
                            Nu = grid[i,]$nu,
                            Eps = grid[i,]$epsilon, M = 1E5)
  })
)
comp.time

all.res.dt <- do.call(rbind, all.res)
write.csv(all.res.dt,
          file = "data/COMP_comparisons.csv",
          row.names = FALSE)