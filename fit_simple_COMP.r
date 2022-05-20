library(COMPoissonReg)
library(cmdstanr)
library(rstan)
####
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())
compress_counts <- function(x){
  tt <- table(x)
  return(
    list(
      n = as.numeric(tt),
      k = as.numeric(names(tt)),
      K = length(tt)
    )
  )
}
####
Mu <- 10
Nu <- .4
epsilon <- 1E-16
MaxIter <- 1E4

##
nobs <- 1000
# set.seed(666)
Y <- COMPoissonReg::rcmp(n = nobs, lambda = Mu, nu = Nu)
hist(Y, probability = TRUE)

cY <- compress_counts(Y)

stan.data <- list(
  K = cY$K,
  N = nobs,
  n = cY$n,
  y = cY$k,
  s_mu = .01,
  r_mu = .01,
  nu_sd = 1,
  eps = epsilon,
  M = MaxIter
)
iterations <- 500
############

## EB
fpath.eb <- "stan/MCMC_COMP_ErrorBounding.stan"
EB_impl <- cmdstanr::cmdstan_model(fpath.eb,
                                   include_paths = "./stan/")

opt_EB <- EB_impl$optimize(data = stan.data)
opt_EB$mle()

EB.raw <-
  EB_impl$sample(data = stan.data,
                 refresh = floor(iterations/5),
                 chains = 4, parallel_chains = 4,
                 iter_warmup = iterations, iter_sampling = iterations,
                 adapt_delta = .90, max_treedepth = 12,
                 show_messages = FALSE)
EB.mcmc <- stanfit(EB.raw)

## brms 
fpath.thresh <- "stan/MCMC_COMP_Threshold.stan"
threshold_impl <- cmdstanr::cmdstan_model(fpath.thresh, include_paths = "./stan/")

opt_threshold <- threshold_impl$optimize(data = stan.data)
opt_threshold$mle()

threshold.raw <-
  threshold_impl$sample(data = stan.data,
                   refresh = floor(iterations/5),
                   chains = 4, parallel_chains = 4,
                   iter_warmup = iterations, iter_sampling = iterations,
                   adapt_delta = .90, max_treedepth = 12,
                   show_messages = FALSE)
threshold.mcmc <- stanfit(threshold.raw)


#### 

EB.mcmc
threshold.mcmc

EB.raw$time()
threshold.raw$time()

check_hmc_diagnostics(EB.mcmc)
check_hmc_diagnostics(threshold.mcmc)


#### 
all.mcmcs <- list(
  EB = EB.mcmc,
  ST = threshold.mcmc
)

nits.df <- reshape2::melt(
  lapply(all.mcmcs, function(x) extract(x, 'n_iter')$n_iter)
)

library(ggplot2)

ggplot(nits.df, aes(y = value, colour = L1, fill = L1)) +
  geom_boxplot(alpha = 0.4) +
  scale_y_continuous("Number of iterations to approximate logZ") +
  theme_bw(base_size = 16)
