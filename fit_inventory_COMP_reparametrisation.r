library(COMPoissonReg)
library(cmdstanr)
library(rstan)
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())
####
epsilon <-  1E-16 # .Machine$double.eps
MaxIter <- 1E4
Ncutoff <- 100

inventory <- read.csv("data/Shmuelli_2005.csv")

nobs <- sum(inventory$frequency)
cY <- list(
  k = inventory$count,
  n = inventory$frequency,
  K = nrow(inventory)
)

stan.data <- list(
  K = cY$K,
  N = nobs,
  n = cY$n,
  y = cY$k,
  s_mu = .01,
  r_mu = .01,
  s_nu = 0.0625,
  r_nu = 0.25,
  eps = epsilon,
  M = MaxIter,
  batchSize = 10,
  N_cap = Ncutoff
)
iterations <- 5000
############

## EB
fpath.eb <- "stan/MCMC_COMP_2_ErrorBounding.stan"
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

## Threshold 
fpath.thresh <- "stan/MCMC_COMP_2_Threshold.stan"
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

## Fixed Cap 
fpath.fixed <- "stan/MCMC_COMP_2_fixedCap.stan"
fixed_impl <- cmdstanr::cmdstan_model(fpath.fixed, include_paths = "./stan/")

opt_fixed <- fixed_impl$optimize(data = stan.data)
opt_fixed$mle()

fixed.raw <-
  fixed_impl$sample(data = stan.data,
                    refresh = floor(iterations/5),
                    chains = 4, parallel_chains = 4,
                    iter_warmup = iterations, iter_sampling = iterations,
                    adapt_delta = .90, max_treedepth = 12,
                    show_messages = FALSE)

## Batches 
fpath.batch <- "stan/MCMC_COMP_2_Batches.stan"
batches_impl <- cmdstanr::cmdstan_model(fpath.batch, include_paths = "./stan/")

opt_batches <- batches_impl$optimize(data = stan.data)
opt_batches$mle()

batches.raw <-
  batches_impl$sample(data = stan.data,
                      refresh = floor(iterations/5),
                      chains = 4, parallel_chains = 4,
                      iter_warmup = iterations, iter_sampling = iterations,
                      adapt_delta = .90, max_treedepth = 12,
                      show_messages = FALSE)


#### 

all.raws <- list(
  EB = EB.raw,
  ST = threshold.raw,
  BA = batches.raw,
  FC = fixed.raw
)

all.mcmcs <- lapply(all.raws, stanfit)

lapply(all.mcmcs, function(x){
  print(x, pars = c("mu", "nu", "n_iter"),
        digits_summary = 3)
})

lapply(all.mcmcs, check_hmc_diagnostics)

lapply(all.raws, function(x) x$time())

efficiencies <- lapply(seq_along(all.raws), function(j){
  ESS <- monitor(all.mcmcs[[j]], print = FALSE)$n_eff[c(1, 2)]
  ttime <- all.raws[[j]]$time()$total
  return(ESS/ttime)
})
names(efficiencies) <- names(all.mcmcs)

efficiencies

#### 
nits.df <- reshape2::melt(
  lapply(all.mcmcs,
         function(x) extract(x, 'n_iter')$n_iter)
)

library(ggplot2)

ggplot(nits.df,
       aes(y = value, colour = L1, fill = L1)) +
  geom_boxplot(alpha = 0.4) +
  scale_y_continuous("Iterations to approximate logZ") +
  theme_bw(base_size = 16)