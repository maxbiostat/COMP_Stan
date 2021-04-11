library(COMPoissonReg)
library(cmdstanr)
library(rstan)
#### Stan stuff
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())
get_diagnostics <- function (object) {
  require(rstan)
  out <- data.frame(
    n_divergences = get_num_divergent(object),
    n_maxtreedepth = get_num_max_treedepth(object),
    n_low_bfmi = sum(as.numeric((get_bfmi(object) < 0.2)))
  )
  return(out)
}
parse_performance <- function(object){
  moni <- monitor(object, print = FALSE)
  summy <- as.data.frame(moni)[, c("mean", "2.5%",
                                   "97.5%", "Rhat",
                                   "n_eff", "Bulk_ESS",
                                   "Tail_ESS")]
  summy$parameter <- rownames(summy)
  return( reshape2::melt(summy) )
}
cv <- function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
get_times <- function(obj){
  times <- obj$time()
  out <- data.frame(
    total = times$total,
    warmup = sum(times$chains$warmup),
    sampling = sum(times$chains$sampling),
    cv_warmup = cv(times$chains$warmup),
    cv_sampling = cv(times$chains$sampling)
  )
}
### Data gen and parsing stuff
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
############
adaptive_impl <- cmdstanr::cmdstan_model("stan/MCMC_COMP_adaptive.stan",
                                         include_paths = "./stan/")
brms_impl <- cmdstanr::cmdstan_model("stan/MCMC_COMP_brms.stan",
                                     include_paths = "./stan/")

simulate_fit_report <- function(rep, Mu, Nu, nobs = 1000, 
                                epsilon = 1E-16,
                                MaxIter = 1E4, iterations = 500){
  ## Prep
  simu.info <- data.frame(
    mu = Mu,
    nu = Nu,
    nobs = nobs,
    epsilon = epsilon,
    max_iter =  MaxIter,
    n_iter = iterations,
    replicate = rep
  )
  
  ## Simulate
  cY <- compress_counts(COMPoissonReg::rcmp(n = nobs, lambda = Mu, nu = Nu))
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
  ## Fit
  b <- capture.output(
    brms.raw <-
      brms_impl$sample(data = stan.data, refresh = 0, chains = 4,
                       parallel_chains = 4, iter_warmup = iterations,
                       adapt_delta = .90,  max_treedepth = 10,
                       iter_sampling = iterations, show_messages = FALSE)
  )
  brms.mcmc <- stanfit(brms.raw)
  
  a <- capture.output(
    adaptive.raw <-
      adaptive_impl$sample(data = stan.data, refresh = 0, chains = 4,
                           parallel_chains = 4, iter_warmup = iterations,
                           adapt_delta = .90,  max_treedepth = 10,
                           iter_sampling = iterations, show_messages = FALSE) 
  )
  adaptive.mcmc <- stanfit(adaptive.raw)
  ## Report
  diag.dt <- rbind(data.frame(get_diagnostics(brms.mcmc),
                              implementation = "brms"),
                   data.frame(get_diagnostics(adaptive.mcmc),
                              implementation = "adaptive"))
  diag.dt <- data.frame(diag.dt, simu.info)
  #
  performance.dt <- rbind(data.frame(parse_performance(brms.mcmc),
                                     implementation = "brms"),
                          data.frame(parse_performance(adaptive.mcmc),
                                     implementation = "adaptive"))
  performance.dt <- data.frame(performance.dt, simu.info)
  #
  time.dt <- rbind(data.frame(get_times(brms.raw),
                              implementation = "brms"),
                   data.frame(get_times(adaptive.raw),
                              implementation = "adaptive"))
  time.dt <- data.frame(time.dt, simu.info)
  #########
  result <- list(
    data = cY,
    diagnostics = diag.dt,
    ess_stuff = performance.dt,
    timings = time.dt
  )
  return(result)
}
do_batch <- function(n_reps, Mu, Nu, nobs = 1000, 
                     epsilon = 1E-16, MaxIter = 1E4,
                     iterations = 500, verbose = TRUE){
  all.simus <- lapply(1:n_reps,
                      function(i){
                        if(verbose) cat("Doing replicate:", i, "\n")
                        simulate_fit_report(rep = i,
                                            Mu = Mu, Nu = Nu, nobs = nobs,
                                            epsilon = epsilon,
                                            MaxIter = MaxIter,
                                            iterations = iterations)  
                      }
  )
  return(all.simus)
}
parse_simu <- function(simu){
  list(
    diagnostics = do.call(rbind, lapply(simu, function(x) x$diagnostics )),
    performance = do.call(rbind, lapply(simu, function(x) x$ess_stuff )),
    timings = do.call(rbind, lapply(simu, function(x) x$timings))
  )
}

true.Mu <- 10
true.Nu <- .5
simu <-  do_batch(n_reps = 10,
                  Mu = true.Mu, Nu = true.Nu, nobs = 1000,
                    epsilon = 1E-16,
                  MaxIter = 1E4, iterations = 500)

save(simu, file = paste0("data/Mu_", true.Mu, "_Nu_", true.Nu, ".RData"))

parsed <- parse_simu(simu)

parsed$diagnostics

ggplot(parsed$timings, aes(x = implementation, y = total)) + geom_boxplot()

ggplot(subset(parsed$performance, variable == "mean"),
       aes(x = implementation, y = value,
           colour = implementation, fill = implementation)) +
  geom_boxplot(alpha = 0.4) +
  scale_x_discrete("") +
  scale_y_continuous("") +
  facet_wrap(.~parameter, scales = "free_y") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") 

ESS <- subset(parsed$performance, variable == "n_eff")
ESS <- merge(ESS, parsed$timings, by = c("implementation", "replicate"))

ggplot(ESS, 
       aes(x = implementation, y = value/(warmup + sampling),
           fill = implementation) ) +
  geom_boxplot(alpha = .4) +
  theme_bw(base_size = 16) + 
  scale_y_continuous("ESS per second")