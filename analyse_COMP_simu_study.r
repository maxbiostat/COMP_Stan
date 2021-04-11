parse_simu <- function(simu){
  list(
    diagnostics = do.call(rbind, lapply(simu, function(x) x$diagnostics )),
    performance = do.call(rbind, lapply(simu, function(x) x$ess_stuff )),
    timings = do.call(rbind, lapply(simu, function(x) x$timings))
  )
}
get_ess <- function(sim){
  parsed <- parse_simu(sim)
  ESS <- subset(parsed$performance, variable == "n_eff")
  ESS <- merge(ESS, parsed$timings, all.x = TRUE)
  return(ESS)
}
read_ess_rdata <- function(fname){
  load(fname)
  ss <- strsplit(fname, "_")
  info <- data.frame(
    mu = as.numeric(ss[[1]][2]),
    nu = as.numeric(gsub(".RData", "", ss[[1]][4]))
  )
  parsed <- parse_simu(simu)
  out <- subset(parsed$performance, variable == "n_eff")
  out <- merge(out, parsed$timings, all.x = TRUE)
  out$mu <- info$mu
  out$nu <- info$nu
  return(out)
}
##########
rdatas <- system("ls data/*.RData", intern = TRUE)
# fname <- rdatas[1]
all.res <- do.call(rbind, lapply(rdatas, read_ess_rdata))
head(all.res)
all.res$parameterComb <- paste0("mu_", all.res$mu, "_nu_", all.res$nu)
library(ggplot2)

ggplot(all.res, 
       aes(x = implementation, y = value/(warmup + sampling),
           fill = implementation) ) +
  geom_boxplot(alpha = 0.4) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")+
  facet_wrap(parameterComb~., scales = "free_y") +
  scale_y_continuous("ESS per second")
