library(ggplot2)
library(dplyr)
library(glue)

gr <- (sqrt(5) + 1) / 2
stdW <- 8

comp <- read.csv("data/COMP_comparisons.csv")
head(comp)
names(comp)


## Methods: Gaunt (NA values), Naive, Old_brms, New_brms, Adaptive

for (mtd in c("ST", "errorBounding")){
  g <- ggplot(na.omit(comp) %>% filter(method == mtd),
              aes(y = as.factor(mu), x = as.factor(nu))) + theme_bw()+
    geom_tile(aes(fill = as.factor(below_tolerance),
                  color = as.factor(below_tolerance))) +
    facet_wrap(vars(as.factor(log10(eps)))) +
    labs(x = "Nu", y = "Mu", fill = "Beat tol.") + guides(color = FALSE) +
    theme(axis.text = element_text(size = 5))
  ggsave(glue("figs/{mtd}_beatTol.pdf"), plot = g, device="pdf",
         width = stdW, height = stdW / gr)
}

## percent comparisons
comp %>%
  group_by(method) %>%
  summarise(perc_beat_tol = mean(below_tolerance)) %>%
  ungroup()

comp %>%
  group_by(eps) %>%
  summarise(gaunt = mean(below_tolerance[method == "Gaunt"]),
            naive = mean(below_tolerance[method == "ST"]),
            adaptive = mean(below_tolerance[method == "errorBounding"])) %>%
  ungroup()