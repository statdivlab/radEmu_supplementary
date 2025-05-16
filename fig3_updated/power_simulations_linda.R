library(radEmu)
library(maaslin3)
library(MicrobiomeStat)
library(LOCOM)
library(tidyverse)

# run script to get function "generate_test_data" within folder "radEmu/scripts"
source("radEmu/scripts/generate_test_data.R")

# get seed
args <- commandArgs(trailingOnly = FALSE)
if (length(args) == 0) {
  batch <- 1
} else {
  arg <- args[length(args)]
  batch <- abs(readr::parse_number(arg))
}

# parallelize 500 trials for 18 settings across 900 tasks
nsim <- 500
sim_per_script <- 10
script_per_setting <- 50
seeds <- ((batch - 1) %% script_per_setting) * sim_per_script + 1:sim_per_script
setting <- floor((batch - 1) / script_per_setting) + 1

# set directory to save results
sim_dir_name <- "radEmu/results/power_sims/"

#parameters determining form of simulations
sample_sizes <- c(250, 50)
Js <- c(10)
distns <- c("ZINB")
alt <- c(0, 0.5, 1:5)
sim_settings <- expand.grid(sample_sizes, Js, distns, alt)
colnames(sim_settings) <- c("n", "J", "distn", "alt")
sim_settings <- as.data.frame(sim_settings)

# generate B matrices
set.seed(9443)
B_list <- lapply(Js,
                 function(j)
                   radEmu:::get_sim_bs(j))

# function to run simulations
sim_function <- function(sim_no, seeds) {
  
  n <- sim_settings[sim_no,"n"]
  J <- sim_settings[sim_no,"J"]
  which_J <- which(Js ==J)
  distn <- sim_settings[sim_no, "distn"]
  alt <- sim_settings[sim_no, "alt"]
  b0 <- B_list[[which_J]]$b0
  b1 <- B_list[[which_J]]$b1
  b1[J/2] <- alt
  X <- cbind(intercept = 1, cov = rep(0:1, each = n/2))
  rezzies <- data.frame("n" = n,
                        "J" = J,
                        "distn" = distn,
                        "alt" = alt,
                        "seed" = seeds,
                        "linda_est" = NA,
                        "linda_p" = NA)
  
  for (sim in 1:length(seeds)) {
    print(sim)
    set.seed(seeds[sim])
    Y <- generate_test_data(n = n,
                            J = J,
                            b0 = b0,
                            b1 = b1,
                            distn = distn,
                            zinb_size = 5,
                            zinb_zero_prop = 0.6,
                            mean_count_before_ZI = 50)
    
    # linda
    features <- t(Y)
    meta <- data.frame(cov = X[, 2])
    linda_res <- try({linda(feature.dat = features, meta.dat = meta, formula = '~cov',
                            feature.dat.type = "count", is.winsor = FALSE)})
    if (!inherits(linda_res, "try-error")) {
      rezzies[sim, "linda_est"] <- linda_res$output$cov$log2FoldChange[J/2]
      rezzies[sim, "linda_p"] <- linda_res$output$cov$pvalue[J/2]
    }
  }
  
  return(rezzies)
}

results <- sim_function(setting, seeds) 
saveRDS(results, paste0(sim_dir_name, "/linda_additional_setting", setting, "_first_seed", seeds[1], ".rds"))