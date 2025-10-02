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
sim_dir_name <- "radEmu/results/t1e_sims/"

#parameters determining form of simulations
sample_sizes <- c(250, 50, 10)
Js <- c(250, 50, 10)
distns <- c("ZINB", "Poisson")
sim_settings <- expand.grid(sample_sizes, Js, distns)
colnames(sim_settings) <- c("n", "J", "distn")
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
  distn <- sim_settings[sim_no,"distn"]
  b0 <- B_list[[which_J]]$b0
  b1 <- B_list[[which_J]]$b1
  X <- cbind(intercept = 1, cov = rep(0:1, each = n/2))
  rezzies <- data.frame("n" = n,
                        "J" = J,
                        "distn" = distn,
                        "seed" = seeds,
                        "linda_est" = NA,
                        "linda_p" = NA,
                        "locom_est" = NA,
                        "locom_p" = NA,
                        "maaslin_est" = NA,
                        "maaslin_p" = NA)
  
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
                       feature.dat.type = "count")})
    if (!inherits(linda_res, "try-error")) {
      rezzies[sim, "linda_est"] <- linda_res$output$cov$log2FoldChange[J/2]
      rezzies[sim, "linda_p"] <- linda_res$output$cov$pvalue[J/2]
    }
    
    # locom
    colnames(Y) <- paste0("taxon", 1:ncol(Y))
    locom_res <- try({locom(otu.table = Y, Y = X[, 2], filter.thresh = 0)})
    if (!inherits(locom_res, "try-error")) {
      ind <- which(colnames(locom_res$effect.size) == paste0("taxon", J/2))
      
      if (length(ind) == 1) {
        rezzies[sim, "locom_est"] <- locom_res$effect.size[ind]
        rezzies[sim, "locom_p"] <- locom_res$p.otu[ind]
      }
    }
    
    # maaslin
    dat <- data.frame(cov = X[, 2])
    rownames(Y) <- paste0("sample", 1:nrow(Y))
    rownames(dat) <- rownames(Y)
    file <- paste0("~/radEmu/results/t1e_sims/maaslin_output/maaslin_output_setting", setting, "_seed", seeds[sim])
    #file <- paste0("maaslin_output_setting", setting, "_seed", seeds[sim])
    maaslin_res <- try({maaslin3(input_data = Y, input_metadata = dat,
                            formula = ~cov, output = file, 
                            plot_summary_plot = FALSE, plot_associations = FALSE,
                            save_models = FALSE)})
    if (!inherits(maaslin_res, "try-error") & file.exists(paste0(file, "/all_results.tsv"))) {
      maas_res <- read.delim(paste0(file, "/all_results.tsv")) %>%
        filter(feature == paste0("taxon", J/2), model == "abundance")
      if (nrow(maas_res) > 0) {
        rezzies[sim, "maaslin_est"] <- maas_res$coef
        rezzies[sim, "maaslin_p"] <- maas_res$pval_individual
      }
    }
  }
  
  return(rezzies)
}

results <- sim_function(setting, seeds) 
saveRDS(results, paste0(sim_dir_name, "/additional_setting", setting, "_first_seed", seeds[1], ".rds"))