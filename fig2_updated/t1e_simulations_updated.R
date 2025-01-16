library(radEmu)
library(TreeSummarizedExperiment)
library(ANCOMBC)
library(ALDEx2)
library(DESeq2)
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
                        "aldex_est" = NA,
                        "aldex_p" = NA,
                        "ancom_est" = NA,
                        "ancom_p" = NA,
                        "clr_est" = NA,
                        "clr_p" = NA,
                        "deseq_est" = NA,
                        "deseq_p" = NA,
                        "rademu_est" = NA,
                        "wald_p" = NA,
                        "score_p" = NA)
  
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
    
    # ALDEx2 
    clr_data <- aldex.clr(t(Y), X) 
    aldex_model <- try({
      suppressMessages(
        aldex.glm(clr = clr_data, X)
      )
    })
    if (!inherits(aldex_model, "try-error")) {
      rezzies[sim, "aldex_est"] <- aldex_model$`cov:Est`[J/2]
      rezzies[sim, "aldex_p"] <- aldex_model$`cov:pval`[J/2]
    }
    
    # ANCOM-BC2
    rownames(Y) <- paste0("sample", 1:nrow(Y))
    colnames(Y) <- paste0("category", 1:ncol(Y))
    tse_data <- TreeSummarizedExperiment(assays = list("counts" = t(Y)),
                                         colData = data.frame(cov = X[, 2]))
    ancom_model <- try({
      suppressMessages(
        ancombc2(data = tse_data,
                 assay_name = "counts", 
                 p_adj_method = "none",
                 fix_formula = "cov",
                 prv_cut = 0,
                 alpha = 0.05,
                 pseudo_sens = FALSE)
      )
    })
    if (!inherits(ancom_model, "try-error")) {
      rezzies[sim, "ancom_est"] <- ancom_model$res$lfc_cov[J/2]
      rezzies[sim, "ancom_p"] <- ancom_model$res$p_cov[J/2]
    }
    
    # CLR t-test with pseudocount of 1 for all counts
    Y_pseudo <- Y + 1
    log_Y <- log(Y_pseudo)
    mean_log_Y <- rowMeans(log_Y)
    clr_data <- log_Y - matrix(mean_log_Y, nrow = n, ncol = J, byrow = FALSE)
    est <- sapply(1:J, function(x) {mean(clr_data[X[, 2] == 1, x]) - 
        mean(clr_data[X[, 2] == 0, x])})
    rezzies[sim, "clr_est"] <- est[J/2]
    ttest_res <- try({t.test(clr_data[X[, 2] == 1, J/2], 
                        clr_data[X[, 2] == 0, J/2])})
    if (!inherits(ttest_res, "try-error")) {
      rezzies[sim, "clr_p"] <- ttest_res$p.value
    }
    
    # DESeq2
    dds <- DESeqDataSetFromMatrix(countData = t(Y_pseudo),
                                  colData = data.frame(X = as.factor(X[, 2])),
                                  design = ~ X)
    dds <- try({DESeq(dds)})
    if (!inherits(dds, "try-error")) {
      deseq_res <- results(dds, name = "X_1_vs_0")
      rezzies[sim, "deseq_est"] <- deseq_res$log2FoldChange[J/2]
      rezzies[sim, "deseq_p"] <- deseq_res$pvalue[J/2]
    }
    
    # radEmu 
    emu_model <- try({
      emuFit(Y = Y,
             X = X,
             test_kj = data.frame(k = 2, j = J/2),
             tau = 2,
             B_null_tol = 0.005,
             tolerance = 0.005,
             constraint_tol = 0.001,
             return_wald_p = TRUE,
             unobserved_taxon_error = FALSE)})
    if (!inherits(emu_model, "try-error")) {
      rezzies[sim, "rademu_est"] <- emu_model$coef$estimate[J/2]
      rezzies[sim, "wald_p"] <- emu_model$coef$wald_p[J/2]
      rezzies[sim, "score_p"] <- emu_model$coef$pval[J/2]
    }
  }
  
  return(rezzies)
}

results <- sim_function(setting, seeds) 
saveRDS(results, paste0(sim_dir_name, "/setting", setting, "_first_seed", seeds[1], ".rds"))