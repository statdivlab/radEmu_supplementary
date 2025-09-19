library(radEmu)
library(TreeSummarizedExperiment)
library(ANCOMBC)
library(ALDEx2)
library(DESeq2)
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

# parallelize 504 trials for 126 settings across 882 tasks
nsim <- 504
sim_per_script <- 72
script_per_setting <- 7
seeds <- ((batch - 1) %% script_per_setting) * sim_per_script + 1:sim_per_script
setting <- floor((batch - 1) / script_per_setting) + 1

# set directory to save results
sim_dir_name <- "radEmu/results/sims/"

#parameters determining form of simulations
sample_sizes <- c(250, 50, 10)
Js <- c(250, 50, 10)
distns <- c("ZINB", "Poisson")
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
  b2 <- sample(b1, size = J, replace = FALSE)
  test_ind <- c(ceiling(J / 4), J / 2, ceiling(3 * J / 4))
  b1[test_ind[1]] <- -alt
  b1[test_ind[2]] <- alt
  b1[test_ind[3]] <- alt 
  X <- cbind(intercept = 1, 
             cov1 = runif(n, min = 0, max = 1),
             cov2 = rnorm(n, mean = 0, sd = 2))
  rezzies <- data.frame("n" = n,
                        "J" = J,
                        "distn" = distn,
                        "alt" = alt,
                        "seed" = seeds,
                        "aldex_est" = NA,
                        "aldex_p" = NA,
                        "ancom_est" = NA,
                        "ancom_p" = NA,
                        "clr_est" = NA,
                        "clr_p" = NA,
                        "deseq_est" = NA,
                        "deseq_p" = NA,
                        "linda_est" = NA,
                        "linda_p" = NA,
                        "locom_est" = NA,
                        "locom_p" = NA,
                        "maaslin_est" = NA,
                        "maaslin_p" = NA,
                        "rademu_est" = NA,
                        "wald_p" = NA,
                        "score_p" = NA)
  
  for (sim in 1:length(seeds)) {
    print(sim)
    set.seed(seeds[sim])
    Y <- generate_test_data(n = n,
                            J = J,
                            B = rbind(b0, b1, b2), 
                            X = X,
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
      rezzies[sim, "aldex_est"] <- aldex_model$`cov1:Est`[J/2]
      rezzies[sim, "aldex_p"] <- aldex_model$`cov1:pval`[J/2]
    }
    
    # ANCOM-BC2
    rownames(Y) <- paste0("sample", 1:nrow(Y))
    colnames(Y) <- paste0("category", 1:ncol(Y))
    tse_data <- TreeSummarizedExperiment(assays = list("counts" = t(Y)),
                                         colData = data.frame(cov1 = X[, 2],
                                                              cov2 = X[, 3]))
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
      rezzies[sim, "ancom_est"] <- ancom_model$res$lfc_cov1[J/2]
      rezzies[sim, "ancom_p"] <- ancom_model$res$p_cov1[J/2]
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
    
    # linda
    features <- t(Y)
    meta <- data.frame(cov = X[, 2])
    linda_res <- try({linda(feature.dat = features, meta.dat = meta, formula = '~cov',
                            feature.dat.type = "count", is.winsor = FALSE)})
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
    file <- paste0("~/radEmu/results/power_sims/maaslin_output/maaslin_output_setting", setting, "_seed", seeds[sim])
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