library(radEmu) # version 2.1.1.1
library(mia) # version 1.14.0
library(TreeSummarizedExperiment) # version 2.14.0
library(ANCOMBC) # version 2.8.1
library(ALDEx2) # version 1.38.0
library(DESeq2) # version 1.46.0
library(maaslin3) # version 0.99.16
library(MicrobiomeStat) # version 1.2
library(LOCOM) # version 2.0
library(tidyverse) # version 2.0.0

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

# remaining jobs that we need to run
remain_df <- data.frame(setting = rep(c(27, 45, 45, 99, 117, 117), each = 8), 
                        seed = c(441:448, 33:40, 209:216, 113:120, 57:64, 113:120))
setting <- remain_df$setting[batch]
seeds <- remain_df$seed[batch]

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

sim_no <- setting

n <- sim_settings[sim_no,"n"]
J <- sim_settings[sim_no,"J"]
which_J <- which(Js ==J)
distn <- sim_settings[sim_no, "distn"]

print(paste0("n = ", n, ", J = ", J, " distn = ", distn))

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
rezzies1 <- data.frame("n" = n,
                       "J" = J,
                       "distn" = distn,
                       "alt" = -alt,
                       "delta_level" = "low", 
                       "seed" = seeds,
                       "aldex_est" = NA,
                       "aldex_p" = NA,
                       "aldex_time" = NA,
                       "ancom_est" = NA,
                       "ancom_p" = NA,
                       "ancom_time" = NA,
                       "clr_est" = NA,
                       "clr_p" = NA,
                       "clr_time" = NA,
                       "deseq_est" = NA,
                       "deseq_p" = NA,
                       "deseq_time" = NA,
                       "linda_est" = NA,
                       "linda_p" = NA,
                       "linda_time" = NA,
                       "locom_est" = NA,
                       "locom_p" = NA,
                       "locom_time" = NA,
                       "maaslin_est" = NA,
                       "maaslin_p" = NA,
                       "maaslin_time" = NA,
                       "rademu_est" = NA,
                       "wald_p" = NA,
                       "score_p" = NA,
                       "rademu_time" = NA)
rezzies2 <- data.frame("n" = n,
                       "J" = J,
                       "distn" = distn,
                       "alt" = alt,
                       "delta_level" = "mid", 
                       "seed" = seeds,
                       "aldex_est" = NA,
                       "aldex_p" = NA,
                       "aldex_time" = NA,
                       "ancom_est" = NA,
                       "ancom_p" = NA,
                       "ancom_time" = NA,
                       "clr_est" = NA,
                       "clr_p" = NA,
                       "clr_time" = NA,
                       "deseq_est" = NA,
                       "deseq_p" = NA,
                       "deseq_time" = NA,
                       "linda_est" = NA,
                       "linda_p" = NA,
                       "linda_time" = NA,
                       "locom_est" = NA,
                       "locom_p" = NA,
                       "locom_time" = NA,
                       "maaslin_est" = NA,
                       "maaslin_p" = NA,
                       "maaslin_time" = NA,
                       "rademu_est" = NA,
                       "wald_p" = NA,
                       "score_p" = NA,
                       "rademu_time" = NA)
rezzies3 <- data.frame("n" = n,
                       "J" = J,
                       "distn" = distn,
                       "alt" = alt,
                       "delta_level" = "high", 
                       "seed" = seeds,
                       "aldex_est" = NA,
                       "aldex_p" = NA,
                       "aldex_time" = NA,
                       "ancom_est" = NA,
                       "ancom_p" = NA,
                       "ancom_time" = NA,
                       "clr_est" = NA,
                       "clr_p" = NA,
                       "clr_time" = NA,
                       "deseq_est" = NA,
                       "deseq_p" = NA,
                       "deseq_time" = NA,
                       "linda_est" = NA,
                       "linda_p" = NA,
                       "linda_time" = NA,
                       "locom_est" = NA,
                       "locom_p" = NA,
                       "locom_time" = NA,
                       "maaslin_est" = NA,
                       "maaslin_p" = NA,
                       "maaslin_time" = NA,
                       "rademu_est" = NA,
                       "wald_p" = NA,
                       "score_p" = NA,
                       "rademu_time" = NA)

sim <- 1
print(paste0("setting: ", sim_no, ", seed: ", seeds))
set.seed(seeds[sim])
Y <- generate_test_data_include_delta(n = n,
                                      J = J,
                                      B = rbind(b0, b1, b2), 
                                      X = X,
                                      distn = distn,
                                      zinb_size = 5,
                                      zinb_zero_prop = 0.6,
                                      mean_count_before_ZI = 50)$Y

# ALDEx2 
aldex_start <- proc.time()
clr_data <- aldex.clr(t(Y), X) 
aldex_model <- try({
  suppressMessages(
    aldex.glm(clr = clr_data, X)
  )
})
aldex_end <- proc.time() - aldex_start
rezzies1[sim, "aldex_time"] <- aldex_end[3]
if (!inherits(aldex_model, "try-error")) {
  rezzies1[sim, "aldex_est"] <- aldex_model$`cov1:Est`[test_ind[1]]
  rezzies1[sim, "aldex_p"] <- aldex_model$`cov1:pval`[test_ind[1]]
  rezzies2[sim, "aldex_est"] <- aldex_model$`cov1:Est`[test_ind[2]]
  rezzies2[sim, "aldex_p"] <- aldex_model$`cov1:pval`[test_ind[2]]
  rezzies3[sim, "aldex_est"] <- aldex_model$`cov1:Est`[test_ind[3]]
  rezzies3[sim, "aldex_p"] <- aldex_model$`cov1:pval`[test_ind[3]]
}

print("aldex done")
rezzies <- rbind(rezzies1, rezzies2, rezzies3)
saveRDS(rezzies, paste0(sim_dir_name, "/remaining_setting", setting, "_first_seed", seeds[1], ".rds"))

# ANCOM-BC2
rownames(Y) <- paste0("sample", 1:nrow(Y))
colnames(Y) <- paste0("category", 1:ncol(Y))
ancom_start <- proc.time()
tse_data <- TreeSummarizedExperiment(assays = list("counts" = t(Y)),
                                     colData = data.frame(cov1 = X[, 2],
                                                          cov2 = X[, 3]))
ancom_model <- try({
  suppressMessages(
    ancombc2(data = tse_data,
             assay_name = "counts", 
             p_adj_method = "none",
             fix_formula = "cov1 + cov2",
             prv_cut = 0,
             alpha = 0.05,
             pseudo_sens = FALSE)
  )
})
ancom_end <- proc.time() - ancom_start
rezzies1[sim, "ancom_time"] <- ancom_end[3]
if (!inherits(ancom_model, "try-error")) {
  rezzies1[sim, "ancom_est"] <- ancom_model$res$lfc_cov1[test_ind[1]]
  rezzies1[sim, "ancom_p"] <- ancom_model$res$p_cov1[test_ind[1]]
  rezzies2[sim, "ancom_est"] <- ancom_model$res$lfc_cov1[test_ind[2]]
  rezzies2[sim, "ancom_p"] <- ancom_model$res$p_cov1[test_ind[2]]
  rezzies3[sim, "ancom_est"] <- ancom_model$res$lfc_cov1[test_ind[3]]
  rezzies3[sim, "ancom_p"] <- ancom_model$res$p_cov1[test_ind[3]]
}

print("ancom done")

rezzies <- rbind(rezzies1, rezzies2, rezzies3)
saveRDS(rezzies, paste0(sim_dir_name, "/remaining_setting", setting, "_first_seed", seeds[1], ".rds"))

# CLR t-test with pseudocount of 1 for all counts
clr_start <- proc.time()
Y_pseudo <- Y + 1
log_Y <- log(Y_pseudo)
mean_log_Y <- rowMeans(log_Y)
clr_data <- log_Y - matrix(mean_log_Y, nrow = n, ncol = J, byrow = FALSE)
mod1 <- lm(clr_Y ~ cov1 + cov2, 
           data = data.frame(clr_Y = clr_data[, test_ind[1]], 
                             cov1 = X[, 2], cov2 = X[, 3]))
rezzies1[sim, "clr_est"] <- summary(mod1)$coef[2, 1]
rezzies1[sim, "clr_p"] <- summary(mod1)$coef[2, 4]
mod2 <- lm(clr_Y ~ cov1 + cov2, 
           data = data.frame(clr_Y = clr_data[, test_ind[2]], 
                             cov1 = X[, 2], cov2 = X[, 3]))
rezzies2[sim, "clr_est"] <- summary(mod2)$coef[2, 1]
rezzies2[sim, "clr_p"] <- summary(mod2)$coef[2, 4]
mod3 <- lm(clr_Y ~ cov1 + cov2, 
           data = data.frame(clr_Y = clr_data[, test_ind[3]], 
                             cov1 = X[, 2], cov2 = X[, 3]))
rezzies3[sim, "clr_est"] <- summary(mod3)$coef[2, 1]
rezzies3[sim, "clr_p"] <- summary(mod3)$coef[2, 4]
clr_end <- proc.time() - clr_start
rezzies1[sim, "clr_time"] <- clr_end[3]

print("clr done")

rezzies <- rbind(rezzies1, rezzies2, rezzies3)
saveRDS(rezzies, paste0(sim_dir_name, "/remaining_setting", setting, "_first_seed", seeds[1], ".rds"))

# DESeq2
deseq_start <- proc.time()
dds <- DESeqDataSetFromMatrix(countData = t(Y_pseudo),
                              colData = data.frame(cov1 = X[, 2],
                                                   cov2 = X[, 3]),
                              design = ~ cov1 + cov2)
dds <- try({DESeq(dds, sfType = "poscounts")})
deseq_end <- proc.time() - deseq_start
rezzies1[sim, "deseq_time"] <- deseq_end[3]
if (!inherits(dds, "try-error")) {
  deseq_res <- results(dds, name = "cov1")
  rezzies1[sim, "deseq_est"] <- deseq_res$log2FoldChange[test_ind[1]]
  rezzies1[sim, "deseq_p"] <- deseq_res$pvalue[test_ind[1]]
  rezzies2[sim, "deseq_est"] <- deseq_res$log2FoldChange[test_ind[2]]
  rezzies2[sim, "deseq_p"] <- deseq_res$pvalue[test_ind[2]]
  rezzies3[sim, "deseq_est"] <- deseq_res$log2FoldChange[test_ind[3]]
  rezzies3[sim, "deseq_p"] <- deseq_res$pvalue[test_ind[3]]
}

print("deseq done")

rezzies <- rbind(rezzies1, rezzies2, rezzies3)
saveRDS(rezzies, paste0(sim_dir_name, "/remaining_setting", setting, "_first_seed", seeds[1], ".rds"))

# linda
features <- t(Y)
meta <- data.frame(cov1 = X[, 2], cov2 = X[, 3])
linda_start <- proc.time()
linda_res <- try({linda(feature.dat = features, meta.dat = meta, formula = '~cov1 + cov2',
                        feature.dat.type = "count", is.winsor = FALSE)})
linda_end <- proc.time() - linda_start
rezzies1[sim, "linda_time"] <- linda_end[3]
if (!inherits(linda_res, "try-error")) {
  rezzies1[sim, "linda_est"] <- linda_res$output$cov1$log2FoldChange[test_ind[1]]
  rezzies1[sim, "linda_p"] <- linda_res$output$cov1$pvalue[test_ind[1]]
  rezzies2[sim, "linda_est"] <- linda_res$output$cov1$log2FoldChange[test_ind[2]]
  rezzies2[sim, "linda_p"] <- linda_res$output$cov1$pvalue[test_ind[2]]
  rezzies3[sim, "linda_est"] <- linda_res$output$cov1$log2FoldChange[test_ind[3]]
  rezzies3[sim, "linda_p"] <- linda_res$output$cov1$pvalue[test_ind[3]]
}

print("linda done")

rezzies <- rbind(rezzies1, rezzies2, rezzies3)
saveRDS(rezzies, paste0(sim_dir_name, "/remaining_setting", setting, "_first_seed", seeds[1], ".rds"))

# maaslin
dat <- data.frame(cov1 = X[, 2], cov2 = X[, 3])
rownames(Y) <- paste0("sample", 1:nrow(Y))
rownames(dat) <- rownames(Y)
file <- paste0("~/radEmu/results/sims/maaslin_output/maaslin_output_setting", setting, "_seed", seeds[sim])
#file <- "maaslin_temp"
maaslin_start <- proc.time()
maaslin_res <- try({maaslin3(input_data = Y, input_metadata = dat,
                             formula = ~cov1 + cov2, output = file, 
                             plot_summary_plot = FALSE, plot_associations = FALSE,
                             save_models = FALSE)})
maaslin_end <- proc.time() - maaslin_start
rezzies1[sim, "maaslin_time"] <- maaslin_end[3]
if (!inherits(maaslin_res, "try-error") & file.exists(paste0(file, "/all_results.tsv"))) {
  maas_res <- read.delim(paste0(file, "/all_results.tsv")) 
  res1 <- maas_res %>%
    filter(feature == paste0("taxon", test_ind[1]), model == "abundance", value == "cov1")
  if (nrow(res1) > 0) {
    rezzies1[sim, "maaslin_est"] <- res1$coef
    rezzies1[sim, "maaslin_p"] <- res1$pval_individual
  }
  res2 <- maas_res %>%
    filter(feature == paste0("taxon", test_ind[2]), model == "abundance", value == "cov1")
  if (nrow(res2) > 0) {
    rezzies2[sim, "maaslin_est"] <- res2$coef
    rezzies2[sim, "maaslin_p"] <- res2$pval_individual
  }
  res3 <- maas_res %>%
    filter(feature == paste0("taxon", test_ind[3]), model == "abundance", value == "cov1")
  if (nrow(res3) > 0) {
    rezzies3[sim, "maaslin_est"] <- res3$coef
    rezzies3[sim, "maaslin_p"] <- res3$pval_individual
  }
}

print("maaslin done")

rezzies <- rbind(rezzies1, rezzies2, rezzies3)
saveRDS(rezzies, paste0(sim_dir_name, "/remaining_setting", setting, "_first_seed", seeds[1], ".rds"))

# radEmu 
emu_start <- proc.time()
emu_model <- try({
  emuFit(Y = Y,
         X = X,
         test_kj = data.frame(k = 2, j = test_ind),
         tau = 2,
         B_null_tol = 0.005,
         tolerance = 0.005,
         constraint_tol = 0.001,
         return_wald_p = TRUE,
         match_row_names = FALSE, 
         unobserved_taxon_error = FALSE)})
emu_end <- proc.time() - emu_start
rezzies1[sim, "rademu_time"] <- emu_end[3]
if (!inherits(emu_model, "try-error")) {
  rezzies1[sim, "rademu_est"] <- emu_model$coef$estimate[test_ind[1]]
  rezzies1[sim, "wald_p"] <- emu_model$coef$wald_p[test_ind[1]]
  rezzies1[sim, "score_p"] <- emu_model$coef$pval[test_ind[1]]
  rezzies2[sim, "rademu_est"] <- emu_model$coef$estimate[test_ind[2]]
  rezzies2[sim, "wald_p"] <- emu_model$coef$wald_p[test_ind[2]]
  rezzies2[sim, "score_p"] <- emu_model$coef$pval[test_ind[2]]
  rezzies3[sim, "rademu_est"] <- emu_model$coef$estimate[test_ind[3]]
  rezzies3[sim, "wald_p"] <- emu_model$coef$wald_p[test_ind[3]]
  rezzies3[sim, "score_p"] <- emu_model$coef$pval[test_ind[3]]
}

print("radEmu done")

rezzies <- rbind(rezzies1, rezzies2, rezzies3)
saveRDS(rezzies, paste0(sim_dir_name, "/remaining_setting", setting, "_first_seed", seeds[1], ".rds"))

# locom
colnames(Y) <- paste0("taxon", 1:ncol(Y))
locom_start <- proc.time()
locom_res <- try({locom(otu.table = Y, Y = X[, 2:3], filter.thresh = 0)})
locom_end <- proc.time() - locom_start
rezzies1[sim, "locom_time"] <- locom_end[3]
if (!inherits(locom_res, "try-error")) {
  
  ind <- try({which(colnames(locom_res$effect.size) == paste0("taxon", test_ind[1]))})
  if (!inherits(ind, "try-error") && length(ind) == 1) {
    rezzies1[sim, "locom_est"] <- locom_res$effect.size[ind]
    rezzies1[sim, "locom_p"] <- locom_res$p.otu[ind]
  }
  
  ind <- try({which(colnames(locom_res$effect.size) == paste0("taxon", test_ind[2]))})
  if (!inherits(ind, "try-error") && length(ind) == 1) {
    rezzies2[sim, "locom_est"] <- locom_res$effect.size[ind]
    rezzies2[sim, "locom_p"] <- locom_res$p.otu[ind]
  }
  
  ind <- try({which(colnames(locom_res$effect.size) == paste0("taxon", test_ind[3]))})
  if (!inherits(ind, "try-error") && length(ind) == 1) {
    rezzies3[sim, "locom_est"] <- locom_res$effect.size[ind]
    rezzies3[sim, "locom_p"] <- locom_res$p.otu[ind]
  }
}

print("locom done")

rezzies <- rbind(rezzies1, rezzies2, rezzies3)
saveRDS(rezzies, paste0(sim_dir_name, "/remaining_setting", setting, "_first_seed", seeds[1], ".rds"))
