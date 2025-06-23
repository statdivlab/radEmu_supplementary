library(tidyverse)
library(IFAA)
library(SummarizedExperiment)

# get permutation
args <- commandArgs(trailingOnly = FALSE)
if (length(args) == 0) {
  perm <- 1
} else {
  arg <- args[length(args)]
  perm <- abs(readr::parse_number(arg))
}

# get X and Y matrices
Y <- readRDS("radEmu/wirbel_permute/data/Y.rds")
X_perm <- readRDS(paste0("radEmu/wirbel_permute/data/X_perm", perm, ".rds"))

# save results 
df <- data.frame(perm = perm,
                 ind = 1:ncol(Y),
                 tax = colnames(Y),
                 ifaa_est = NA,
                 ifaa_p = NA,
                 ifaa_tax = NA)

# IFAA
experimental_data <- SummarizedExperiment(assays = list("shotgun" = t(Y)),
                                          colData= X_perm) 
set.seed(123) # For full reproducibility
ifaa_results <- IFAA(experiment_dat = experimental_data,
                       testCov="groupCRC",
                       ctrlCov=colnames(X_perm)[-c(1,2)])
df$ifaa_tax <- ifaa_results$full_results$taxon
df$ifaa_est <- ifaa_results$full_results$estimate
df$ifaa_p <- ifaa_results$full_results$unadj.p.value

# save results 
saveRDS(df, paste0("radEmu/wirbel_permute/results/ifaa_perm", perm, ".rds"))
