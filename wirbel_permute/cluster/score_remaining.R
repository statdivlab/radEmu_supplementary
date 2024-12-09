library(tidyverse)
library(radEmu)

args <- commandArgs(trailingOnly = FALSE)
if (length(args) == 0) {
  batch <- 1
} else {
  arg <- args[length(args)]
  batch <- abs(readr::parse_number(arg))
}

print(batch)
df_run_again <- readRDS("radEmu/wirbel_permute/results/df_run_again.rds")
df <- df_run_again[batch, ]

# get Y matrix 
Y <- readRDS("radEmu/wirbel_permute/data/Y.rds")

# get X matrix 
X_perm <- readRDS(paste0("radEmu/wirbel_permute/data/X_perm", df$perm, ".rds"))

# get estimates under the alternate
B <- readRDS(paste0("radEmu/wirbel_permute/results/radEmu_fit_B_perm", df$perm, ".rds"))

emuRes <- emuFit(Y = Y,
                 refit = FALSE,
                 B = B, 
                 X = X_perm,
                 tolerance = 0.01,
                 test_kj = data.frame(k = 2,
                                      j = df$tests),
                 return_wald_p = FALSE,
                 compute_cis = FALSE,
                 run_score_tests = TRUE,
                 verbose = TRUE,
                 constraint_tol = 1e-3)

df$score_p <- emuRes$coef$pval[df$tests]

# save results
saveRDS(df, paste0("radEmu/wirbel_permute/results/run_again_score_batch", batch, ".rds"))