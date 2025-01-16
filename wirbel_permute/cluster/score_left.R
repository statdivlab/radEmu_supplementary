library(tidyverse)
library(radEmu)

args <- commandArgs(trailingOnly = FALSE)
if (length(args) == 0) {
  batch <- 1
} else {
  arg <- args[length(args)]
  batch <- abs(readr::parse_number(arg))
}

test_per_script <- 5
inds <- (batch - 1) * test_per_script + 1:test_per_script
df_remain <- readRDS("radEmu/wirbel_permute/results/df_remain.rds")
tests_left <- df_remain[inds, ]

df <- data.frame(perm = tests_left$perm,
                 tests = tests_left$tests,
                 tax = tests_left$tax,
                 score_p = NA)

# get Y matrix 
Y <- readRDS("radEmu/wirbel_permute/data/Y.rds")

print(paste0("batch number ", batch))
# run robust score tests 
for (i in 1:nrow(df)) {
  # get X matrix 
  X_perm <- readRDS(paste0("radEmu/wirbel_permute/data/X_perm", df$perm[i], ".rds"))
  
  # get estimates under the alternate
  B <- readRDS(paste0("radEmu/wirbel_permute/results/radEmu_fit_B_perm", df$perm[i], ".rds"))
  
  print(i)
  
  emuRes <- emuFit(Y = Y,
                   refit = FALSE,
                   B = B, 
                   X = X_perm,
                   tolerance = 0.01,
                   test_kj = data.frame(k = 2,
                                        j = df$tests[i]),
                   return_wald_p = FALSE,
                   compute_cis = FALSE,
                   run_score_tests = TRUE,
                   verbose = TRUE,
                   constraint_tol = 1e-3)
  
  df$score_p[i] <- emuRes$coef$pval[df$tests[i]]
}

# save results
saveRDS(df, paste0("radEmu/wirbel_permute/results/rerun_score_batch", batch, ".rds"))