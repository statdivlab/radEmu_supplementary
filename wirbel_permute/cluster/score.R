library(tidyverse)
library(radEmu)

args <- commandArgs(trailingOnly = FALSE)
if (length(args) == 0) {
  batch <- 1
} else {
  arg <- args[length(args)]
  batch <- abs(readr::parse_number(arg))
}

test_per_script <- 20
script_per_perm <- 43
tests <- ((batch - 1) %% script_per_perm) * test_per_script + 1:test_per_script
perm <- floor((batch - 1) / script_per_perm) + 1

if (tests[1] == 841) {
  tests <- 841:845
}

# get X and Y matrices
Y <- readRDS("radEmu/wirbel_permute/data/Y.rds")
X_perm <- readRDS(paste0("radEmu/wirbel_permute/data/X_perm", perm, ".rds"))

# get estimates under the alternate
B <- readRDS(paste0("radEmu/wirbel_permute/results/radEmu_fit_B_perm", perm, ".rds"))

# save results
df <- data.frame(perm = perm,
                 tests = tests,
                 tax = colnames(Y)[tests],
                 score_p = NA)

print(paste0("permutation number ", perm))
# run robust score tests 
for (i in 1:nrow(df)) {
  print(i)
  
  emuRes <- emuFit(Y = Y,
         refit = FALSE,
         B = B, 
         X = X_perm,
         tolerance = 0.01,
         test_kj = data.frame(k = 2,
                              j = tests[i]),
         return_wald_p = FALSE,
         compute_cis = FALSE,
         run_score_tests = TRUE,
         verbose = TRUE,
         constraint_tol = 1e-3)
  
  df$score_p[i] <- emuRes$coef$pval[tests[i]]
}

# save results
saveRDS(df, paste0("radEmu/wirbel_permute/results/score_perm", perm, "first_test", tests[1], ".rds"))

