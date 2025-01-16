library(tidyverse)
library(radEmu)

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

# radEmu estimation and Wald tests
full_fit_wald <- emuFit(Y = Y,
                        X = X_perm,
                        tolerance = 0.01,
                        test_kj = data.frame(k = 2,
                                             j = 1:ncol(Y)),
                        return_wald_p = TRUE,
                        compute_cis = TRUE,
                        run_score_tests = FALSE,
                        verbose = TRUE)
saveRDS(full_fit_wald$coef, paste0("radEmu/wirbel_permute/results/radEmu_fit_wald_perm", perm, ".rds"))
saveRDS(full_fit_wald$B, paste0("radEmu/wirbel_permute/results/radEmu_fit_B_perm", perm, ".rds"))
