library(tidyverse)

source("fig456/0-process.R")

saveRDS(Y, "wirbel_permute/data/Y.rds")
saveRDS(X, "wirbel_permute/data/X.rds")

# permute rows of X 100 times
n_perm <- 100
for (i in 1:n_perm) {
  set.seed(i)
  X_perm <- X[sample(1:nrow(X)), ]
  saveRDS(X_perm, paste0("wirbel_permute/data/X_perm", i, ".rds"))
}
