library(tidyverse)
library(MicrobiomeStat)
library(maaslin3)
library(DESeq2)

files <- list.files("fig2_updated/additional_results", full.names=T)
files_other <- list.files("fig2_updated/results", full.names=T)

results <- vector(length(files), mode = "list")
results_other <- vector(length(files), mode = "list")
counter <- 1
for(file in files){
  res_so_far <- readRDS(file)
  results[[counter]] <- res_so_far
  counter <- counter + 1
}
counter <- 1
for (file in files_other) {
  res_so_far <- readRDS(file)
  results_other[[counter]] <- res_so_far
  counter <- counter + 1
}

results <- do.call(bind_rows, results) %>% as_tibble
results_other <- do.call(bind_rows, results_other) %>% as_tibble
results <- full_join(results, results_other)

results %>% group_by(n, J, distn) %>% 
  summarise(mean(is.na(deseq_p)), mean(is.na(linda_p)), mean(is.na(maaslin_p))) %>%
  filter(`mean(is.na(deseq_p))` >= 0.05 |  
           `mean(is.na(linda_p))` >= 0.05 | 
           `mean(is.na(maaslin_p))` >= 0.05) 

# investigate maaslin NAs compared to sum(Y[, J/2] > 0)
Y_test_sums <- data.frame(seed = 1:500,
                      prev10 = NA,
                      prev50 = NA,
                      prev250 = NA)
#parameters determining form of simulations
n <- 10
Js <- c(250, 50, 10)
distns <- "ZINB"

# generate B matrices
set.seed(9443)
B_list <- lapply(Js,
                 function(j)
                   radEmu:::get_sim_bs(j))
n <- 10
X <- cbind(intercept = 1, cov = rep(0:1, each = n/2))
for (seed in 1:500) {
  print(seed)
  set.seed(seed)
  Y <- generate_test_data(n = n,
                          J = 10,
                          b0 = B_list[[3]]$b0,
                          b1 = B_list[[3]]$b1,
                          distn = distn,
                          zinb_size = 5,
                          zinb_zero_prop = 0.6,
                          mean_count_before_ZI = 50)
  Y_test_sums[seed, 2] <- sum(Y[, 10/2] > 0)
  
  set.seed(seed)
  Y <- generate_test_data(n = n,
                          J = 50,
                          b0 = B_list[[2]]$b0,
                          b1 = B_list[[2]]$b1,
                          distn = distn,
                          zinb_size = 5,
                          zinb_zero_prop = 0.6,
                          mean_count_before_ZI = 50)
  Y_test_sums[seed, 3] <- sum(Y[, 50/2] > 0)
  
  set.seed(seed)
  Y <- generate_test_data(n = n,
                          J = 250,
                          b0 = B_list[[1]]$b0,
                          b1 = B_list[[1]]$b1,
                          distn = distn,
                          zinb_size = 5,
                          zinb_zero_prop = 0.6,
                          mean_count_before_ZI = 50)
  Y_test_sums[seed, 4] <- sum(Y[, 250/2] > 0)
}
colSums(Y_test_sums[, 2:4] == 0)
colSums(Y_test_sums[, 2:4] < 2)
results %>% filter(n == 10, distn == "ZINB") %>% group_by(J) %>%
  summarise(sum(is.na(maaslin_p)))
# this does not explain it! instead look at specific cases... 

results %>% filter(is.na(maaslin_p)) %>% group_by(J, is.na(maaslin_est)) %>%
  summarise(n())
results %>% filter(is.na(maaslin_p), J == 10) %>% head()
Y_test_sums[c(1, 10, 12, 15, 112, 117), 2]

set.seed(12)
Y <- generate_test_data(n = 10,
                        J = 10,
                        b0 = B_list[[3]]$b0,
                        b1 = B_list[[3]]$b1,
                        distn = "ZINB",
                        zinb_size = 5,
                        zinb_zero_prop = 0.6,
                        mean_count_before_ZI = 50)
dat <- data.frame(cov = X[, 2])
rownames(Y) <- paste0("sample", 1:nrow(Y))
rownames(dat) <- rownames(Y)
file <- paste0("maaslin_test")
maaslin_res <- try({maaslin3(input_data = Y, input_metadata = dat,
                             formula = ~cov, output = file, 
                             plot_summary_plot = FALSE, plot_associations = FALSE,
                             save_models = FALSE)})
# maaslin either excludes feature because it appears in 0-1 samples
# or it fails to converge, with message "Fitting error (NA p-value returned from fitting procedure)"

# LinDA
results %>% group_by(n, J, distn) %>% summarise(sum(is.na(linda_p)))
results %>% filter(is.na(linda_p), n == 250, J == 10, distn == "ZINB") %>% head()
X <- cbind(intercept = 1, cov = rep(0:1, each = 250/2))
set.seed(1)
Y <- generate_test_data(n = 250,
                        J = 10,
                        b0 = B_list[[3]]$b0,
                        b1 = B_list[[3]]$b1,
                        distn = "ZINB",
                        zinb_size = 5,
                        zinb_zero_prop = 0.6,
                        mean_count_before_ZI = 50)
features <- t(Y)
meta <- data.frame(cov = X[, 2])
linda_res <- try({linda(feature.dat = features, meta.dat = meta, formula = '~cov',
                        feature.dat.type = "count", is.winsor = FALSE)})
linda_res$output$cov[10/2, ]
# for linda, set is.winsor = FALSE and rerun

# DESeq2
results %>% group_by(n, J, distn) %>% summarise(mean(is.na(deseq_p)))

results %>% filter(is.na(deseq_p)) %>% head()
X <- cbind(intercept = 1, cov = rep(0:1, each = 250/2))
set.seed(1)
Y <- generate_test_data(n = 250,
                        J = 10,
                        b0 = B_list[[3]]$b0,
                        b1 = B_list[[3]]$b1,
                        distn = "Poisson",
                        zinb_size = 5,
                        zinb_zero_prop = 0.6,
                        mean_count_before_ZI = 50)
Y_pseudo <- Y + 1
dds <- DESeqDataSetFromMatrix(countData = t(Y_pseudo),
                              colData = data.frame(X = as.factor(X[, 2])),
                              design = ~ X)
dds <- try({DESeq(dds)})

results %>% filter(is.na(deseq_p), n == 10, J == 10, distn == "Poisson") %>% head()
X <- cbind(intercept = 1, cov = rep(0:1, each = 10/2))

set.seed(8)
Y <- generate_test_data(n = 10,
                        J = 10,
                        b0 = B_list[[3]]$b0,
                        b1 = B_list[[3]]$b1,
                        distn = "Poisson",
                        zinb_size = 5,
                        zinb_zero_prop = 0.6,
                        mean_count_before_ZI = 50)
Y_pseudo <- Y + 1
dds <- DESeqDataSetFromMatrix(countData = t(Y_pseudo),
                              colData = data.frame(X = as.factor(X[, 2])),
                              design = ~ X)
dds <- try({DESeq(dds)})

results %>% filter(is.na(deseq_p), n == 10, J == 250, distn == "ZINB") %>% head()
X <- cbind(intercept = 1, cov = rep(0:1, each = 10/2))

set.seed(7)
Y <- generate_test_data(n = 10,
                        J = 250,
                        b0 = B_list[[1]]$b0,
                        b1 = B_list[[1]]$b1,
                        distn = "ZINB",
                        zinb_size = 5,
                        zinb_zero_prop = 0.6,
                        mean_count_before_ZI = 50)
Y_pseudo <- Y + 1
dds <- DESeqDataSetFromMatrix(countData = t(Y_pseudo),
                              colData = data.frame(X = as.factor(X[, 2])),
                              design = ~ X)
dds <- try({DESeq(dds)})
deseq_res <- results(dds, name = "X_1_vs_0")
deseq_res[J/2, ]
