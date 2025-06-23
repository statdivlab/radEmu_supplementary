library(radEmu)
library(tidyverse)

# run script to get function "generate_test_data" within folder "radEmu/scripts"
source("simulation-data-generation/generate_test_data.R")

# set hyperparams
n <- 250
distn <- "ZINB"
Js <- c(250, 50)
nsim <- 50

# results
alt_time <- data.frame(seed = 1:50, time50 = NA, time250 = NA)
null_time <- data.frame(seed = 1:50, time50 = NA, time250 = NA)

# generate B matrices
set.seed(9443)
B_list <- lapply(Js,
                 function(j)
                   radEmu:::get_sim_bs(j))
X <- cbind(intercept = 1, cov = rep(0:1, each = n/2))

# iterate over seeds
for (seed in 1:nsim) {
  
  print(seed)
  
  # J = 50
  set.seed(seed)
  Y50 <- generate_test_data(n = n,
                          J = Js[2],
                          b0 = B_list[[2]]$b0,
                          b1 = B_list[[2]]$b1,
                          distn = distn,
                          zinb_size = 5,
                          zinb_zero_prop = 0.6,
                          mean_count_before_ZI = 50)
  t0 <- proc.time()
  emu_model <- try({
    emuFit(Y = Y50,
           X = X,
           tau = 2,
           B_null_tol = 0.005,
           tolerance = 0.005,
           constraint_tol = 0.001,
           unobserved_taxon_error = FALSE,
           run_score_tests = FALSE,
           compute_cis = FALSE,
           match_row_names = FALSE)})
  t1 <- proc.time() - t0
  alt_time[seed, "time50"] <- t1[3]
  t0 <- proc.time()
  emu_score <- try({
    emuFit(Y = Y50,
           X = X,
           test_kj = data.frame(k = 2, j = 25),
           tau = 2,
           B_null_tol = 0.005,
           tolerance = 0.005,
           constraint_tol = 0.001,
           unobserved_taxon_error = FALSE,
           fitted_model = emu_model,
           refit = FALSE,
           compute_cis = FALSE,
           match_row_names = FALSE)})
  t1 <- proc.time() - t0
  null_time[seed, "time50"] <- t1[3]
  
  # J = 250
  set.seed(seed)
  Y250 <- generate_test_data(n = n,
                            J = Js[1],
                            b0 = B_list[[1]]$b0,
                            b1 = B_list[[1]]$b1,
                            distn = distn,
                            zinb_size = 5,
                            zinb_zero_prop = 0.6,
                            mean_count_before_ZI = 50)
  t0 <- proc.time()
  emu_model <- try({
    emuFit(Y = Y250,
           X = X,
           tau = 2,
           B_null_tol = 0.005,
           tolerance = 0.005,
           constraint_tol = 0.001,
           unobserved_taxon_error = FALSE,
           run_score_tests = FALSE,
           compute_cis = FALSE,
           match_row_names = FALSE)})
  t1 <- proc.time() - t0
  alt_time[seed, "time250"] <- t1[3]
  t0 <- proc.time()
  emu_score <- try({
    emuFit(Y = Y250,
           X = X,
           test_kj = data.frame(k = 2, j = 125),
           tau = 2,
           B_null_tol = 0.005,
           tolerance = 0.005,
           constraint_tol = 0.001,
           unobserved_taxon_error = FALSE,
           fitted_model = emu_model,
           refit = FALSE,
           compute_cis = FALSE,
           match_row_names = FALSE)})
  t1 <- proc.time() - t0
  null_time[seed, "time250"] <- t1[3]
}

median(alt_time$time50)
median(alt_time$time250)
median(null_time$time50)
median(null_time$time250)
