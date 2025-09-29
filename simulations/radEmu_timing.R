library(tidyverse)

#parameters determining form of simulations
n <- 250
Js <- c(250, 50)
distn <- "ZINB"
alt <- 0

# generate B matrices
set.seed(9443)
B_list <- lapply(Js,
                 function(j)
                   radEmu:::get_sim_bs(j))
b0 <- B_list[[2]]$b0
b1 <- B_list[[2]]$b1
b2 <- sample(b1, size = J, replace = FALSE)
B50 <- rbind(b0, b1, b2)
b0 <- B_list[[1]]$b0
b1 <- B_list[[1]]$b1
b2 <- sample(b1, size = J, replace = FALSE)
B250 <- rbind(b0, b1, b2)

X <- cbind(intercept = 1, 
           cov1 = runif(n, min = 0, max = 1),
           cov2 = rnorm(n, mean = 0, sd = 2))

nsim <- 50

time_sim <- data.frame(seed = 1:nsim,
                       est_time50 = NA,
                       score_time50 = NA,
                       est_time250 = NA,
                       score_time250 = NA)

# run simulations
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)
  
  # J = 50
  J <- 50
  Y <- generate_test_data_include_delta(n = n,
                                        J = J,
                                        B = B50, 
                                        X = X,
                                        distn = distn,
                                        zinb_size = 5,
                                        zinb_zero_prop = 0.6,
                                        mean_count_before_ZI = 50)$Y
  
  est_start <- proc.time()
  est <- emuFit(Y = Y,
         X = X,
         test_kj = data.frame(k = 2, j = J / 2),
         tau = 2,
         B_null_tol = 0.005,
         tolerance = 0.005,
         constraint_tol = 0.001,
         return_wald_p = TRUE,
         match_row_names = FALSE, 
         unobserved_taxon_error = FALSE,
         run_score_tests = FALSE)
  est_end <- proc.time() - est_start
  score_start <- proc.time()
  test <- emuFit(Y = Y,
                X = X,
                test_kj = data.frame(k = 2, j = J / 2),
                tau = 2,
                B_null_tol = 0.005,
                tolerance = 0.005,
                constraint_tol = 0.001,
                match_row_names = FALSE, 
                unobserved_taxon_error = FALSE,
                run_score_tests = TRUE,
                fitted_model = est,
                refit = FALSE,
                compute_cis = FALSE)
  score_end <- proc.time() - score_start
  time_sim$est_time50[sim] <- est_end[3]
  time_sim$score_time50[sim] <- score_end[3]
  
  # J = 250
  J <- 250
  Y <- generate_test_data_include_delta(n = n,
                                        J = J,
                                        B = B250, 
                                        X = X,
                                        distn = distn,
                                        zinb_size = 5,
                                        zinb_zero_prop = 0.6,
                                        mean_count_before_ZI = 50)$Y
  est_start <- proc.time()
  est <- emuFit(Y = Y,
                X = X,
                test_kj = data.frame(k = 2, j = J / 2),
                tau = 2,
                B_null_tol = 0.005,
                tolerance = 0.005,
                constraint_tol = 0.001,
                match_row_names = FALSE, 
                unobserved_taxon_error = FALSE,
                run_score_tests = FALSE)
  est_end <- proc.time() - est_start
  score_start <- proc.time()
  test <- emuFit(Y = Y,
                 X = X,
                 test_kj = data.frame(k = 2, j = J / 2),
                 tau = 2,
                 B_null_tol = 0.005,
                 tolerance = 0.005,
                 constraint_tol = 0.001,
                 match_row_names = FALSE, 
                 unobserved_taxon_error = FALSE,
                 run_score_tests = TRUE,
                 fitted_model = est,
                 refit = FALSE,
                 compute_cis = FALSE)
  score_end <- proc.time() - score_start
  time_sim$est_time250[sim] <- est_end[3]
  time_sim$score_time250[sim] <- score_end[3]
  
}

median(time_sim$est_time50)
median(time_sim$est_time250)
median(time_sim$score_time50)
median(time_sim$score_time250)
