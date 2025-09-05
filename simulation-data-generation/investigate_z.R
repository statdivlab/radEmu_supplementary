library(radEmu)
library(tidyverse)

source("simulation-data-generation/generate_test_data.R")

n <- 50
J <- 250

Bs <- radEmu:::get_sim_bs(J)

# Poisson data 
m <- 50
pois_res <- matrix(nrow = 500, ncol = 50)
pois_z_res <- matrix(nrow = 500, ncol = 50)

for (i in 1:500) {
  set.seed(i)
  #print(i)
  
  res <- generate_test_data_testing(n = n,
                                    J = J,
                                    b0 = Bs$b0,
                                    b1 = Bs$b1,
                                    distn = "Poisson",
                                    zinb_size = 5,
                                    zinb_zero_prop = 0.6,
                                    mean_count_before_ZI = m)
  Y <- res$Y
  z <- res$z
  
  pois_res[i, ] <- rowSums(Y)
  pois_z_res[i, ] <- z
  
}
mean(pois_res)
J * m * exp(1/2)
hist(pois_z_res)
pois_z_df <- data.frame(z = as.vector(pois_z_res),
                        it = rep(1:500, 50))
ggplot(pois_z_df, aes(x = it, y = z)) + 
  geom_point(alpha = 0.2) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle("simulated z's over 500 iterations, Poisson")

# ZINB data 
m <- 50
zinb_res <- matrix(nrow = 500, ncol = 50)
zinb_z_res <- matrix(nrow = 500, ncol = 50)

for (i in 1:500) {
  set.seed(i)
  #print(i)
  
  res <- generate_test_data_testing(n = n,
                            J = J,
                            b0 = Bs$b0,
                            b1 = Bs$b1,
                            distn = "ZINB",
                            zinb_size = 5,
                            zinb_zero_prop = 0.6,
                            mean_count_before_ZI = m)
  Y <- res$Y
  z <- res$z
  
  zinb_res[i, ] <- rowSums(Y)
  zinb_z_res[i, ] <- z
  
}
mean(zinb_res)
J * m * exp(1/2) * (1 - 0.6)
zinb_z_df <- data.frame(z = as.vector(zinb_z_res),
                        it = rep(1:500, 50))
ggplot(zinb_z_df, aes(x = it, y = z)) + 
  geom_point(alpha = 0.2) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle("simulated z's over 500 iterations, ZINB")
