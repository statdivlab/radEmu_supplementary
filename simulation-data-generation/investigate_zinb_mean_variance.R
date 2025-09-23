library(tidyverse)

source("simulation-data-generation/generate_test_data.R")

n <- 250
J <- 250
set.seed(9443)
B <- radEmu:::get_sim_bs(J)
B <- rbind(B$b0, B$b1, sample(B$b1, size = J, replace = FALSE))

res_logmean <- matrix(NA, nrow = 500, ncol = n * J)
res_Y <- matrix(NA, nrow = 500, ncol = n * J)

for (sim in 1:500) {
  print(sim)
  
  X <- cbind(intercept = 1, 
             cov1 = runif(n, min = 0, max = 1),
             cov2 = rnorm(n, mean = 0, sd = 2))
  
  res <- generate_test_data_include_delta(n = n,
                                          J = J,
                                          B = B, 
                                          X = X,
                                          distn = "ZINB",
                                          zinb_size = 5,
                                          zinb_zero_prop = 0.6,
                                          mean_count_before_ZI = 50)
  
  res_logmean[sim, ] <- as.vector(res$log_means)
  res_Y[sim, ] <- as.vector(res$Y)
  
}
hist(as.vector(res_logmean))
summary(exp(as.vector(res_logmean)))
summary(as.vector(res_logmean))

hist(log(as.vector(res_Y)))
summary(as.vector(res_Y))

# theoretically, expect var(Y) = .4mu + .32mu^2 

fn <- function(x) {
  return(.4 * x + .32 * x^2)
}
x_vals <- exp(seq(from = -100, to = 13, length.out = 5000))
df <- data.frame(x = rep(x_vals, 2), 
                 y = c(fn(x_vals), x_vals), 
                 Distribution = rep(c("ZINB", "Poisson"), each = 5000))
ggplot(df, aes(x = x, y = y, color = Distribution)) +
  geom_line() + 
  theme_bw() + 
  labs(x = "Mean", y = "Variance")
ggsave("simulation-data-generation/full.png")
df %>% filter(x < 1280) %>% # 99% of mu fall in this interval 
  ggplot(aes(x = x, y = y, color = Distribution)) +
  geom_line() + 
  theme_bw() + 
  labs(x = "Mean (up to 99th quantile of means)", y = "Variance")
df %>% filter(x < 155) %>% # 95% of mu fall in this interval 
  ggplot(aes(x = x, y = y, color = Distribution)) +
  geom_line() + 
  theme_bw() + 
  labs(x = "Mean (up to 95th quantile of means)", y = "Variance") 
ggsave("simulation-data-generation/95th_quantile.png")
df %>% filter(x < 8) %>% # 75% of mu fall in this interval 
  ggplot(aes(x = x, y = y, color = Distribution)) +
  geom_line() + 
  theme_bw() + 
  labs(x = "Mean (up to 75th quantile of means)", y = "Variance") 
ggsave("simulation-data-generation/75th_quantile.png")
df %>% filter(x < .7) %>% # 50% of mu fall in this interval 
  ggplot(aes(x = x, y = y, color = Distribution)) +
  geom_line() + 
  theme_bw() + 
  labs(x = "Mean (up to 50th quantile of means)", y = "Variance") 
ggsave("simulation-data-generation/50th_quantile.png")

# check on variance of simulated data 
x_sims <- exp(seq(from = -100, to = 13, length.out = 5000))
sim_res <- matrix(NA, 5000, 100)
for (i in 1:5000) {
  sim_res[i, ] <- stats::rnbinom(100, mu = x_sims[i], size= 5)*(1 - stats::rbinom(100,1,prob = .6))
}
sim_vars <- apply(sim_res, 1, var)
sim_df <- data.frame(x = x_sims, y = sim_vars)
ggplot(sim_df, aes(x = x, y = y)) + geom_point() + 
  geom_line(data = df %>% filter(Distribution == "ZINB"), aes(x = x, y = y), color = "red")
sim_df %>% filter(x < 1280) %>% 
  ggplot(aes(x = x, y = y)) + geom_point() + 
  geom_line(data = df %>% filter(Distribution == "ZINB", x < 1280), aes(x = x, y = y), color = "red")
sim_df %>% filter(x < 155) %>% 
  ggplot(aes(x = x, y = y)) + geom_point() + 
  geom_line(data = df %>% filter(Distribution == "ZINB", x < 155), aes(x = x, y = y), color = "red")
sim_df %>% filter(x < 8) %>% 
  ggplot(aes(x = x, y = y)) + geom_point() + 
  geom_line(data = df %>% filter(Distribution == "ZINB", x < 8), aes(x = x, y = y), color = "red")
sim_df %>% filter(x < .7) %>% 
  ggplot(aes(x = x, y = y)) + geom_point() + 
  geom_line(data = df %>% filter(Distribution == "ZINB", x < .7), aes(x = x, y = y), color = "red")
