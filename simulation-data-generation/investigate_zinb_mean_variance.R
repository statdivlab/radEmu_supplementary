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
mean(exp(as.vector(res_logmean)) <= 1.875) # 1.875 is the value at which the two functions cross

hist(log(as.vector(res_Y)))
summary(as.vector(res_Y))

# theoretically, expect var(Y) = .4mu + .32mu^2 
library(latex2exp)

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
  labs(x = TeX("$\\mu_{ij}:\\; X_i \\beta_j + z_i + \\delta_j$"), 
       y = "Variance")
ggsave("simulation-data-generation/full.png")
df %>% filter(x < 1280) %>% # 99% of mu fall in this interval 
  ggplot(aes(x = x, y = y, color = Distribution)) +
  geom_line() + 
  theme_bw() + 
  labs(x = expression(atop(mu[ij] * ": " ~ X[i] * beta[j] + z[i] + delta[j],
                           "0–99th quantiles shown")), 
       y = "Variance")
p95 <- df %>% filter(x < 155) %>% # 95% of mu fall in this interval 
  ggplot(aes(x = x, y = y, color = Distribution)) +
  geom_line() + 
  theme_bw() + 
  labs(x = expression(atop(mu[ij] * ": " ~ X[i] * beta[j] + z[i] + delta[j],
                           "0–95th quantiles shown")), 
       y = "Variance")
ggsave("simulation-data-generation/95th_quantile.png")
p75 <- df %>% filter(x < 8) %>% # 75% of mu fall in this interval 
  ggplot(aes(x = x, y = y, color = Distribution)) +
  geom_line() + 
  theme_bw() + 
  labs(x = expression(atop(mu[ij] * ": " ~ X[i] * beta[j] + z[i] + delta[j],
                           "0–75th quantiles shown")), 
       y = "Variance")
ggsave("simulation-data-generation/75th_quantile.png")
p50 <- df %>% filter(x < .7) %>% # 50% of mu fall in this interval 
  ggplot(aes(x = x, y = y, color = Distribution)) +
  geom_line() + 
  theme_bw() + 
  labs(x = expression(atop(mu[ij] * ": " ~ X[i] * beta[j] + z[i] + delta[j],
                           "0–50th quantiles shown")), 
       y = "Variance")
ggsave("simulation-data-generation/50th_quantile.png")

library(ggpubr)
ggarrange(p95, p75, p50, nrow = 1, common.legend = TRUE, legend = "bottom")
ggsave("simulation-data-generation/mean_variance.pdf", width = 8, height = 5)

# fn value at quantiles
quants <- quantile(exp(as.vector(res_logmean)), c(0.95, 0.75, 0.5))
fn(quants[1])/quants[1]
fn(quants[2])/quants[2]
fn(quants[3])/quants[3]

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
