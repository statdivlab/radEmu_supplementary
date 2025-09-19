library(tidyverse)

files <- list.files("old_simulations/fig3_updated/results", full.names=T)

results <- vector(length(files), mode = "list")
counter <- 1
for(file in files){
  res_so_far <- readRDS(file)
  results[[counter]] <- res_so_far
  counter <- counter + 1
}

results <- do.call(bind_rows, results) %>% as_tibble
results

### confirm n_sim=500
results %>%
  tibble %>%
  group_by(n, J, distn, alt) %>%
  summarise(n())

# add in type I error results 
power_res <- results 
files <- list.files("fig2_updated/results", full.names=T)
results <- vector(length(files), mode = "list")
counter <- 1
for(file in files){
  res_so_far <- readRDS(file)
  results[[counter]] <- res_so_far
  counter <- counter + 1
}

results <- do.call(bind_rows, results) %>% as_tibble
results <- results %>% mutate(alt = 0)
full_res <- rbind(power_res, results)

# add additional results 

files <- list.files("fig3_updated/additional_results", full.names=T)

results <- vector(length(files), mode = "list")
counter <- 1
for(file in files){
  res_so_far <- readRDS(file)
  results[[counter]] <- res_so_far
  counter <- counter + 1
}

results <- do.call(bind_rows, results) %>% as_tibble
results

### confirm n_sim=500
results %>%
  tibble %>%
  group_by(n, J, distn, alt) %>%
  summarise(n())

# add in type I error results 
power_res <- results 
files <- list.files("fig2_updated/additional_results", full.names=T)
results <- vector(length(files), mode = "list")
counter <- 1
for(file in files){
  res_so_far <- readRDS(file)
  results[[counter]] <- res_so_far
  counter <- counter + 1
}

results <- do.call(bind_rows, results) %>% as_tibble
results <- results %>% mutate(alt = 0)
full_res_add <- rbind(power_res, results)

fuller_res <- full_join(full_res, full_res_add)
fuller_res <- fuller_res %>% dplyr::select(-(contains(c("aldex", "ancom", "clr", "deseq"))))

# linda files 
files <- list.files("fig3_updated/linda_results", full.names=T)

linda_results <- vector(length(files), mode = "list")
counter <- 1
for(file in files){
  res_so_far <- readRDS(file)
  linda_results[[counter]] <- res_so_far
  counter <- counter + 1
}

linda_results <- do.call(bind_rows, linda_results) %>% as_tibble
full_results <- full_join(fuller_res, linda_results, by = c("n", "J", "seed", "distn", "alt"))
sum(is.na(full_results$linda_p.x))
sum(is.na(full_results$linda_p.y))
all_results <- full_results %>%
  mutate(linda_est.x = ifelse(is.na(linda_est.x), linda_est.y, linda_est.x),
         linda_p.x = ifelse(is.na(linda_p.x), linda_p.y, linda_p.x)) %>%
  dplyr::select(-linda_est.y, -linda_p.y)
sum(is.na(all_results$linda_p.x))

plot_res <- all_results %>%
  dplyr::select(-contains("est")) %>%
  pivot_longer(6:10, names_to = "method", values_to = "pval") %>%
  mutate(n = factor(n, levels = c(10,50,250))) %>%
  mutate(J = factor(J, levels = c(10,50,250))) %>%
  mutate(distn_J = paste("Y ~ ", distn,"\n",J," Taxa",sep = "")) %>%
  mutate(distn_J = factor(distn_J,
                          levels =
                            c("Y ~ Poisson\n10 Taxa",
                              "Y ~ Poisson\n50 Taxa",
                              "Y ~ Poisson\n250 Taxa",
                              "Y ~ ZINB\n10 Taxa",
                              "Y ~ ZINB\n50 Taxa",
                              "Y ~ ZINB\n250 Taxa")),
         n = paste0("n = ", n)) %>%
  mutate(n = factor(n, levels = 
                      c("n = 10", "n = 50", "n = 250"))) %>%
  mutate(test = ifelse(method == "linda_p.x", "LinDA", 
                       ifelse(method == "locom_p", "LOCOM", 
                              ifelse(method == "maaslin_p", "MaAsLin3",
                                     ifelse(method == "wald_p",
                                            "radEmu robust Wald test",
                                            "radEmu robust score test"))))) %>%
  group_by(n, distn_J, test, alt) %>%
  summarise(power = mean(pval <= 0.05, na.rm = T),
            prop_na = mean(is.na(pval)))

# see which methods control type 1 error rate for each setting
upper_lim <- 0.05 + qnorm(0.975) * sqrt(0.05 * 0.95 / 500)
t1e_res <- plot_res %>%
  filter(alt == 0) %>%
  mutate(valid = ifelse(power <= upper_lim, T, F))
plot_res$valid <- NA
for (row in 1:nrow(plot_res)) {
  n <- plot_res$n[row]
  distn_J <- plot_res$distn_J[row]
  test <- plot_res$test[row]
  plot_res$valid[row] <- t1e_res$valid[t1e_res$n == n & t1e_res$distn_J == distn_J &
                           t1e_res$test == test]
}

# see which methods have fewer than 5% of NA's in each setting
na_res <- plot_res %>%
  group_by(test, n, distn_J) %>%
  summarise(prop_na = mean(prop_na)) %>%
  mutate(too_many_na = ifelse(prop_na > 0.05, TRUE, FALSE))
plot_res$too_many_na <- NA
for (row in 1:nrow(plot_res)) {
  n <- plot_res$n[row]
  distn_J <- plot_res$distn_J[row]
  test <- plot_res$test[row]
  plot_res$too_many_na[row] <- na_res$too_many_na[na_res$n == n & na_res$distn_J == distn_J &
                                                    na_res$test == test]
}

plot_res %>% filter(valid, !too_many_na) %>%
ggplot() +
  geom_line(aes(x = alt, y = power, color = test, linetype = test)) + 
  facet_grid(n ~ distn_J) + 
  #ggtitle("Power simulations") + 
  labs(x = expression(paste(beta[1])),
       y = "Power",
       color = "Test", linetype = "Test") + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  guides(color = guide_legend(position = "bottom", nrow = 2),
         linetype = guide_legend(position = "bottom", nrow = 2)) + 
  #scale_linetype_manual(values=c(6, 5, 4, 4, 1, 2)) + 
  scale_color_manual(values = c("#E69F00", "#CC79A7", 
                       "#661100",
                       "#3446eb",  "#56B4E9")) + 
  scale_linetype_manual(values = c("longdash", "dotdash", "dotted", "solid", "twodash")) + 
  xlim(c(0, 5))
  NULL
ggsave("fig3_updated/power_additional.pdf", height = 8, width = 12)

# estimates
est_res <- all_results %>%
  dplyr::select(-contains("_p")) %>%
  pivot_longer(6:9, names_to = "method", values_to = "est") %>%
  mutate(n = factor(n, levels = c(10,50,250))) %>%
  mutate(J = factor(J, levels = c(10,50,250))) %>%
  mutate(distn_J = paste("Y ~ ", distn,"\n",J," Taxa",sep = "")) %>%
  mutate(distn_J = factor(distn_J,
                          levels =
                            c("Y ~ Poisson\n10 Taxa",
                              "Y ~ Poisson\n50 Taxa",
                              "Y ~ Poisson\n250 Taxa",
                              "Y ~ ZINB\n10 Taxa",
                              "Y ~ ZINB\n50 Taxa",
                              "Y ~ ZINB\n250 Taxa")),
         n = paste0("n = ", n)) %>%
  mutate(n = factor(n, levels = 
                      c("n = 10", "n = 50", "n = 250"))) %>% 
  mutate(test = ifelse(method == "linda_est.x", "LinDA", 
                       ifelse(method == "locom_est", "LOCOM", 
                              ifelse(method == "maaslin_est", "MaAsLin3", "radEmu")))) %>%
  mutate(corrected_est = ifelse(test %in% c("LinDA"),
                                est*log(2), est))
est_res %>% filter(test == "MaAsLin3") %>%
  ggplot(aes(x = n, y = corrected_est)) + 
  geom_hline(aes(yintercept = alt), color = "red") + 
  geom_boxplot() + 
  facet_grid(distn_J~alt) + 
  theme_bw(base_size = 18) + 
  labs(x = "Sample size",
       y = "Estimate") + 
  theme(axis.text.x = element_text(size = 12))
