library(tidyverse)

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
results_other <- do.call(bind_rows, results_other) %>% as_tibble %>%
  dplyr::select(-(contains(c("aldex", "ancom", "clr", "deseq"))))
results <- full_join(results, results_other)

# linda files 
files <- list.files("fig3_updated/linda_results", full.names=T)

linda_results <- vector(length(files), mode = "list")
counter <- 1
for(file in files){
  res_so_far <- readRDS(file)
  linda_results[[counter]] <- res_so_far
  counter <- counter + 1
}

linda_results <- do.call(bind_rows, linda_results) %>% as_tibble %>%
  filter(alt == 0) %>%
  dplyr::select(-alt)
full_results <- full_join(results, linda_results, by = c("n", "J", "seed", "distn"))
sum(is.na(full_results$linda_p.x))
sum(is.na(full_results$linda_p.y))
all_results <- full_results %>%
  mutate(linda_est.x = ifelse(is.na(linda_est.x), linda_est.y, linda_est.x),
         linda_p.x = ifelse(is.na(linda_p.x), linda_p.y, linda_p.x)) %>%
  dplyr::select(-linda_est.y, -linda_p.y)
sum(is.na(all_results$linda_p.x))

### confirm n_sim=500 (yes)
all_results %>%
  tibble %>%
  group_by(n, J, distn) %>%
  summarise(n())

plot_res <- all_results %>%
  dplyr::select(-contains("est")) %>%
  pivot_longer(5:9, names_to = "method", values_to = "pval") %>%
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
                                                   "radEmu robust score test")))))

# see which methods have fewer than 5% of NA's in each setting
na_res <- plot_res %>%
  group_by(test, n, distn_J) %>%
  summarise(prop_na = mean(is.na(pval))) %>%
  mutate(too_many_na = ifelse(prop_na > 0.05, TRUE, FALSE))
plot_res$too_many_na <- NA
for (row in 1:nrow(plot_res)) {
  n <- plot_res$n[row]
  distn_J <- plot_res$distn_J[row]
  test <- plot_res$test[row]
  plot_res$too_many_na[row] <- na_res$too_many_na[na_res$n == n & na_res$distn_J == distn_J &
                                                    na_res$test == test]
}
plot_res %>% filter(!too_many_na) %>% 
ggplot() +
  geom_qq(distribution = stats::qunif,
          aes(sample = pval,
              color = test, linetype = test
          ), geom="line",
          linewidth = 0.5) +
  #geom_abline(aes(intercept = 0, slope = 1),linetype= "dotted") +
  geom_abline(aes(intercept = 0, slope = 1)) +
  facet_grid(n~distn_J) +
  labs(color = "Test", linetype = "Test") +
  xlab("Theoretical p-value quantiles") +
  ylab("Empirical p-value quantiles") +
  theme_bw(base_size = 20) +
  theme(panel.grid.minor = element_blank()) +
  #scale_color_manual(values=c("#3446eb",  "#56B4E9")) + # "#208f1a",
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12, angle = 60, hjust= 1),
        axis.text.y = element_text(size = 12)) +
  scale_x_sqrt(breaks=c(0,0.01,0.05,0.25,0.5,1)) +
  scale_y_sqrt(breaks = c(0,0.01,0.05,0.25,0.5,1)) +
  coord_equal() +
  scale_color_manual(values = c("#E69F00", "#CC79A7", 
                                "#661100",
                                "#3446eb",  "#56B4E9")) + 
  guides(color = guide_legend(position = "bottom", nrow = 2),
         linetype = guide_legend(position = "bottom", nrow = 2)) + 
  scale_linetype_manual(values = c("longdash", "dotdash", "dotted", "solid", "twodash"))
  NULL
ggsave("fig2_updated/t1e_additional.pdf", height = 8, width = 12)

plot_res %>% group_by(n, distn_J, test) %>%
 # filter(test == "radEmu robust score test") %>%
  summarise(mean(pval <= 0.05, na.rm = T)) %>%
  print(n = 90)

# estimates
plot_res <- all_results %>%
  dplyr::select(-contains("_p")) %>%
  pivot_longer(5:8, names_to = "method", values_to = "est") %>%
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
ggplot(plot_res, aes(x = n, y = corrected_est)) + 
  geom_hline(aes(yintercept = 0), color = "red") + 
  geom_boxplot() + 
  facet_grid(distn_J~test) + 
  theme_bw(base_size = 18) + 
  labs(x = "Sample size",
       y = "Estimate") + 
  theme(axis.text.x = element_text(size = 12))
ggplot(plot_res, aes(x = n, y = corrected_est)) + 
  geom_hline(aes(yintercept = 0), color = "red") + 
  geom_boxplot() + 
  facet_grid(distn_J~test) + 
  theme_bw(base_size = 18) + 
  labs(x = "Sample size",
       y = "Estimate") + 
  theme(axis.text.x = element_text(size = 12)) + 
  ylim(c(-5, 5))
ggsave("fig2_updated/bias_additional.pdf", height = 8, width = 12)

ggplot(plot_res, aes(x = n, y = corrected_est)) + 
  geom_hline(aes(yintercept = 0), color = "red") + 
  geom_boxplot() + 
  facet_grid(distn_J~test) + 
  theme_bw(base_size = 18) + 
  labs(x = "Sample size",
       y = "Estimate") + 
  theme(axis.text.x = element_text(size = 12)) + 
  ylim(c(-1, 1))
