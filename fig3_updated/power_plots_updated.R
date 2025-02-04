library(tidyverse)

files <- list.files("fig3_updated/results", full.names=T)

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

plot_res <- full_res %>%
  dplyr::select(-contains("est")) %>%
  pivot_longer(6:11, names_to = "method", values_to = "pval") %>%
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
  mutate(test = ifelse(method == "aldex_p", "ALDEx2", 
                       ifelse(method == "ancom_p", "ANCOM-BC2", 
                              ifelse(method == "clr_p", "CLR t-test",
                                     ifelse(method == "deseq_p", "DESeq2",
                                            ifelse(method == "wald_p",
                                                   "radEmu robust Wald test",
                                                   "radEmu robust score test")))))) %>%
  group_by(n, distn_J, test, alt) %>%
  summarise(power = mean(pval <= 0.05, na.rm = T))

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

ggplot(plot_res %>% filter(valid)) +
  geom_line(aes(x = alt, y = power, color = test, linetype = test)) + 
  facet_grid(n ~ distn_J) + 
  ggtitle("Power simulations") + 
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
  scale_linetype_manual(values = c("dashed", "longdash", "dotdash", "dotted", "solid", "twodash")) + 
  xlim(c(0,5)) + 
  scale_color_manual(values = c("#E69F00", "#CC79A7", 
                                "#661100", "#009E73", 
                                "#3446eb",  "#56B4E9")) + 
  NULL
ggsave("fig3_updated/power.pdf", height = 8, width = 12)

