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

### confirm n_sim=10000 (yes)
results %>%
  tibble %>%
  group_by(n, J, distn, alt) %>%
  summarise(n())

plot_res <- results %>%
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
                              "Y ~ ZINB\n250 Taxa"))) %>%
  mutate(test = ifelse(method == "aldex_p", "ALDEx2", 
                       ifelse(method == "ancom_p", "ANCOM-BC2", 
                              ifelse(method == "clr_p", "CLR t-test",
                                     ifelse(method == "deseq_p", "DESeq2",
                                            ifelse(method == "wald_p",
                                                   "radEmu robust Wald test",
                                                   "radEmu robust score test")))))) %>%
  group_by(n, distn_J, test, alt) %>%
  summarise(power = mean(pval <= 0.05, na.rm = T))

ggplot(plot_res) +
  geom_line(aes(x = alt, y = power, color = test)) + 
  facet_grid(distn_J ~ n) + 
  ggtitle("Power simulations") + 
  labs(x = expression(paste(beta[1])),
       y = "Power",
       linetype = "Sample size",
       color = "Test") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(position = "bottom", nrow = 3),
         linetype = "none") + 
  #scale_linetype_manual(values=c(2, 2, 2, 4, 2)) + 
  xlim(c(0,5)) + 
  coord_equal() +
  NULL
ggsave("fig3_updated/power.pdf", height = 8, width = 12)

# estimates
plot_res <- results %>%
  dplyr::select(-contains("_p")) %>%
  pivot_longer(5:9, names_to = "method", values_to = "est") %>%
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
                              "Y ~ ZINB\n250 Taxa"))) %>%
  mutate(test = ifelse(method == "aldex_est", "ALDEx2", 
                       ifelse(method == "ancom_est", "ANCOM-BC2", 
                              ifelse(method == "clr_est", "CLR t-test",
                                     ifelse(method == "deseq_est", "DESeq2", "radEmu"))))) 
ggplot(plot_res, aes(x = n, y = est)) + 
  geom_hline(aes(yintercept = 0), color = "red") + 
  geom_boxplot() + 
  facet_grid(distn_J~test) + 
  theme_bw() + 
  labs(x = "Sample size",
       y = "Estimate")
ggsave("fig2_updated/bias.pdf", height = 8, width = 12)
