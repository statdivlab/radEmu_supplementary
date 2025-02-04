library(tidyverse)

files <- list.files("fig2_updated/results", full.names=T)

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
  group_by(n, J, distn) %>%
  summarise(n())

plot_res <- results %>%
  dplyr::select(-contains("est")) %>%
  pivot_longer(5:10, names_to = "method", values_to = "pval") %>%
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
                                                   "radEmu robust score test")))))) 
ggplot(plot_res) +
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
                                "#661100", "#009E73", 
                                "#3446eb",  "#56B4E9")) + 
  scale_linetype_manual(values = c("dashed", "longdash", "dotdash", "dotted", "solid", "twodash"))
  NULL
ggsave("fig2_updated/t1e.pdf", height = 8, width = 12)

plot_res %>% group_by(n, distn_J, test) %>%
 # filter(test == "radEmu robust score test") %>%
  summarise(mean(pval <= 0.05, na.rm = T)) %>%
  print(n = 108)

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
                              "Y ~ ZINB\n250 Taxa")),
         n = paste0("n = ", n)) %>%
  mutate(n = factor(n, levels = 
                      c("n = 10", "n = 50", "n = 250"))) %>% 
  mutate(test = ifelse(method == "aldex_est", "ALDEx2", 
                       ifelse(method == "ancom_est", "ANCOM-BC2", 
                              ifelse(method == "clr_est", "CLR t-test",
                                     ifelse(method == "deseq_est", "DESeq2", "radEmu"))))) %>%
  mutate(corrected_est = ifelse(test %in% c("ALDEx2", "DESeq2"),
                                est*log(2), est))
ggplot(plot_res, aes(x = n, y = corrected_est)) + 
  geom_hline(aes(yintercept = 0), color = "red") + 
  geom_boxplot() + 
  facet_grid(distn_J~test) + 
  theme_bw(base_size = 18) + 
  labs(x = "Sample size",
       y = "Estimate") + 
  theme(axis.text.x = element_text(size = 12))
ggsave("fig2_updated/bias.pdf", height = 8, width = 12)

summary(plot_res %>% filter(test == "ALDEx2") %>% pull(pval))
