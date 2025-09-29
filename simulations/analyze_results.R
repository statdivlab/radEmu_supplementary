library(tidyverse)

files <- list.files("simulations/results", full.names=T)

results <- vector(length(files), mode = "list")
counter <- 1
for(file in files){
  res_so_far <- readRDS(file)
  results[[counter]] <- res_so_far
  counter <- counter + 1
}

results <- do.call(bind_rows, results) %>% as_tibble
results

# see where we're missing trials (should be 504 * 3 = 1512 per setting)
results %>%
  tibble %>%
  group_by(n, J, distn, abs(alt)) %>%
  summarise(n()) %>%
  print(n = 126)

# see which methods failed
results %>% 
  tibble %>%
  group_by(n, J, distn, abs(alt)) %>%
  summarise(ald_na = mean(is.na(aldex_p)),
            anc_na = mean(is.na(ancom_p)),
            clr_na = mean(is.na(clr_p)),
            des_na = mean(is.na(deseq_p)),
            lind_na = mean(is.na(linda_p)),
            loc_na = mean(is.na(locom_p)),
            maas_na = mean(is.na(maaslin_p)),
            wald_na = mean(is.na(wald_p)),
            score_na = mean(is.na(score_p))) %>%
  print(n = 126)

# type I error plots 
plot_res <- results %>%
  filter(alt == 0) %>%
  dplyr::select(-contains("est"), -contains("time"), -delta_level, -alt) %>%
  pivot_longer(5:13, names_to = "method", values_to = "pval") %>%
  mutate(n = factor(n, levels = c(10,50,250))) %>%
  mutate(J = factor(J, levels = c(10,50,250))) %>%
  mutate(distn_J = paste("Y ~ ", distn,"\n",J," categories",sep = "")) %>%
  mutate(distn_J = factor(distn_J,
                          levels =
                            c("Y ~ Poisson\n10 categories",
                              "Y ~ Poisson\n50 categories",
                              "Y ~ Poisson\n250 categories",
                              "Y ~ ZINB\n10 categories",
                              "Y ~ ZINB\n50 categories",
                              "Y ~ ZINB\n250 categories")),
         n = paste0("n = ", n)) %>%
  mutate(n = factor(n, levels = 
                      c("n = 10", "n = 50", "n = 250"))) %>% 
  mutate(test = ifelse(method == "aldex_p", "ALDEx2", 
                       ifelse(method == "ancom_p", "ANCOM-BC2", 
                              ifelse(method == "clr_p", "CLR linear model",
                                     ifelse(method == "deseq_p", "DESeq2",
                                            ifelse(method == "linda_p", "LinDA", 
                                                   ifelse(method == "locom_p", "LOCOM", 
                                                          ifelse(method == "maaslin_p", "MaAsLin3",
                                                                 ifelse(method == "wald_p",
                                                                        "radEmu robust Wald test",
                                                                        "radEmu robust score test")))))))))

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

# original methods
plot_res %>% filter(!too_many_na) %>%
  filter(test %in% c("ALDEx2", "ANCOM-BC2", "CLR linear model", "DESeq2", 
                     "radEmu robust Wald test", "radEmu robust score test")) %>%
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
  theme_bw(base_size = 19) +
  theme(panel.grid.minor = element_blank()) +
  #scale_color_manual(values=c("#3446eb",  "#56B4E9")) + # "#208f1a",
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12, angle = 60, hjust= 1),
        axis.text.y = element_text(size = 12)) +
  scale_x_sqrt(breaks=c(0,0.01,0.05,0.25,0.5,1)) +
  scale_y_sqrt(breaks = c(0,0.01,0.05,0.25,0.5,1)) +
  coord_equal() +
  scale_linetype_manual(values = c("dashed", "longdash", "dotdash", "dotted", "solid", "twodash")) + 
  scale_color_manual(values = c("#E69F00", "#CC79A7", 
                                "#661100", "#009E73", 
                                "#3446eb",  "#56B4E9"))
NULL
ggsave("simulations/revised_t1e.pdf", height = 8, width = 12)

# additional methods
plot_res %>% filter(!too_many_na) %>%
  filter(test %in% c("LinDA", "LOCOM", "MaAsLin3", 
                     "radEmu robust Wald test", "radEmu robust score test")) %>%
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
  theme_bw(base_size = 19) +
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
  scale_linetype_manual(values = c("longdash", "dotdash", "dotted", "solid", "twodash")) + 
  NULL
ggsave("simulations/revised_t1e_additional.pdf", height = 8, width = 12)

lower <- 0.05 - qnorm(0.975) * sqrt(.05 * .95 / 500)
upper <- 0.05 + qnorm(0.975) * sqrt(.05 * .95 / 500)

# conservative 
plot_res %>% group_by(n, distn_J, test) %>%
  summarise(t1e = mean(pval <= 0.05, na.rm = T)) %>%
  filter(t1e < lower) %>%
  arrange(test) %>%
  print(n = 150)

# anticonservative
plot_res %>% group_by(n, distn_J, test) %>%
  summarise(t1e = mean(pval <= 0.05, na.rm = T)) %>%
  filter(t1e > upper) %>%
  arrange(test) %>%
  print(n = 150)

# power plots
plot_res <- results %>%
  dplyr::select(-contains("est"), -contains("time"), -delta_level) %>%
  pivot_longer(6:14, names_to = "method", values_to = "pval") %>%
  mutate(n = factor(n, levels = c(10,50,250))) %>%
  mutate(J = factor(J, levels = c(10,50,250))) %>%
  mutate(distn_J = paste("Y ~ ", distn,"\n",J," categories",sep = "")) %>%
  mutate(distn_J = factor(distn_J,
                          levels =
                            c("Y ~ Poisson\n10 categories",
                              "Y ~ Poisson\n50 categories",
                              "Y ~ Poisson\n250 categories",
                              "Y ~ ZINB\n10 categories",
                              "Y ~ ZINB\n50 categories",
                              "Y ~ ZINB\n250 categories")),
         n = paste0("n = ", n)) %>%
  mutate(n = factor(n, levels = 
                      c("n = 10", "n = 50", "n = 250"))) %>% 
  mutate(test = ifelse(method == "aldex_p", "ALDEx2", 
                       ifelse(method == "ancom_p", "ANCOM-BC2", 
                              ifelse(method == "clr_p", "CLR linear model",
                                     ifelse(method == "deseq_p", "DESeq2",
                                            ifelse(method == "linda_p", "LinDA", 
                                                   ifelse(method == "locom_p", "LOCOM", 
                                                          ifelse(method == "maaslin_p", "MaAsLin3",
                                                                 ifelse(method == "wald_p",
                                                                        "radEmu robust Wald test",
                                                                        "radEmu robust score test"))))))))) %>%
  mutate(alt = abs(alt)) %>%
  group_by(n, distn_J, test, alt) %>%
  summarise(power = mean(pval <= 0.05, na.rm = T),
            prop_na = mean(is.na(pval)))

# see which methods control type 1 error rate for each setting
t1e_res <- plot_res %>%
  filter(alt == 0) %>%
  mutate(valid = ifelse(power <= upper, T, F))
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

# original methods
plot_res %>% filter(valid, !too_many_na) %>% 
  filter(test %in% c("ALDEx2", "ANCOM-BC2", "CLR linear model", "DESeq2", 
                     "radEmu robust Wald test", "radEmu robust score test")) %>%
  ggplot() +
  geom_line(aes(x = alt, y = power, color = test, linetype = test)) + 
  facet_grid(n ~ distn_J) + 
  #ggtitle("Power simulations") + 
  labs(x = expression(paste(beta[1])),
       y = "Power",
       color = "Test", linetype = "Test") + 
  theme_bw(base_size = 19) + 
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_linetype_manual(values = c("dashed", "longdash", "dotdash", "dotted", "solid", "twodash")) + 
  scale_color_manual(values = c("#E69F00", "#CC79A7", 
                                "#661100", "#009E73", 
                                "#3446eb",  "#56B4E9")) + 
  theme(legend.position = "bottom") + 
  NULL
ggsave("simulations/revised_power.pdf", height = 8, width = 12)

# additional methods
plot_res %>% filter(valid, !too_many_na) %>% 
  filter(test %in% c("LinDA", "LOCOM", "MaAsLin3", 
                     "radEmu robust Wald test", "radEmu robust score test")) %>%
  ggplot() +
  geom_line(aes(x = alt, y = power, color = test, linetype = test)) + 
  facet_grid(n ~ distn_J) + 
  #ggtitle("Power simulations") + 
  labs(x = expression(paste(beta[1])),
       y = "Power",
       color = "Test", linetype = "Test") + 
  theme_bw(base_size = 19) + 
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_color_manual(values = c("#E69F00", #"#CC79A7", 
                                "#661100",
                                "#3446eb",  "#56B4E9")) + 
  scale_linetype_manual(values = c("longdash", #"dotdash", 
                                   "dotted", "solid", "twodash")) + 
  theme(legend.position = "bottom") + 
  NULL
ggsave("simulations/revised_power_additional.pdf", height = 8, width = 12)


