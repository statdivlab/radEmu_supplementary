library(tidyverse)

files <- list.files("t1e_sims_tmp/", full.names=T)
# files <- list.files("../../research/radEmu/sim_results_full_info_10k_dec_19/", full.names=T)

is_rds <- sapply(files,function(x) grepl(".rds",x,fixed = TRUE))

files <- files[is_rds]

results <- vector(length(files),mode = "list")
counter <- 1
for(file in files){
  res_so_far <- readRDS(file)
  get_datas <- file %>% strsplit("n_") %>% extract2(1) %>% strsplit("_") %>% extract2(2)
  res_so_far$n <- get_datas[1]
  res_so_far$J <- get_datas[3]
  res_so_far$distn <- file %>% strsplit("n_") %>% extract2(1)  %>% strsplit(".rds") %>% extract2(3)
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

results %>%
  pivot_longer(c(score_p_null_info,wald_p)) %>%
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
  mutate(test = sapply(name,
                       function(x)
                         ifelse(x == "score_p_null_info","Robust Score Test",
                                "Robust Wald Test"))) %>%
  ggplot() +
  geom_qq(distribution = stats::qunif,
          aes(sample = value,
              color = test
          ), geom="line",
          linewidth = 0.5) +
  geom_abline(aes(intercept = 0, slope = 1),linetype= "dotted") +
  facet_grid(n~distn_J) +
  labs(color = "Test") +
  xlab("Theoretical p-value quantiles") +
  ylab("Empirical p-value quantiles") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_color_manual(values=c("#3446eb",  "#56B4E9")) + # "#208f1a",
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 60, size=6,hjust= 1),
        axis.text.y = element_text(size = 6)) +
  scale_x_sqrt(breaks=c(0,0.01,0.05,0.25,0.5,1)) +
  scale_y_sqrt(breaks = c(0,0.01,0.05,0.25,0.5,1)) +
  coord_equal() +
  NULL

# ggsave("images/null_10000.pdf", width=7, height = 4)
