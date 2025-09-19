library(tidyverse)

#############################################
##### load in weak alternative data
#############################################
files <- list.files("power_weak/", full.names=T)
# files <- list.files("../../research/radEmu/sim_weak_alt_full_info_500_jan_4", full.names=T)

is_rds <- sapply(files,function(x) grepl(".rds",x,fixed = TRUE))

files <- files[is_rds]

results <- vector(length(files),mode = "list")
counter <- 1
for(file in files){
  res_so_far <- readRDS(file)
  get_datas <- file %>% strsplit("n_") %>% extract2(1) %>% strsplit("_") %>% extract2(2)
  # get_datas <- file %>% strsplit("n_") %>% extract2(1) %>% strsplit("_") %>% extract2(3)
  res_so_far$n <- get_datas[1]
  res_so_far$J <- get_datas[3]
  res_so_far$distn <- file %>% strsplit("n_") %>% extract2(1)  %>% strsplit(".rds") %>% extract2(3)
  # res_so_far$distn <- file %>% strsplit("n_") %>% extract2(1)  %>% strsplit(".rds") %>% extract2(4)
  results[[counter]] <- res_so_far
  counter <- counter + 1
}
weak_results <- do.call(rbind,results) %>% as_tibble

#############################################
##### load in strong alternative data
#############################################
files <- list.files("power_strong/", full.names=T)
# files <- list.files("../../research/radEmu/sim_strong_alt_full_info_500_jan_4/", full.names=T)


is_rds <- sapply(files,function(x) grepl(".rds",x,fixed = TRUE))

files <- files[is_rds]

results <- vector(length(files),mode = "list")
counter <- 1
for(file in files){
  res_so_far <- readRDS(file)
  get_datas <- file %>% strsplit("n_") %>% extract2(1) %>% strsplit("_") %>% extract2(2)
  # get_datas <- file %>% strsplit("n_") %>% extract2(1) %>% strsplit("_") %>% extract2(3)
  res_so_far$n <- get_datas[1]
  res_so_far$J <- get_datas[3]
  res_so_far$distn <- file %>% strsplit("n_") %>% extract2(1)  %>% strsplit(".rds") %>% extract2(3)
  # res_so_far$distn <- file %>% strsplit("n_") %>% extract2(1)  %>% strsplit(".rds") %>% extract2(4)
  results[[counter]] <- res_so_far
  counter <- counter + 1
}

strong_results <- do.call(rbind,results) %>% as_tibble

results <- bind_rows(mutate(weak_results, alt = "weak"),
                     mutate(strong_results, alt = "strong"))

pd <- position_dodge(0.5)

weak_fig <- weak_results %>%
  pivot_longer(c(score_p_null_info,wald_p)) %>%
  mutate(name = ifelse(name == "score_p_null_info", "Robust Score Test",
                       ifelse(name == "wald_p", "Robust Wald Test", "um what?!"))) %>%
  group_by(n,J,distn,name) %>%
  summarize(emp_type1 = mean(value<=0.05,na.rm = TRUE),
            lower = binom.test(x = sum(value <=0.05,na.rm = TRUE),
                               n = sum(!is.na(value)))$conf.int[1],
            upper = binom.test(x = sum(value <=0.05,na.rm = TRUE),
                               n = sum(!is.na(value)))$conf.int[2],
            nsims = length(value),

            nfailed = sum(is.na(value))) %>%
  mutate(n = factor(n, levels = c(10,50,250))) %>%
  mutate(J = paste(J, "taxa")) %>%
  mutate(J = factor(J, levels = c("10 taxa","50 taxa","250 taxa"))) %>%
  rename("Test" = name) %>%
  mutate(distn = paste("Y ~", distn)) %>%
  ggplot(aes(x = n,
             color = Test)) +
  geom_point(aes(y = emp_type1),
             position = pd,
             size= 0.5) +
  geom_errorbar(aes(ymin = lower, ymax= upper),
                width = 0.25,
                position = pd) +
  facet_grid(distn~J) +
  geom_abline(aes(slope = 0, intercept = 0.05),linetype = "dotted") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("Power under weak alternative") +
  xlab("") +
  ylim(0,1) +
  xlab("Sample size") +
  scale_color_manual(values=c("#3446eb",  "#56B4E9")) + # "#208f1a",
  NULL

strong_fig <- strong_results %>%
  pivot_longer(c(score_p_null_info,wald_p)) %>%
  mutate(name = ifelse(name == "score_p_null_info", "Robust Score Test",
                       ifelse(name == "wald_p", "Robust Wald Test", "um what?!"))) %>%
  group_by(n,J,distn,name) %>%
  summarize(emp_type1 = mean(value<=0.05,na.rm = TRUE),
            lower = binom.test(x = sum(value <=0.05,na.rm = TRUE),
                               n = sum(!is.na(value)))$conf.int[1],
            upper = binom.test(x = sum(value <=0.05,na.rm = TRUE),
                               n = sum(!is.na(value)))$conf.int[2],
            nsims = length(value),

            nfailed = sum(is.na(value))) %>%
  mutate(n = factor(n, levels = c(10,50,250))) %>%
  mutate(J = paste(J, "taxa")) %>%
  mutate(J = factor(J, levels = c("10 taxa","50 taxa","250 taxa"))) %>%
  rename("Test" = name) %>%
  mutate(distn = paste("Y ~", distn)) %>%
  ggplot(aes(x = n,
             color = Test)) +
  geom_point(aes(y = emp_type1),
             position = pd,
             size= 0.5) +
  geom_errorbar(aes(ymin = lower, ymax= upper),
                width = 0.25,
                position = pd) +
  facet_grid(distn~J) +
  geom_abline(aes(slope = 0, intercept = 0.05),linetype = "dotted") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylim(0,1) +
  ylab("Power under strong alternative") +
  xlab("Sample size") +
  scale_color_manual(values=c("#3446eb",  "#56B4E9")) + # "#208f1a",
  NULL

ggpubr::ggarrange(weak_fig, strong_fig, ncol=2, common.legend=T, legend="bottom")
# ggsave("images/power-wide.pdf", width=8.5, height = 4.0)

