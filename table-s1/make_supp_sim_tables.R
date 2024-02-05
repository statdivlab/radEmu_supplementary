#### Scripts to reproduce the data analyses in the `radEmu` paper
####
#### David Clausen and Amy Willis
#### Dec 2023 -- Jan 2024

#### Reproduce the tables in the supplement

files <- list.files("t1e_sims_tmp/", full.names=T)
# files <- list.files("../../research/radEmu/sim_results_full_info_10k_dec_19/", full.names=T)

is_rds <- sapply(files,function(x) grepl(".rds",x,fixed = TRUE))

files <- files[is_rds]

# results <- vector(length(files),mode = "list")
# counter <- 1
# for(file in files){
#   res_so_far <- readRDS(file)
#   res_so_far$n <- strsplit(file,"_")[[1]][2]
#   res_so_far$J <- strsplit(file,"_")[[1]][4]
#   res_so_far$distn <- strsplit( strsplit(file,"_")[[1]][6],".",fixed = TRUE)[[1]][1]
#   results[[counter]] <- res_so_far
#   counter <- counter + 1
# }

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

results <- do.call(rbind,results)

##Null type 1 for score
results %>%
  mutate(n = sapply(n, function(x) paste("n = ",x,sep = "", collapse = ""))) %>%
  mutate(n = factor(n, levels = c("n = 10","n = 50","n = 250")),
         J = factor(J, levels = c(10,50,250)),
         distn = factor(distn, levels = c("Poisson","ZINB"))) %>%
  group_by(n,distn,J) %>%
  summarize(emp_05_score = mean(score_p_null_info<=0.05, na.rm = TRUE)) %>%
  pivot_wider(values_from = emp_05_score,names_from= n) %>%
  as.data.frame() %>%
  xtable::xtable(digits = 3) %>%
  (function(x){
    colnames(x)[1] <- "Distn.";
    return(x)
  }) %>%
  xtable::xtable(digits = 3)

##null type 1 for wald
results %>%
  mutate(n = sapply(n, function(x) paste("n = ",x,sep = "", collapse = ""))) %>%
  mutate(n = factor(n, levels = c("n = 10","n = 50","n = 250")),
         J = factor(J, levels = c(10,50,250)),
         distn = factor(distn, levels = c("Poisson","ZINB"))) %>%
  group_by(n,distn,J) %>%
  summarize(emp_05_wald = mean(wald_p<=0.05, na.rm = TRUE)) %>%
  pivot_wider(values_from = emp_05_wald,names_from= n) %>%
  as.data.frame() %>%
  (function(x){
    colnames(x)[1] <- "Distn.";
    return(x)
  }) %>%
  xtable::xtable(digits = 3)

### Weak alternative summaries

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
weak_results

##weak alt type 1 for score
weak_results %>%
  mutate(n = sapply(n, function(x) paste("n = ",x,sep = "", collapse = ""))) %>%
  mutate(n = factor(n, levels = c("n = 10","n = 50","n = 250")),
         J = factor(J, levels = c(10,50,250)),
         distn = factor(distn, levels = c("Poisson","ZINB"))) %>%
  group_by(n,distn,J) %>%
  summarize(emp_05_score = mean(score_p_null_info<=0.05, na.rm = TRUE)) %>%
  pivot_wider(values_from = emp_05_score,names_from= n) %>%
  as.data.frame() %>%
  xtable::xtable(digits = 3) %>%
  (function(x){
    colnames(x)[1] <- "Distn.";
    return(x)
  }) %>%
  xtable::xtable(digits = 3)

##weak alt type 1 for wald
weak_results %>%
  mutate(n = sapply(n, function(x) paste("n = ",x,sep = "", collapse = ""))) %>%
  mutate(n = factor(n, levels = c("n = 10","n = 50","n = 250")),
         J = factor(J, levels = c(10,50,250)),
         distn = factor(distn, levels = c("Poisson","ZINB"))) %>%
  group_by(n,distn,J) %>%
  summarize(emp_05_wald = mean(wald_p<=0.05, na.rm = TRUE)) %>%
  pivot_wider(values_from = emp_05_wald,names_from= n) %>%
  as.data.frame() %>%
  (function(x){
    colnames(x)[1] <- "Distn.";
    return(x)
  }) %>%
  xtable::xtable(digits = 3)

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

##strong alt type 1 for score
strong_results %>%
  mutate(n = sapply(n, function(x) paste("n = ",x,sep = "", collapse = ""))) %>%
  mutate(n = factor(n, levels = c("n = 10","n = 50","n = 250")),
         J = factor(J, levels = c(10,50,250)),
         distn = factor(distn, levels = c("Poisson","ZINB"))) %>%
  group_by(n,distn,J) %>%
  summarize(emp_05_score = mean(score_p_null_info<=0.05, na.rm = TRUE)) %>%
  pivot_wider(values_from = emp_05_score,names_from= n) %>%
  as.data.frame() %>%
  xtable::xtable(digits = 3) %>%
  (function(x){
    colnames(x)[1] <- "Distn.";
    return(x)
  }) %>%
  xtable::xtable(digits = 3)

##strong alt type 1 for wald
strong_results %>%
  mutate(n = sapply(n, function(x) paste("n = ",x,sep = "", collapse = ""))) %>%
  mutate(n = factor(n, levels = c("n = 10","n = 50","n = 250")),
         J = factor(J, levels = c(10,50,250)),
         distn = factor(distn, levels = c("Poisson","ZINB"))) %>%
  group_by(n,distn,J) %>%
  summarize(emp_05_wald = mean(wald_p<=0.05, na.rm = TRUE)) %>%
  pivot_wider(values_from = emp_05_wald,names_from= n) %>%
  as.data.frame() %>%
  (function(x){
    colnames(x)[1] <- "Distn.";
    return(x)
  }) %>%
  xtable::xtable(digits = 3)
