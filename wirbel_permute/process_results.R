library(tidyverse)
library(qvalue)

# aldex, ancom, clr t-test
files <- list.files("wirbel_permute/results/", full.names = T)
is_other <- sapply(files, function(x) grepl("other", x, fixed = TRUE))
files <- files[is_other]
df_list <- vector(mode = "list", length = length(files))

for (ind in 1:length(files)) {
  df_list[[ind]] <- readRDS(files[ind])
}
other_methods <- df_list

# deseq 
files <- list.files("wirbel_permute/results/", full.names = T)
is_deseq <- sapply(files, function(x) grepl("deseq", x, fixed = TRUE))
files <- files[is_deseq]
df_list <- vector(mode = "list", length = length(files))

for (ind in 1:length(files)) {
  df_list[[ind]] <- readRDS(files[ind])
}
# check on number na's for each row in each permutation
deseq_na <- unlist(lapply(df_list, function(x) {sum(is.na(x[, 5]))}))

# check on significant qvals
# aldex
qval <- lapply(other_methods, function(x) {qvalue::qvalue(x[, 5], pi0 = 1)$qvalues})
aldex_sum_small <- unlist(lapply(qval, function(x) {sum(x <= 0.05, na.rm = T)}))
# ancom
qval <- lapply(other_methods, function(x) {cbind(x, qvalue::qvalue(x[, 7])$qvalues)})
ancom_sum_small <- unlist(lapply(qval, function(x) {sum(x[, 13] <= 0.05, na.rm = T)}))
ancom_sum_small_ss <- unlist(lapply(qval, function(x) {sum(x[, 13] <= 0.05 & x[, 8] == T, na.rm = T)}))
# clr
qval <- lapply(other_methods, function(x) {qvalue::qvalue(x[, 10])$qvalues})
clr_sum_small <- unlist(lapply(qval, function(x) {sum(x <= 0.05, na.rm = T)}))
# deseq 
qval <- lapply(df_list, function(x) {qvalue::qvalue(x[, 5])$qvalues})
deseq_sum_small <- unlist(lapply(qval, function(x) {sum(x <= 0.05, na.rm = T)}))

# score tests, see what has run
Y <- readRDS("wirbel_permute/data/Y.rds")
files <- list.files("wirbel_permute/results/", full.names = T)
is_score <- sapply(files, function(x) grepl("score_perm", x, fixed = TRUE))
files <- files[is_score]
score_list <- vector(mode = "list", length = 20)
for (ind in 1:length(files)) {
  perm <- parse_number(files[ind])
  score_list[[perm]] <- rbind(score_list[[perm]], readRDS(files[ind]))
}
ran <- lapply(score_list, function(x) {x[, 2]})
list_remain <- lapply(ran, function(x) 
  {(1:ncol(Y))[!((1:ncol(Y)) %in% x)]})
num_remain <- unlist(lapply(list_remain, length))
sum(num_remain)
df_remain <- data.frame(tests = unlist(list_remain),
                        perm = unlist(sapply(1:20, function(x) {rep(x, num_remain[x])})))
df_remain$tax <- colnames(Y)[df_remain$tests]
saveRDS(df_remain, "wirbel_permute/results/df_remain.rds")

# add additional re-run score tests 
files <- list.files("wirbel_permute/results/", full.names = T)
is_rerun <- sapply(files, function(x) grepl("rerun", x, fixed = TRUE))
files <- files[is_rerun]
df_list <- vector(mode = "list", length = length(files))
for (ind in 1:length(files)) {
  df_list[[ind]] <- readRDS(files[ind])
}
rerun_res <- do.call(bind_rows, df_list) %>% as_tibble
df_rerun <- full_join(df_remain, rerun_res)
df_run_again <- df_rerun %>% filter(is.na(score_p))
saveRDS(df_run_again, "wirbel_permute/results/df_run_again.rds")

# add more score tests 
files <- list.files("wirbel_permute/results/", full.names = T)
is_rerun <- sapply(files, function(x) grepl("run_again_score", x, fixed = TRUE))
files <- files[is_rerun]
df_list <- vector(mode = "list", length = length(files))
for (ind in 1:length(files)) {
  df_list[[ind]] <- readRDS(files[ind])
}
run_again_res <- do.call(bind_rows, df_list) %>% as_tibble
all_rerun <- full_join(run_again_res, df_rerun, by = c("tests", "perm", "tax")) %>%
  mutate(score_p = ifelse(!is.na(score_p.x), score_p.x,
                          ifelse(!is.na(score_p.y), score_p.y, NA))) %>%
  dplyr::select(-score_p.x, -score_p.y)

# add rerun score tests to original score test results 
for (ind in 1:nrow(all_rerun)) {
  perm <- all_rerun$perm[ind]
  score_list[[perm]] <- rbind(score_list[[perm]], all_rerun[ind, ])
}
lapply(score_list, dim)

# compare q-values, including 23 NA values across 20 permutations 
qval <- lapply(score_list, function(x) {qvalue::qvalue(x[, 4])$qvalues})
score_sum_small <- unlist(lapply(qval, function(x) {sum(x <= 0.05, na.rm = T)}))
score_list_0 <- lapply(score_list, function(x) {
  x %>% mutate(score_p = ifelse(is.na(score_p), 0, score_p))
})
qval <- lapply(score_list_0, function(x) {qvalue::qvalue(x[, 4])$qvalues})
score_sum_small_0 <- unlist(lapply(qval, function(x) {sum(x <= 0.05, na.rm = T)}))
unlist(lapply(score_list, function(x) {sum(is.na(x[, 4]))}))
