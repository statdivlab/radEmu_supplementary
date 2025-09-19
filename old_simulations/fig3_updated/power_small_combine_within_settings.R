library(tidyverse)

files <- list.files("radEmu/results/power_sims/", full.names = T)
is_small <- sapply(files, function(x) grepl("small", x, fixed = TRUE))
files <- files[is_small]
df_list <- vector(mode = "list", length = 18)

for (file in files) {
  res <- readRDS(file)
  setting <- readr::parse_number(stringr::str_remove(file, "radEmu/results/power_sims/small_setting"))
  if (is.null(df_list[[setting]])) {
    df_list[[setting]] <- res
  } else {
    df_list[[setting]] <- rbind(df_list[[setting]], res)
  }
}

for (s in 1:18) {
  saveRDS(df_list[[s]], paste0("radEmu/results/power_sims/small_full_setting", s, ".RDS"))
}