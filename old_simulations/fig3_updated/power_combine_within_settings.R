library(tidyverse)

files <- list.files("radEmu/results/power_sims/", full.names = T)
df_list <- vector(mode = "list", length = 90)

for (file in files) {
  res <- readRDS(file)
  setting <- readr::parse_number(stringr::str_remove(file, "radEmu/results/power_sims/setting"))
  if (is.null(df_list[[setting]])) {
    df_list[[setting]] <- res
  } else {
    df_list[[setting]] <- rbind(df_list[[setting]], res)
  }
}

for (s in 1:90) {
  saveRDS(df_list[[s]], paste0("radEmu/results/power_sims/full_setting", s, ".RDS"))
}