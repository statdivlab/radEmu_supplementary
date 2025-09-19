library(tidyverse)

files <- list.files("radEmu/results/power_sims/", full.names = T)
keep <- str_detect(files, "additional")
files <- files[keep]
df_list <- vector(mode = "list", length = 108)

for (file in files) {
  res <- readRDS(file)
  setting <- readr::parse_number(stringr::str_remove(file, "radEmu/results/power_sims/additional_setting"))
  if (is.null(df_list[[setting]])) {
    df_list[[setting]] <- res
  } else {
    df_list[[setting]] <- rbind(df_list[[setting]], res)
  }
}

for (s in 1:108) {
  saveRDS(df_list[[s]], paste0("radEmu/results/power_sims/additional_full_setting", s, ".RDS"))
}