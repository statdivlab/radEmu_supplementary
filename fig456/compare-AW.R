
#### TODO remove

# all.equal(readRDS("crc_reanalysis_full_fit_wald.RDS")$coef,
#           readRDS("../../research/radEmu/re_wirb_j.RDS")$coef)

# summary(readRDS("crc_reanalysis_full_fit_wald.RDS")$coef$estimate -
# # Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# # -1.065e-01 -1.800e-06 -4.000e-08  1.190e-05  1.490e-06  1.065e-01


new <- readRDS("crc_reanalysis_full_fit_wald.RDS")$coef
orig <- readRDS("../../research/radEmu/re_wirb_j.RDS")$coef

max(abs(new$estimate - orig$estimate))
which.max(abs(new$estimate - orig$estimate))
max(abs(new$wald_p - orig$wald_p))
which.max(abs(new$wald_p - orig$wald_p))

