library(tidyverse)
library(TreeSummarizedExperiment)
library(ANCOMBC)
library(ALDEx2)
library(DESeq2)
library(phyloseq)

# get permutation
args <- commandArgs(trailingOnly = FALSE)
if (length(args) == 0) {
  perm <- 1
} else {
  arg <- args[length(args)]
  perm <- abs(readr::parse_number(arg))
}

# get X and Y matrices
Y <- readRDS("radEmu/wirbel_permute/data/Y.rds")
X_perm <- readRDS(paste0("radEmu/wirbel_permute/data/X_perm", perm, ".rds"))

# save results 
df <- data.frame(perm = perm,
                 ind = 1:ncol(Y),
                 tax = colnames(Y),
                 aldex_est = NA,
                 aldex_p = NA,
                 ancom_est = NA,
                 ancom_p = NA,
                 ancom_pass_ss = NA,
                 clr_lm_est = NA,
                 clr_lm_p = NA,
                 deseq_est = NA,
                 deseq_p = NA)

# ALDEx2 
y_ald <- t(Y) 
x_ald <- as.matrix(X_perm)
ald_clr <- aldex.clr(y_ald, x_ald, mc.samples = 128*4)
ald_glm <- aldex.glm(clr = ald_clr, x_ald)
df$aldex_est <- ald_glm$`groupCRC:Est`
df$aldex_p <- ald_glm$`groupCRC:pval`
saveRDS(df, paste0("radEmu/wirbel_permute/results/other_methods_perm", perm, ".rds"))

# ANCOM-BC
ancom_dat <- TreeSummarizedExperiment(assays = list("shotgun" = t(Y)),
                                      colData= X_perm[, -1])
ancom_mod <- ancombc2(data = ancom_dat,
                      assay_name = "shotgun",
                      p_adj_method = "none",
                      fix_formula= "groupCRC + agespline1 + agespline2 + bmispline1 + bmispline2 + genderM + samplingAfter + studyUS + studyCN + studyDE + studyAT",
                      prv_cut = 0,
                      alpha = 0.05, 
                      pseudo_sens = TRUE)
df$ancom_est <- ancom_mod$res$lfc_groupCRC
df$ancom_p <- ancom_mod$res$p_groupCRC 
df$ancom_pass_ss <- ancom_mod$res$passed_ss_groupCRC
saveRDS(df, paste0("radEmu/wirbel_permute/results/other_methods_perm", perm, ".rds"))

# CLR linear model
Y_pseudo <- Y + 1
log_Y <- log(Y_pseudo)
mean_log_Y <- rowMeans(log_Y)
clr_data <- log_Y - matrix(mean_log_Y, nrow = nrow(Y), ncol = ncol(Y), byrow = FALSE)
for (j in 1:ncol(Y)) {
  print(j)
  mod <- lm(y ~ groupCRC + agespline1 + agespline2 + bmispline1 + bmispline2 + genderM + samplingAfter + studyUS + studyCN + studyDE + studyAT,
            data = data.frame(y = clr_data[, j], X_perm))
  df$clr_lm_est[j] <- coef(summary(mod))[2, 1]
  df$clr_lm_p[j] <- coef(summary(mod))[2, 4]
}
saveRDS(df, paste0("radEmu/wirbel_permute/results/other_methods_perm", perm, ".rds"))

# DESeq2 
otus <- otu_table(Y, taxa_are_rows=F)
xx <- X_perm %>% as.data.frame
rownames(xx) <- otus %>% sample_names
sd <- sample_data(xx)
ps <- phyloseq(sd, otus)
psdeseq <- phyloseq_to_deseq2(ps,
                              design = ~ agespline1 + agespline2 + bmispline1 + 
                                bmispline2 + genderM + samplingAfter + studyUS + 
                                studyCN + studyDE + studyAT + groupCRC)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(psdeseq), 1, gm_mean)
psdeseq <- estimateSizeFactors(psdeseq, geoMeans = geoMeans)
deseq_res <- DESeq(psdeseq,
                   full = ~ agespline1 + agespline2 + bmispline1 + bmispline2 +
                     genderM + samplingAfter + studyUS + studyCN + studyDE + studyAT + groupCRC,
                   reduced = ~ agespline1 + agespline2 + bmispline1 + bmispline2 +
                     genderM + samplingAfter + studyUS + studyCN + studyDE + studyAT,
                   fitType="local",
                   test="LRT",
                   parallel=FALSE)
df$deseq_est <- results(deseq_res)$log2FoldChange
df$deseq_p <- results(deseq_res)$pvalue
saveRDS(df, paste0("radEmu/wirbel_permute/results/other_methods_perm", perm, ".rds"))
