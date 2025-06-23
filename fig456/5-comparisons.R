#### Scripts to reproduce the data analyses in the `radEmu` paper
####
#### David Clausen and Amy Willis
#### Dec 2023 -- Jan 2024

#### Part 5 -- construct figure 6

####################################################
################ make figures
####################################################


library(tidyverse)
library(phyloseq)
library(corncob)
library(ALDEx2)
library(DESeq2)
library(magrittr)
# BiocManager::install("IFAA")
library(IFAA)
# BiocManager::install("ANCOMBC")
library(ANCOMBC)
# BiocManager::install("ExperimentSubset")
library(ExperimentSubset)

#########################
##### corncob #####
#########################
otus <- otu_table(Y, taxa_are_rows=F)
xx <- X %>% as.data.frame
rownames(xx) <- otus %>% sample_names
sd <- sample_data(xx)
ps <- phyloseq(sd, otus)
xx %>% colnames

# system.time({
  # dt_cc <- differentialTest(formula = ~ groupCRC + agespline1 + agespline2 + bmispline1 + bmispline2 +
  #                             genderM + samplingAfter + studyUS + studyCN + studyDE + studyAT,
  #                           formula_null = ~ agespline1 + agespline2 + bmispline1 + bmispline2 +
  #                             genderM + samplingAfter + studyUS + studyCN + studyDE + studyAT,
  #                           phi.formula= ~ samplingAfter + studyUS + studyCN + studyDE + studyAT,
  #                           phi.formula_null= ~ samplingAfter + studyUS + studyCN + studyDE + studyAT,
  #                           data=ps,
  #                           test="LRT")
# })
# system("say done")
# saveRDS(dt_cc, "wirbel-corncob.RDS")
dt_cc <- readRDS("wirbel-corncob.RDS")

#########################
##### radEmu
#########################
ps_compare <- readRDS("our_ps_compare.RDS")

#########################
##### ALDEx2
#########################

y_ald <- t(Y) ## J x n
x_ald <- as.matrix(xx) ## n x p
dim(x_ald)
dim(y_ald)

## "one way anova" lol
# aldex(y_ald, xx$groupCRC,
#       mc.samples=16, test="t", effect=TRUE,
#       include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL)

set.seed(231212)
# system.time({
#   ald_clr <- aldex.clr(y_ald, x_ald, mc.samples=128*4) ## works, 23 minutes
#   ald_glm <- aldex.glm(clr=ald_clr, x_ald, fdr.method="BH")
# })
# user   system  elapsed
# 1449.554   17.323 1483.525
# system("say done")
# saveRDS(ald_glm, "wirbel-aldex.RDS")

ald_glm <- readRDS("wirbel-aldex.RDS")

ald_ps <- ald_glm$`groupCRC:pval`
names(ald_ps) <- rownames(ald_glm)

#########################
##### DeSeq2
#########################

psdeseq <- phyloseq_to_deseq2(ps,
                              design = ~ agespline1 + agespline2 + bmispline1 + bmispline2 + genderM + samplingAfter + studyUS + studyCN + studyDE + studyAT + groupCRC)
# # calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(psdeseq), 1, gm_mean)
psdeseq <- estimateSizeFactors(psdeseq, geoMeans = geoMeans)

system.time({
  deseq_res <- DESeq(psdeseq,
                     full = ~ agespline1 + agespline2 + bmispline1 + bmispline2 +
                       genderM + samplingAfter + studyUS + studyCN + studyDE + studyAT + groupCRC,
                     reduced = ~ agespline1 + agespline2 + bmispline1 + bmispline2 +
                       genderM + samplingAfter + studyUS + studyCN + studyDE + studyAT,
                     fitType="local",
                     test="LRT",
                     parallel=TRUE)
})
# saveRDS(deseq_res, "wirbel-deseq2.RDS")
# dds <- readRDS("wirbel-deseq2.RDS")

# user  system elapsed
# 169.206  15.243  35.548
system("say done")
res <- results(deseq_res)


# res <- results(dds)

deseq_ps <- res$pvalue
names(deseq_ps) <- rownames(res)

#########################
##### IFAA #####
#########################

experimental_data <- SummarizedExperiment(assays = list("shotgun" = t(Y)),
                                          colData= X) # colData= X[, -1])????
experimental_data

# run IFAA
set.seed(123) # For full reproducibility
system.time({
  ifaa_results <- IFAA(experiment_dat = experimental_data,
                       testCov="groupCRC",
                       ctrlCov=colnames(X)[-c(1,2)]) ### otherwise missing
})
# user  system elapsed
# 19.194   2.402 279.032

# saveRDS(ifaa_results, "wirbel-ifaa.RDS")
ifaa_results <- readRDS("wirbel-ifaa.RDS")

ifaa_results %>% names
ifaa_results$full_results %>% names
cbind(Y %>% colnames,
      ifaa_results$full_results$taxon)

ifaa_results$full_results %>%
  as_tibble %>%
  filter(is.na(unadj.p.value))

ifaa_ps <- tibble("ifaa_taxon" = ifaa_results$full_results$taxon,
                  "p_ifaa" = ifaa_results$full_results$unadj.p.value) %>%
  mutate(which = stringdist::amatch(ifaa_taxon, colnames(Y), maxDist=8)) %>%
  mutate("taxon" = colnames(Y)[which])

ifaa_ps %>%
  filter(is.na(taxon))

ifaa_ps %>%
  filter(is.na(p_ifaa))

ifaa_ps %>%
  filter(is.na(p_ifaa))
colnames(Y)[209]

tibble("ifaa_taxon" = ifaa_results$full_results$taxon,
       "est" = ifaa_results$full_results$estimate,
       "p_ifaa" = ifaa_results$full_results$unadj.p.value) %>%
  mutate(which = stringdist::amatch(ifaa_taxon, colnames(Y), maxDist=8)) %>%
  mutate("category" = colnames(Y)[which]) %>%
  inner_join(ps_compare %>% dplyr::select(category, estimate, pval_score)) %>%
  ggplot(aes(x = estimate, y = est)) +
  geom_point() +
  ylab("IFAA effect size ests") + xlab("radEmu effect size ests")
# woah! not a good look, IFAA

#########################
##### ANCOMBC
#########################
experimental_data %>% class
tse %>% class

experimental_data2 <- TreeSummarizedExperiment(assays = list("shotgun" = t(Y)),
                                               colData= X[, -1])
# system.time({
#   ancom_results <- ancom(data = experimental_data2,
#                          assay_name = "shotgun",
#                          p_adj_method = "none",
#                          phyloseq = NULL,
#                          tax_level="Species",
#                          main_var = "groupCRC",
#                          # prv_cut = ### TODO,
#                          adj_formula = "agespline1 + agespline2 + bmispline1 + bmispline2 + genderM + samplingAfter + studyUS + studyCN + studyDE + studyAT",
#                          alpha = 0.05, n_cl = 6)
# })
# #   user  system elapsed
# # 2.444   1.961 840.977
# system("say done")
# 
# ancom_results %>% names
# ancom_results$p_data %>% class
# ancom_results$p_data %>% dim # 552 x 552
# ancom_results$res %>% as_tibble
# ancom_results$res %>% names
# ancom_results$beta_data %>% dim # 552 x 552
# 
# ancom_results$q_data %>% dim
# ancom_results$q_data[1:5, 1:5]

### I have no idea what to do with 552 p-values per taxon.

system.time({
  ancombc_results <- ancombc2(data = experimental_data2,
                              assay_name = "shotgun",
                              p_adj_method = "none",
                              tax_level="Species",
                              fix_formula= "groupCRC + agespline1 + agespline2 + bmispline1 + bmispline2 + genderM + samplingAfter + studyUS + studyCN + studyDE + studyAT",
                              prv_cut = 0,
                              alpha = 0.05,
                              n_cl = 6)
})
# user  system elapsed
# 31.220   1.616 297.254
saveRDS(ancombc_results, "wirbel-ancombc2.RDS")

ancombc_results %>% names


ancombc_ps <- tibble("category" = ancombc_results$res$taxon,
                     "p_ancombc2" = ancombc_results$res$p_groupCRC)

#########################
##### all
#########################

### find way to visualize that corncob has some NA's



ps_df <- ps_compare %>%
  inner_join(dt_cc$p %>% enframe(name="category", value = "p_corncob")) %>%
  inner_join(ald_ps %>% enframe(name="category", value = "p_aldex")) %>%
  inner_join(deseq_ps %>% enframe(name="category", value = "p_deseq2")) %>%
  inner_join(ancombc_ps) %>%
  inner_join(ifaa_ps %>% dplyr::rename(category = taxon)) %>%
  dplyr::inner_join(wirb_orig_small, by = c("category" = "species")) %>%
  mutate(t_r = -log10(pval_score),
         t_c = -log10(p_corncob),
         t_al = -log10(p_aldex),
         t_d = -log10(p_deseq2),
         t_an = -log10(p_ancombc2),
         t_i = -log10(p_ifaa),
         t_w = -log10(pval.meta))


rad_corn <- ps_df %>%
  mutate(`Missing p-value` = is.na(t_c)) %>%
  mutate(t_c = ifelse(is.na(t_c), 25, t_c)) %>%
  ggplot(aes(x = t_r, y = t_c)) +
  geom_point(cex = 0.1, aes(col = `Missing p-value`)) +
  geom_abline() +
  theme_bw() +
  ylim(0, 25) +
  ylab(expression(paste("-lo", g[10], "(", p[corncob], ")"), parse = TRUE)) +
  xlab(expression(paste("-lo", g[10], "(", p[radEmu], ")"), parse = TRUE)) +
  scale_color_manual(values=c("black",  "blue")) +
  NULL

rad_aldex <- ps_df %>%
  mutate(`Missing p-value` = is.na(t_al)) %>%
  mutate(t_al = ifelse(is.na(t_al), 25, t_al)) %>%
  ggplot(aes(x = t_r, y = t_al)) +
  geom_point(cex = 0.1, aes(col = `Missing p-value`)) +
  geom_abline() +
  theme_bw() +
  ylim(0, 25) +
  xlab(expression(paste("-lo", g[10], "(", p[radEmu], ")"), parse = TRUE)) +
  ylab(expression(paste("-lo", g[10], "(", p[ALDEx2], ")"), parse = TRUE)) +
  scale_color_manual(values=c("black",  "blue")) +
  NULL
summary(ps_df$p_aldex)

rad_deseq <- ps_df %>%
  mutate(`Missing p-value` = is.na(t_d)) %>%
  mutate(t_d = ifelse(is.na(t_d), 25, t_d)) %>%
  mutate(t_d = ifelse(t_d > 25, 25, t_d)) %>% 
  ggplot(aes(x = t_r, y = t_d)) +
  geom_point(cex = 0.1, aes(col = `Missing p-value`)) +
  geom_abline() +
  theme_bw() +
  ylim(0, 25) +
  xlab(expression(paste("-lo", g[10], "(", p[radEmu], ")"), parse = TRUE)) +
  ylab(expression(paste("-lo", g[10], "(", p[DeSeq2], ")"), parse = TRUE)) +
  scale_color_manual(values=c("black",  "blue")) +
  NULL
ps_df %>%
  filter(p_deseq2 < 10^{-25}) %>%
  select(category, p_deseq2, pval_score) %>%
  mutate(log10(pval_score))

rad_wilcoxon <- ps_df %>%
  mutate(`Missing p-value` = is.na(pval.meta)) %>%
  mutate(t_w = ifelse(is.na(t_w), 25, t_w)) %>%
  ggplot(aes(x= t_r, y = t_w)) +
  geom_point(cex = 0.1, aes(col = `Missing p-value`)) +
  geom_abline() +
  theme_bw() +
  ylim(0, 25) +
  xlab(expression(paste("-lo", g[10], "(", p[radEmu], ")"), parse = TRUE)) +
  ylab(expression(paste("-lo", g[10], "(", p[Wirbel~et~al], ")"), parse = TRUE)) +
  scale_color_manual(values=c("black",  "blue")) +
  NULL


rad_ifaa <- ps_df %>%
  mutate(`Missing p-value` = is.na(p_ifaa)) %>%
  mutate(t_i = ifelse(is.na(t_i), 25, t_i)) %>%
  ggplot(aes(x = t_r, y = t_i)) +
  geom_point(cex = 0.1, aes(col = `Missing p-value`)) +
  geom_abline() +
  theme_bw() +
  ylim(0, 25) +
  xlab(expression(paste("-lo", g[10], "(", p[radEmu], ")"), parse = TRUE)) +
  ylab(expression(paste("-lo", g[10], "(", p[IFAA], ")"), parse = TRUE)) +
  scale_color_manual(values=c("black",  "blue")) +
  NULL



rad_ancom <- ps_df %>%
  mutate(`Missing p-value` = is.na(p_ancombc2)) %>%
  mutate(t_an = ifelse(is.na(t_an), 25, t_an)) %>%
  mutate(t_an = ifelse(t_an > 25, 25, t_an)) %>%
  ggplot(aes(x = t_r, y = t_an)) +
  geom_point(cex = 0.1, aes(col = `Missing p-value`)) +
  geom_abline() +
  theme_bw() +
  ylim(0, 25) +
  xlab(expression(paste("-lo", g[10], "(", p[radEmu], ")"), parse = TRUE)) +
  ylab(expression(paste("-lo", g[10], "(", p[ANCOMBC2], ")"), parse = TRUE)) +
  scale_color_manual(values=c("black",  "blue")) +
  NULL


ggpubr::ggarrange(rad_corn, rad_aldex,
                  rad_deseq, rad_ifaa, rad_ancom, rad_wilcoxon, nrow = 2, ncol = 3, common.legend=TRUE)
ggsave("images/compare-da.pdf", width = 7, height = 3)

ps_df %>%
  summarise("cc" = mean(pval_score < p_corncob, na.rm=T) %>% multiply_by(100) %>% signif(2),
            "aldex" = mean(pval_score < p_aldex, na.rm=F) %>% multiply_by(100) %>% signif(2),
            "deseq" = mean(pval_score < p_deseq2, na.rm=T) %>% multiply_by(100) %>% signif(2),
            "ifaa" = mean(pval_score < p_ifaa, na.rm=T) %>% multiply_by(100) %>% signif(2),
            "ancombc2" = mean(pval_score < p_ancombc2, na.rm=T) %>% multiply_by(100) %>% signif(2),
            "wirb" = mean(pval_score < pval.meta, na.rm=T) %>% multiply_by(100) %>% signif(2)) %>%
  knitr::kable(format="latex")

c(ps_df %>% with(cor(p_corncob, pval_score,method = "spearman", use = "complete.obs")),
  ps_df %>% with(cor(p_aldex, pval_score,method = "spearman", use = "complete.obs")),
  ps_df %>% with(cor(p_deseq2, pval_score,method = "spearman", use = "complete.obs")),
  ps_df %>% with(cor(p_ifaa, pval_score,method = "spearman", use = "complete.obs")),
  ps_df %>% with(cor(p_ancombc2, pval_score,method = "spearman", use = "complete.obs")),
  ps_df %>% with(cor(pval.meta, pval_score,method = "spearman", use = "complete.obs"))) %>% round(2)

# ------------------- compare estimates -----------------------------------
est_df <- ps_compare %>%
  #inner_join(dt_cc$p %>% enframe(name="category", value = "p_corncob")) %>%
  inner_join(data.frame(category = rownames(ald_glm), ald_est = ald_glm$`groupCRC:Est`*log(2))) %>%
  inner_join(data.frame(category = rownames(res), deseq_est = res$log2FoldChange*log(2))) %>%
  inner_join(data.frame(category = ancombc_results$res$taxon, ancom_est = ancombc_results$res$lfc_groupCRC)) %>%
  inner_join(data.frame(category = ifaa_ps$taxon, ifaa_est = ifaa_results$full_results$estimate)) 
c(est_df %>% with(cor(ald_est, estimate, method = "spearman", use = "complete.obs")),
  est_df %>% with(cor(deseq_est, estimate, method = "spearman", use = "complete.obs")),
  est_df %>% with(cor(ancom_est, estimate, method = "spearman", use = "complete.obs")),
  est_df %>% with(cor(ifaa_est, estimate, method = "spearman", use = "complete.obs"))) %>% round(2)

est_long <- est_df %>%
  dplyr::select(estimate, ald_est, ancom_est, deseq_est, ifaa_est) %>%
  pivot_longer(cols = -estimate, names_to = "method", values_to = "value")
# Create ggplot
ggplot(est_long, aes(x = estimate, y = value)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~method) + 
  labs(
    x = "Estimate",
    y = "Values",
    color = "Method",
    title = "Comparison of Estimates"
  ) + 
  ylim(-5, 5) +
  geom_abline(aes(slope = 1, intercept = 0), color = "red")
