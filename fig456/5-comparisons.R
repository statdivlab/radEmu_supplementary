library(phyloseq)
library(corncob)
library(ALDEx2)
library(DESeq2)
library(magrittr)
# devtools::install_github("gitlzg/IFAA")
library(IFAA)

#########################
##### corncob #####
#########################
otus <- otu_table(Y, taxa_are_rows=F)
xx <- X %>% as.data.frame
rownames(xx) <- otus %>% sample_names
sd <- sample_data(xx)
ps <- phyloseq(sd, otus)
xx %>% colnames

system.time({
  dt_cc <- differentialTest(formula = ~ groupCRC + agespline1 + agespline2 + bmispline1 + bmispline2 +
                              genderM + samplingAfter + studyUS + studyCN + studyDE + studyAT,
                            formula_null = ~ agespline1 + agespline2 + bmispline1 + bmispline2 +
                              genderM + samplingAfter + studyUS + studyCN + studyDE + studyAT,
                            phi.formula= ~ samplingAfter + studyUS + studyCN + studyDE + studyAT,
                            phi.formula_null= ~ samplingAfter + studyUS + studyCN + studyDE + studyAT,
                            data=ps,
                            test="LRT")
})
system("say done")
# saveRDS(dt_cc, "wirbel-corncob.RDS")
# dt_cc <- readRDS("wirbel-corncob.RDS")

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
system.time({
  ald_clr <- aldex.clr(y_ald, x_ald, mc.samples=128*4) ## works, 23 minutes
  ald_glm <- aldex.glm(clr=ald_clr, x_ald, fdr.method="BH")
})
#     user   system  elapsed
# 1528.093   29.369 2242.442
system("say done")

ald_ps <- ald_glm$`groupCRC:pval`
names(ald_ps) <- rownames(ald_glm)

#########################
##### DeSeq2
#########################

psdeseq <- phyloseq_to_deseq2(ps, ~ groupCRC + agespline1 + agespline2 + bmispline1 + bmispline2 +
                                genderM + samplingAfter + studyUS + studyCN + studyDE + studyAT)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(psdeseq), 1, gm_mean)
psdeseq <- estimateSizeFactors(psdeseq, geoMeans = geoMeans)
dds <- estimateDispersions(psdeseq)
# system.time({dds <- nbinomWaldTest(dds, maxit = 1000)})
system.time({dds <- nbinomWaldTest(dds, maxit = 20000)}) ## 12 minutes for 20000 maxit
# user  system elapsed
# 659.127   1.547 664.093
# found results columns, replacing these
# 27 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest
res <- results(dds)

deseq_ps <- res$pvalue
names(deseq_ps) <- rownames(psdeseq)

#########################
##### all
#########################

### find way to visualize that corncob has some NA's



ps_df <- ps_compare %>%
  inner_join(dt_cc$p %>% enframe(name="category", value = "p_corncob")) %>%
  inner_join(ald_ps %>% enframe(name="category", value = "p_aldex")) %>%
  inner_join(deseq_ps %>% enframe(name="category", value = "p_deseq2")) %>%
  mutate(t_r = -log10(pval_score),
         t_c = -log10(p_corncob),
         t_a = -log10(p_aldex),
         t_d = -log10(p_deseq2))


ps_df
ps_df %>%
  filter(is.na(t_c))
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
  mutate(`Missing p-value` = is.na(t_a)) %>%
  mutate(t_a = ifelse(is.na(t_a), 25, t_a)) %>%
  ggplot(aes(x = t_r, y = t_a)) +
  geom_point(cex = 0.1, aes(col = `Missing p-value`)) +
  geom_abline() +
  theme_bw() +
  ylim(0, 25) +
  xlab(expression(paste("-lo", g[10], "(", p[radEmu], ")"), parse = TRUE)) +
  ylab(expression(paste("-lo", g[10], "(", p[ALDEx2], ")"), parse = TRUE)) +
  scale_color_manual(values=c("black",  "blue")) +
  NULL

rad_deseq <- ps_df %>%
  mutate(`Missing p-value` = is.na(t_d)) %>%
  mutate(t_d = ifelse(is.na(t_d), 25, t_d)) %>%
  ggplot(aes(x = t_r, y = t_d)) +
  geom_point(cex = 0.1, aes(col = `Missing p-value`)) +
  geom_abline() +
  theme_bw() +
  ylim(0, 25) +
  xlab(expression(paste("-lo", g[10], "(", p[radEmu], ")"), parse = TRUE)) +
  ylab(expression(paste("-lo", g[10], "(", p[DeSeq2], ")"), parse = TRUE)) +
  scale_color_manual(values=c("black",  "blue")) +
  NULL

rad_wilcoxon <- dplyr::inner_join(ps_df, wirb_orig_small, by = c("category" = "species")) %>%
  mutate(`Missing p-value` = is.na(pval.meta)) %>%
  ggplot(aes(y = -log10(pval.meta), x= -log10(pval_score))) +
  geom_point(cex = 0.1, aes(col = `Missing p-value`)) +
  geom_abline() +
  theme_bw() +
  ylim(0, 25) +
  xlab(expression(paste("-lo", g[10], "(", p[radEmu], ")"), parse = TRUE)) +
  ylab(expression(paste("-lo", g[10], "(", p[Wirbel~et~al], ")"), parse = TRUE)) +
  scale_color_manual(values=c("black",  "blue")) +
  NULL

ggpubr::ggarrange(rad_wilcoxon, rad_corn, rad_aldex, rad_deseq, nrow = 1, common.legend=TRUE)
# ggsave("images/compare-da.pdf", width = 7, height = 2)


ps_df %>%
  filter(p_deseq2 < 10^{-25}) %>%
  dplyr::select(category, p_deseq2, pval_score)

ps_df %>%
  dplyr::inner_join(wirb_orig_small, by = c("category" = "species")) %>%
  summarise("wirb" = mean(pval_score < pval.meta, na.rm=T) %>% multiply_by(100) %>% signif(2),
            "cc" = mean(pval_score < p_corncob, na.rm=T) %>% multiply_by(100) %>% signif(2),
            "aldex" = mean(pval_score < p_aldex, na.rm=F) %>% multiply_by(100) %>% signif(2),
            "deseq" = mean(pval_score < p_deseq2, na.rm=T) %>% multiply_by(100) %>% signif(2)) %>%
  knitr::kable(format="latex")

#########################
##### IFAA #####
#########################

experimental_data <- SummarizedExperiment(assays = list("shotgun" = t(Y)),
                                          colData= X[, -1])
experimental_data

# run IFAA
set.seed(123) # For full reproducibility
ifaa_results <- IFAA(experiment_dat = experimental_data,
                     testCov="groupCRC",
                     ctrlCov=colnames(X)[-c(1,2)])
### 5.22 mins?
