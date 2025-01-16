library(tidyverse)
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
                 deseq_est = NA,
                 deseq_p = NA)

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
saveRDS(df, paste0("radEmu/wirbel_permute/results/deseq_perm", perm, ".rds"))
