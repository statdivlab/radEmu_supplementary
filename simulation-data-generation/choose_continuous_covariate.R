# copy processing from 0-process.R 

# -----------------------------
library(tidyverse)
library(splines2)

### If you'd like to reproduce the results from the preprint, you should load
### release v1.0.0 (available on GitHub).
### We are going to ensure backwards compatibility of syntax, but if we make
### improvements to the software then future releases may result in different results (probably slight).
### We mention this just because we don't want you to be upset if you run a newer version
### and things are slightly different.

### As the developers, we just loaded it locally as follows.
# devtools::load_all(path="../radEmu/")


if(!("species_profiles_g2_l75_motus2.0.0.tsv" %in% list.files())){
  options(timeout = max(300, getOption("timeout")))
  download.file(url = "https://zenodo.org/records/3517209/files/species_profiles_g2_l75_motus2.0.0.tsv?download=1",
                destfile = "species_profiles_g2_l75_motus2.0.0.tsv",
                method = "libcurl")
}

if(!("meta_all.tsv" %in% list.files())){
  download.file(url = "https://zenodo.org/records/3517209/files/meta_all.tsv?download=1",
                destfile = "meta_all.tsv")
}

if(!("orig_wirbel_results.xlsx" %in% list.files())){
  download.file(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-019-0406-6/MediaObjects/41591_2019_406_MOESM3_ESM.xlsx",
                destfile = "orig_wirbel_results.xlsx")}

#mOTU count data
motu_counts <- read.csv("species_profiles_g2_l75_motus2.0.0.tsv", sep='\t')
motu_counts <- t(motu_counts)
# can immediately discard motus we do not observe at all
motu_counts <- motu_counts[,colSums(motu_counts)>0]

# as well as motu labeled "-1"
motu_counts <- motu_counts[,colnames(motu_counts) != "-1"]

#"metadata"
meta <- read.csv("meta_all.tsv",sep = "\t",header = TRUE)


studies <- c("FR-CRC","US-CRC",
             "CN-CRC","DE-CRC",
             "AT-CRC")


meta <- meta %>%
  filter(Study %in% studies) %>%
  filter(Group %in% c("CRC","CTR"))

meta %>%
  with(table(Sampling_rel_to_colonoscopy,Study))

#get id names as they will appear in rownames of motu_counts
clean_ids <- sapply(meta$Sample_ID,function(x) gsub("-",".",x,fixed = TRUE))

counts_to_keep <- sapply(rownames(motu_counts),function(x) x %in% clean_ids)

motu_counts <- motu_counts[counts_to_keep,]

#ensure that rows of motu_counts are in same order as rows of meta
stopifnot(mean(rownames(motu_counts) == clean_ids)==1)


#check completeness
keep_complete <- complete.cases(meta[,c("Study","Age","Gender","BMI","Sampling_rel_to_colonoscopy")])
meta <- meta[keep_complete,]
motu_counts <- motu_counts[keep_complete,]

age_spline <- splines2::bSpline(meta$Age,degree = 1,knots = median(meta$Age))
age_spline[,1] <- (age_spline[,1] - mean(age_spline[,1]))/sd(age_spline[,1])
age_spline[,2] <- (age_spline[,2] - mean(age_spline[,2]))/sd(age_spline[,2])

bmi_spline <- splines2::bSpline(meta$BMI,degree = 1, knots = median(meta$BMI))
bmi_spline[,1] <- (bmi_spline[,1] - mean(bmi_spline[,1]))/sd(bmi_spline[,1])
bmi_spline[,2] <- (bmi_spline[,2] - mean(bmi_spline[,2]))/sd(bmi_spline[,2])

X <- cbind(1,
           as.numeric(meta$Group == "CRC"),
           age_spline,
           bmi_spline,
           as.numeric(meta$Gender == "M"),
           as.numeric(meta$Sampling_rel_to_colonoscopy == "AFTER"),
           as.numeric(meta$Study == "US-CRC"),
           as.numeric(meta$Study == "CN-CRC"),
           as.numeric(meta$Study == "DE-CRC"),
           as.numeric(meta$Study == "AT-CRC"))

colnames(X) <-
  c("(Intercept)",
    "groupCRC",
    "agespline1",
    "agespline2",
    "bmispline1",
    "bmispline2",
    "genderM",
    "samplingAfter",
    "studyUS",
    "studyCN",
    "studyDE",
    "studyAT")

motu_RA <- motu_counts
motu_RA <- diag(1/rowSums(motu_counts)) %*% motu_counts
study_prevalence <- matrix(nrow = ncol(motu_counts),ncol = length(studies))

for(i in 1:length(studies)){
  study_prevalence[,i] <- colSums(motu_RA[meta$Study == studies[i],] >1e-3)
}

included <- which(rowSums(study_prevalence>0)>=3)

length(included)

Y <- motu_counts[,included]

# ---------------------

# design matrix

X[, "agespline1"] %>% min(); X[, "agespline1"] %>% max()
X[, "agespline2"] %>% min(); X[, "agespline2"] %>% max()
X[, "bmispline1"] %>% min(); X[, "bmispline1"] %>% max()
X[, "bmispline2"] %>% min(); X[, "bmispline2"] %>% max()
# overall, range ~ -4.3 - 5

# model fits
library(radEmu)
tt_ef <- emuFit(Y = Y,
                X = X,
                tolerance = 0.01,
                test_kj = data.frame(k = 2,
                                     j = 1:ncol(Y)),
                run_score_tests = FALSE,
                verbose = "development",
                compute_cis = FALSE,
                return_wald_p = FALSE)
est_as1 <- tt_ef$coef %>% filter(covariate == "agespline1") %>% pull(estimate); summary(est_as1)
est_as2 <- tt_ef$coef %>% filter(covariate == "agespline2") %>% pull(estimate); summary(est_as2)
est_bs1 <- tt_ef$coef %>% filter(covariate == "bmispline1") %>% pull(estimate); summary(est_bs1)
est_bs2 <- tt_ef$coef %>% filter(covariate == "bmispline2") %>% pull(estimate); summary(est_bs2)
# overall, range ~ -136 - 11

X[, "agespline1", drop = FALSE] %*% est_as1 %>% as.vector() %>% summary()
X[, "agespline2", drop = FALSE] %*% est_as2 %>% as.vector() %>% summary()
X[, "bmispline1", drop = FALSE] %*% est_bs1 %>% as.vector() %>% summary()
X[, "bmispline2", drop = FALSE] %*% est_bs2 %>% as.vector() %>% summary()
hist(c(X[, "agespline1", drop = FALSE] %*% est_as1 %>% as.vector(),
       X[, "agespline2", drop = FALSE] %*% est_as2 %>% as.vector(),
       X[, "bmispline1", drop = FALSE] %*% est_bs1 %>% as.vector(),
       X[, "bmispline2", drop = FALSE] %*% est_bs2 %>% as.vector()) %>% abs %>% log)
# range -700 - 91, mostly -30-30

# simulate! 
B2 <- radEmu:::get_sim_bs(250)$b1 * 2
X2 <- rnorm(50, mean = 0, sd = 2)

summary(matrix(B2, nrow = 250, ncol = 1) %*% matrix(X2, nrow = 1, ncol = 50) %>% as.vector)
matrix(B2, nrow = 250, ncol = 1) %*% matrix(X2, nrow = 1, ncol = 50) %>% as.vector %>% abs %>% log %>% hist
