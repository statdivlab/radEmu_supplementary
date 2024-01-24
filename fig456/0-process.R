#### Scripts to reproduce the data analyses in the `radEmu` paper
####
#### David Clausen and Amy Willis
#### Dec 2023 -- Jan 2024

#### Part 0 -- download the dataset and process it to obtain the
#### observed count matrix and covariate matrix


############################
###### Load Libraries ######
############################
library(tidyverse)
library(splines2)

devtools::load_all(path="../radEmu/")

############################
###### Download Data ######
############################

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


############################
###### Load Data ######
############################
#mOTU count data
motu_counts <- read.csv("species_profiles_g2_l75_motus2.0.0.tsv", sep='\t')
motu_counts <- t(motu_counts)
# can immediately discard motus we do not observe at all
motu_counts <- motu_counts[,colSums(motu_counts)>0]

# as well as motu labeled "-1"
motu_counts <- motu_counts[,colnames(motu_counts) != "-1"]

#"metadata"
meta <- read.csv("meta_all.tsv",sep = "\t",header = TRUE)


###### Select Relevant Study Data ######
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

###### Covariates ######

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

############################################################
################## Taxa to Include #########################
############################################################

### used same taxon filtering criteria as Wirbel et al (to be consistent)

motu_RA <- motu_counts
motu_RA <- diag(1/rowSums(motu_counts)) %*% motu_counts
study_prevalence <- matrix(nrow = ncol(motu_counts),ncol = length(studies))

for(i in 1:length(studies)){
  study_prevalence[,i] <- colSums(motu_RA[meta$Study == studies[i],] >1e-3)
}

included <- which(rowSums(study_prevalence>0)>=3)

length(included)

Y <- motu_counts[,included]


### yields Y and X

