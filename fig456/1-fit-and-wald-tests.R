#### Scripts to reproduce the data analyses in the `radEmu` paper
####
#### David Clausen and Amy Willis
#### Dec 2023 -- Jan 2024

#### Part 1 -- fit the model and investigate some of the significant/
#### larger magnitude effect size estimates

####################################################
################ reanalysis with radEmu
####################################################

#### fit model

tt_ef <- system.time({
  full_fit_wald <- emuFit(Y = Y,
                          X = X,
                          tolerance = 0.01,
                          test_kj = data.frame(k = 2,
                                               j = 1:ncol(Y)),
                          return_wald_p = TRUE,
                          compute_cis = TRUE,
                          run_score_tests = FALSE,
                          verbose = TRUE)
})
tt_ef ###
# user   system  elapsed
# 1942.937  131.082 2103.753
## => 35 minutes for J = 845, n = 566, p = 12, J*p = 10140

full_fit_wald %>% names
full_fit_wald$coef %>% head


# saveRDS(full_fit_wald, "crc_reanalysis_full_fit_wald.RDS")
# full_fit_wald <- readRDS("crc_reanalysis_full_fit_wald.RDS")

coef_table <- full_fit_wald$coef %>%
  as_tibble %>%
  filter(covariate == "groupCRC") %>%
  mutate(wald_p = 2*pnorm(abs(estimate/se), lower.tail= FALSE)) %>%
  mutate(wald_p_BY = p.adjust(wald_p, method = "BY")) %>%
  mutate(detection = sapply(category_num,function(x) mean(Y[,x]>0)))

coef_table %>%
  filter(wald_p_BY<=0.005) %>%
  filter(detection > 0.3) %>%
  mutate(category = factor(category,levels = category[order(estimate)])) %>%
  ggplot() + geom_errorbar(aes(x = category,ymin = lower, ymax = upper),
                           width= 0.25) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 6))


###### dig in further
coef_table %>%
  filter(wald_p_BY<=0.005) %>%
  filter(detection > 0.3) %>%
  View()

coef_table %>%
  filter(wald_p_BY<=0.005) %>%
  View()

table(Y[,125]>0,meta$Group)

table(Y[,125]>0,meta$Group,meta$Study) %>% apply(3,function(x) x[2,]/colSums(x))

table(Y[,125]>0,meta$Study)  %>% (function(x) x[2,]/colSums(x))

meta$parv_m <- Y[,192] >0
table(meta$parv_m,meta$Group,meta$Study)

ages <- seq(min(meta$Age),max(meta$Age),1)
studies <- unique(meta$Study)
meta$bmi_cat <- cut(meta$BMI,c(min(meta$BMI) - 1,quantile(meta$BMI ,c(1/3,2/3)),max(meta$BMI)))
study_bmi <- expand.grid(studies,unique(meta$bmi_cat))
smooves <-
  lapply(1:nrow(study_bmi),function(ind){
    Y_fr_crc <- Y[meta$Study == study_bmi[ind,1]&meta$bmi_cat == study_bmi[ind,2]&meta$Group == "CRC",]
    age_fr_crc <- meta$Age[meta$Study == study_bmi[ind,1]&meta$bmi_cat == study_bmi[ind,2]&meta$Group == "CRC"]
    Y_fr_crc_sm <- do.call(cbind,
                           lapply(1:ncol(Y),
                                  function(j)
                                    sapply(ages,
                                           function(x){
                                             w <- dnorm(age_fr_crc - x,0,10);
                                             return(weighted.mean(Y_fr_crc[,j],w = w/sum(w)))
                                           })))
    Y_fr_crc_ps <- apply(Y_fr_crc_sm,1,function(x) pseudohuber_center(x,0.1))

    Y_fr_ctr <- Y[meta$Study == study_bmi[ind,1]&meta$bmi_cat == study_bmi[ind,2]&meta$Group == "CTR",]
    age_fr_ctr <- meta$Age[meta$Study == study_bmi[ind,1]&meta$bmi_cat == study_bmi[ind,2]&meta$Group == "CTR"]
    Y_fr_ctr_sm <- do.call(cbind,
                           lapply(1:ncol(Y),
                                  function(j)
                                    sapply(ages,
                                           function(x){
                                             w <- dnorm(age_fr_ctr - x,0,10);
                                             return(weighted.mean(Y_fr_ctr[,j],w = w/sum(w)))
                                           })))

    Y_fr_ctr_ps <- apply(Y_fr_ctr_sm,1,function(x) pseudohuber_center(x,0.1))
    return(data.frame(age = ages,
                      crc_bif = Y_fr_crc_sm[,125],
                      ctr_bif = Y_fr_ctr_sm[,125],
                      crc_ps = Y_fr_crc_ps,
                      ctr_ps = Y_fr_ctr_ps,
                      study = study_bmi[ind,1],
                      bmi_cat = study_bmi[ind,2]))})

do.call(rbind,smooves) %>%
  ggplot() +
  geom_line(aes(x = age,y = (crc_bif/crc_ps)/(ctr_bif/ctr_ps),group = bmi_cat,
                color= bmi_cat)) +
  scale_y_log10() +
  geom_abline(aes(slope = 0, intercept = 0),linetype = "dotted") +
  facet_wrap(~study) +
  theme_bw()

Y_rel <- diag(1/rowSums(Y)) %*% Y

