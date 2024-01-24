#### Scripts to reproduce the data analyses in the `radEmu` paper
####
#### David Clausen and Amy Willis
#### Dec 2023 -- Jan 2024

#### Part 4 -- construct figure 6

####################################################
################ make figures
####################################################

ps_compare <- readRDS("our_ps_compare.RDS")

wirb_orig <- readxl::read_xlsx("orig_wirbel_results.xlsx",sheet = 2,skip=13)

ps_small <- ps_compare %>%
  dplyr::select(category,pval_score) %>%
  transmute(species = category,
            pval_score = pval_score)
wirb_orig_small <-  wirb_orig %>%
  dplyr::select(species,pval.meta)

# numbers for data analysis section text
ps_small %>%
  mutate(qval = p.adjust(pval_score,method = "BH")) %>%
  filter(qval<0.1) %>%
  mutate(detection = sapply(species, function(x) mean(Y[,x]>0))) %>%
  mutate(estimate = sapply(species, function(x) ps_compare$estimate[ps_compare$category==x])) %>%
  mutate(lower = sapply(species, function(x) ps_compare$lower[ps_compare$category==x])) %>%
  mutate(upper = sapply(species, function(x) ps_compare$upper[ps_compare$category==x])) %>%
  mutate(estimate = exp(estimate),
         lower = exp(lower),
         upper = exp(upper)) %>%
  View()

ps_small %>%
  mutate(qval = p.adjust(pval_score,method = "BH")) %>%
  filter(qval<0.1) %>%
  mutate(detection = sapply(species, function(x) mean(Y[,x]>0))) %>%
  mutate(estimate = sapply(species, function(x) ps_compare$estimate[ps_compare$category==x])) %>%
  mutate(lower = sapply(species, function(x) ps_compare$lower[ps_compare$category==x])) %>%
  mutate(upper = sapply(species, function(x) ps_compare$upper[ps_compare$category==x])) %>%
  mutate(estimate = exp(estimate),
         lower = exp(lower),
         upper = exp(upper)) %>%
  View()

wirb_orig_small$pval.meta[wirb_orig_small$pval.meta==0] <- 10^{-25}
  # wirb_orig_small %>% filter(pval.meta>0) %>% with(min(pval.meta))

dplyr::inner_join(ps_small,wirb_orig_small) %>%
  ggplot() +
  geom_point(aes(x = pval.meta,y= pval_score),alpha = 0.3) +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed",
              color = "#56B4E9") +
  xlab("Blocked Wilcoxon p-Value") +
  ylab("Robust Score p-Value") +
  theme_bw() +
  coord_equal()

#spearman corr
dplyr::inner_join(ps_small,wirb_orig_small) %>%
  filter(!is.na(pval.meta)) %>%
  with(cor(pval.meta,pval_score,method = "spearman"))

#taxon agreement
wirbel_taxa <-
  wirb_orig_small %>%
  mutate(wirb_q = p.adjust(pval.meta,method = "BH")) %>%
  filter(wirb_q < 0.005) %>%
  with(species)

top_94_score <- ps_small %>%
  mutate(p_rank = rank(pval_score)) %>%
  filter(p_rank <=94) %>%
  with(species)

#number in their top 94 we identify
sum(top_94_score %in% wirbel_taxa)

#number of our taxa they identify
dplyr::inner_join(ps_small,wirb_orig_small) %>%
  filter(!is.na(pval.meta)) %>%
  mutate(score_q = p.adjust(pval_score, method = "BH")) %>%
  mutate(wirb_q = p.adjust(pval.meta,method = "BH"))

dplyr::inner_join(ps_small,wirb_orig_small) %>%
  (function(x) x[complete.cases(x),]) %>%
  with(cor(pval.meta,pval_score,method = "spearman"))
