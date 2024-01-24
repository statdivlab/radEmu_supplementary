#### Scripts to reproduce the data analyses in the `radEmu` paper
####
#### David Clausen and Amy Willis
#### Dec 2023 -- Jan 2024

#### Part 3 -- construct figures 4 and 5

####################################################
################ make figures
####################################################


ps_compare <- readRDS("our_ps_compare.RDS")

ps_compare %<>%
  mutate(qvals = p.adjust(pval_score,method = "BH")) %>%
  mutate(significant = qvals<0.1)

ps_compare %>%
  filter(significant) %>%
  arrange(desc(estimate)) %>%
  dplyr::select(category, estimate, pval_score, lower, upper) %>%
  mutate(est_exp = signif(exp(estimate), 2),
         lower_exp = signif(exp(lower), 2),
         upper_exp = signif(exp(upper), 2))


# ps_compare %>%
#   ggplot() +
#   geom_errorbar(aes(x = detection,ymin = lower, ymax = upper,
#                     color = significant),
#                 width = 0) +
#   scale_color_manual(values=c("#56B4E9","#3446eb")) +
#   xlab("Proportion of participants in whom taxon is detected") +
#   ylab(bquote(95~"%"~robust~Wald~confidence~intervals~"for"~{beta^j}[CRC]))+
#   theme_bw() +
#   guides(color = guide_legend(title = "BH-Adjusted\nRobust score \np-value < 0.1")) +
#   theme(legend.box.background = element_rect(color = "black"))

ps_compare %>%
  ggplot() +
  geom_linerange(aes(x = detection,
                     color = significant,
                     alpha = significant,
                     ymin = lower,
                     ymax = upper)) +
  ylab(expression(paste("Estimate and robust Wald confidence interval for ", beta[CRC]^j), parse = TRUE)) +
  xlab(expression(paste("Proportion of participants in which taxon is detected"), parse = TRUE)) +
  labs(color = "Score test significant at FDR < 0.25",
       alpha = "Significant at FDR 0.25 (score p)")+
  theme_bw() +
  # ylab(expression(paste("Estimate and robust Wald\nconfidence interval for", beta[CRC]^j), parse = TRUE)) +
  ylab(bquote(95~"%"~robust~Wald~confidence~intervals~"for"~{beta^j}[CRC]))+
  xlab(expression(paste("Proportion of participants in which taxon is detected"), parse = TRUE)) +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(),
    # legend.position="right",
    legend.position=c(.8,.85),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  ) +
  geom_hline(yintercept=0, lty=2) +
  scale_color_manual(values=c("#93b8f5",  "#000080")) + # "#208f1a",
  scale_alpha_manual(values=c(0.85, 1)) +
  labs(color = expression(paste("BH-adjusted ", p[score], " < 0.1"), parse = TRUE)) +
  labs(alpha = expression(paste("BH-adjusted ", p[score], " < 0.1"), parse = TRUE)) +
  # labs(color = expression(paste(FDR[score], " < 0.1"), parse = TRUE)) +
  # labs(alpha = expression(paste(FDR[score], " < 0.1"), parse = TRUE)) +
  NULL

# ggsave("images/wirbel_v_detection.pdf",
#        width=7, height = 4)

ps_compare %>%
  mutate(detection = sapply(category,
                            function(x) mean(Y[,x]>0))) %>%
  mutate(detection = cut(detection, c(0,0.01,0.05,0.1,0.25,0.5,0.75,1))) %>%
  mutate(qvals = p.adjust(pval_score,method = "BH")) %>%
  filter(qvals<0.1) %>%
  mutate(taxon = sapply(category,
                        function(x)
                          strsplit(x," [",fixed = TRUE)[[1]][1])) %>%
  mutate(taxon_order = sapply(taxon,
                              function(x)
                                ifelse(grepl("unknown",x,fixed = TRUE),
                                       paste(
                                         strsplit(x,"unknown ")[[1]][2],
                                         " (unknown)"
                                       ),
                                       x))) %>%
  (function(x){
    x$taxon[x$taxon == "Bifidobacterium catenulatum/kashiwanohense"] <-
      "Bifidobacterium kashiwanohense";
    return(x)
  }) %>%
  mutate(taxon = factor(taxon,levels = unique(taxon)[order(unique(taxon_order))])) %>%
  ggplot() +
  geom_errorbar(aes(x = taxon,ymin = lower, ymax = upper,
                    color = detection), width = 0,
                position = position_dodge(width = 0.5)) +
  theme_bw() +
  scale_color_viridis_d("Detection \nCategory",
                        option = "viridis",
                        begin = 0.1,
                        end = 0.8,
                        direction = -1) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1, size = 6),
        axis.title.y = element_text(vjust = -1)) +
  xlab(bquote("Taxon:"~j)) +
  ylab(bquote("95% robust Wald CI for"~{beta^j}[CRC])) +
  geom_abline(slope = 0, intercept=0, lty = 2)

# ggsave("images/significant_wirbel.pdf",
#        width=7, height = 3.5)




