library(tidyverse)
library(magrittr)
devtools::load_all(path="../radEmu/")

set.seed(1)
n <- 10
J <- 6
bs <- get_sim_bs(J)
bs$b0 <- rep(0, J)
bs$b1 <- c(10, 1, 1, 1, 1, -3)/8
bs$b1 - mean(bs$b1)
Y_pois <- simulate_data(n = n,
                        J = J,
                        b0 = bs$b0,
                        b1 = bs$b1,
                        distn = "Poisson",
                        mean_count_before_ZI = 250)
X <- cbind(1, rep(0:1,each = n/2))

pois_mod_mean <- emuFit(Y = Y_pois,
                        X = X,
                        return_wald_p = FALSE,
                        run_score_tests = TRUE,
                        verbose = FALSE,
                        use_both_cov = FALSE,
                        constraint_fn = mean,
                        constraint_grad_fn = (function(x) rep(1/length(x),length(x))),
                        constraint_param=NA)

pois_mod_ps <-  emuFit(Y = Y_pois,
                       X = X,
                       return_wald_p = FALSE,
                       run_score_tests = TRUE,
                       verbose = FALSE,
                       use_both_cov = FALSE)

the_tb <- bind_rows(pois_mod_mean$coef %>% mutate(null = "Mean Centering"),
                    pois_mod_ps$coef %>% mutate(null = "Smoothed Median Centering")) %>%
  tibble %>%
  dplyr::select(-c(covariate, category)) %>%
  mutate(true_B = rep(bs$b1,2)) %>%
  mutate(true_B_centered = ifelse(null == "Mean Centering", bs$b1 - mean(bs$b1), bs$b1 - pseudohuber_center(bs$b1,0.1))) %>%
  mutate_if(is.integer, as.character, na.rm = TRUE) %>%
  pivot_longer(c(estimate, true_B, true_B_centered), names_to="Type", values_to="B") %>%
  mutate(score_stat = round(score_stat, 2),
         pval = round(pval, 2)) %>%
  filter(Type == "true_B") %>%
  mutate(category_num = as.numeric(category_num) + ifelse(Type == "estimate", 0.1, -0.1)) %>%
  mutate(Type = ifelse(Type == "estimate", "Estimate", "True"))

the_tb %>%
  mutate(B_mean = B - mean(B)) %>%
  mutate(B_pH = B - pseudohuber_center(B, 0.1)) %>%
  pivot_longer(B:B_pH) %>%
  mutate(name = ifelse(name == "B",
                       "Absolute (Uncentered)",
                       ifelse(name == "B_mean",
                              "Mean Centering", "Smoothed Median Centering"))) %>%
  ggplot(aes(x = category_num, y = value, color = Type, ymin = lower, ymax = upper, label = pval)) +
  facet_grid(~name) +
  geom_text(aes(fontface = ifelse(pval < 0.05, 2, 1)),
            data = . %>% filter(Type == "Estimate"), nudge_x=0.3, show.legend = FALSE) +
  geom_linerange(data = . %>% filter(Type == "Estimate"), show.legend=F) +
  geom_point() +
  theme_bw() +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'), #transparent legend panel
        legend.position= "none") +
  xlab(expression(paste("Taxon index ", j), parse = TRUE)) +
  scale_color_manual(values=c("#0d57d1")) + # "#208f1a",
  ylab(expression(paste("True ", beta[j]), parse = TRUE)) +
  # ylab(expression(atop(paste("True ", beta[j], ": Smoothed"), paste("Median Centering")), parse = TRUE)) +
  scale_x_continuous(breaks=seq(1, 7, 1), minor_breaks=NULL) +
  ylim(-1.25, 1.25) +
  geom_hline(yintercept=0, lty = 2, col = "#49484a", alpha = 0.7) +
  NULL

# ggsave("images/compare-centering.pdf",
#        width=7, height=2.3)
