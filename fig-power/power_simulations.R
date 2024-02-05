devtools::load_all(path="../radEmu/")

sim_dir_name <- "power_curves"

#parameters determining form of simulations
sample_sizes <- c(100, 50, 30, 10)
Js <- c(500, 250)
distns <- c("ZINB")
effect_sizes <- seq(from = 0, to = 2.5, length.out = 11)

nsim <- 100

#parameters determining how to execute simulations
ncores <- 4 #cores to parallelize over
timeout <- 500 # how long (in seconds) to let an individual simulation run before quitting (to guard against stalled jobs)

#generate simulation settings and initialize directory to save results in
sim_settings <- expand.grid(sample_sizes,Js,distns, effect_sizes)
colnames(sim_settings) <- c("n","J","distn", "effect_sizes")
sim_settings <- as.data.frame(sim_settings)

sim_results <- vector(nrow(sim_settings),
                      mode = "list")

# dir.create(sim_dir_name)


my_sim_bs <- function(J){
  evens <- ((1:J)%%2 ==0)
  b0 <- numeric(J)
  b0[!evens] <- seq(-3,3,length.out = sum(!evens))
  b0[evens] <- seq(3,-3,length.out = sum(evens))

  ### alter to make smaller range
  b1 <- 1*sinh(seq(-10,10,length.out= J))/sinh(10)
  b1[(J/2):(J/2 + 1)] <- 0

  return(list(b0 = b0,
              b1 = b1))
}
B_list <- lapply(Js, my_sim_bs)

run_sims <- function(sim_no) {
  set.seed(290062 + sim_no)
  message("Simulation Setting ", sim_no,":")
  n <- sim_settings[sim_no,"n"]
  J <- sim_settings[sim_no,"J"]
  which_J <- which(Js ==J)
  distn <- sim_settings[sim_no,"distn"]
  effect_size <- sim_settings[sim_no, "effect_sizes"]
  message("Simulation setting ", sim_no,": J is ", J,", n is ",n,
          ", and effect size is ",effect_size,".")

  b0 <- B_list[[which_J]]$b0
  b1 <- B_list[[which_J]]$b1
  b1[J/2] <- effect_size

  X <- cbind(1,rep(0:1,each = n/2))
  rezzies <- data.frame("estimate" = numeric(nsim),
                        "wald_p" = numeric(nsim),
                        "score_p" = numeric(nsim))
  sim_name <- paste("n",
                    sim_settings[sim_no,"n"],
                    "J",
                    sim_settings[sim_no,"J"],
                    "beta",
                    sim_settings[sim_no,"effect_sizes"],
                    sep = "_",collapse= "_")
  for(sim in 1:nsim){
    message("Simulation ",sim)
    Y <- simulate_data(n = n,
                       J = J,
                       b0 = b0,
                       b1 = b1,
                       distn = distn,
                       zinb_size = 5,
                       zinb_zero_prop = 0.6,
                       mean_count_before_ZI = 50)

    fitted_model <- try({
      setTimeLimit(timeout)
      # suppressMessages(
      emuFit(Y = Y,
             X = X,
             test_kj = data.frame(k = 2, j = J/2),
             tau = 2,
             B_null_tol = 0.005,
             tolerance = 0.005,
             constraint_tol = 0.001,
             return_wald_p = TRUE,

             verbose = FALSE)

      ### these are defaults
             # use_both_cov = FALSE,
             # use_fullmodel_info = FALSE,
             # return_both_score_pvals = FALSE)
      # )
    }
    )

    if(!inherits(fitted_model,"try-error")){
      rezzies[sim,] <- fitted_model$coef[J/2,c("estimate", "wald_p", "pval")]
    } else{
      rezzies[sim,] <- NA
      if(!dir.exists(paste(sim_dir_name,"/",sim_name,sep = "",collapse = ""))){
        dir.create(paste(sim_dir_name,"/",sim_name,sep = "",collapse = ""))
      }
      saveRDS(list(Y = Y,
                   X = X,
                   b0 = b0,
                   b1 = b1,
                   sim = sim),
              paste(sim_dir_name,"/",sim_name,"/sim",sim,".rds",sep = "",collapse = "")
      )
    }

    if(sim%%10==0){
      saveRDS(rezzies[1:sim,],paste(sim_dir_name,"/",sim_name,".rds",sep = "",collapse= ""))
    }

  }
  return(rezzies)
}

# results <- parallel::mclapply(6:nrow(sim_settings),
#                               run_sims,
#                               mc.cores = 4,
#                               mc.preschedule = FALSE)

# results2 <- parallel::mclapply(nrow(sim_settings):65,
#                               run_sims,
#                               mc.cores = 4,
#                               mc.preschedule = FALSE)
results3 <- parallel::mclapply(c(64L, 61L, 57L, 63L, 59L, 62L, 58L),
                              run_sims,
                              mc.cores = 4,
                              mc.preschedule = FALSE)

system("say done")
# saveRDS(results, "results_sims_curves.RDS")

#### we recommend score_p_null_info

list.files("power_curves/")

to_tib <- function(x) {
  get_name <- x %>% str_split_i("/", 2)
  nn <- get_name %>% str_split_i("_", 2) %>% as.numeric
  jj <- get_name %>% str_split_i("_", 4) %>% as.numeric
  bbeta <- get_name %>% str_split_i("_", 6) %>% str_remove(".rds") %>% as.numeric
  x %>% readRDS %>% as_tibble %>% mutate(n = nn, j = jj, beta = bbeta)
}
"power_curves/n_50_J_250_beta_1.25.rds" %>% str_split_i("/", 2) %>% str_split_i("_", 6) %>% str_remove(".rds")

power_results <- list.files("power_curves", full.names=T, recursive=F) %>%
  sapply(to_tib, simplify=F) %>%
  do.call(bind_rows, .)

power_results %>%
  group_by(n, j, beta) %>%
  summarise(nn = n()) %>%
  arrange(desc(nn))

nrow(sim_settings)

power_results %>%
  mutate(n = paste("n =", n)) %>%
  mutate(`Sample size` = n) %>%
  mutate(j = paste("Number of taxa J =", j)) %>%
  group_by(n, j, beta, `Sample size`) %>%
  summarise(power = mean(score_p < 0.05),
            nn = n()) %>%
  mutate(lower = power - 1.96*sqrt(power*(1-power)/nn),
         upper = power + 1.96*sqrt(power*(1-power)/nn)) %>%
  ggplot(aes(x = beta, y = power, ymin = lower, ymax = upper, group = `Sample size`, col = `Sample size`)) +
  # ggplot(aes(x = beta, y = power, ymin = lower, ymax = upper, group = n, col = n)) +
  geom_linerange() +
  geom_line() +
  facet_wrap(~j) +
  ylab("Power") +
  xlab("Effect size") +
  theme_bw() +
  NULL
# ggsave("images/power_curves.pdf", height = 3, width = 7)

power_results %>%
  group_by(n, j, beta) %>%
  mutate(n = paste("n =", n)) %>%
  mutate(j = paste("j =", j)) %>%
  mutate(reject = score_p < 0.05) %>%
  summarise(power = mean(score_p < 0.05),
            nn = n()) %>%
  mutate(logit_power = log(power/(1-power))) %>%
  mutate(lower = power - 1.96*sqrt(power*(1-power)/nn),
         upper = power + 1.96*sqrt(power*(1-power)/nn)) %>%
  mutate(lower = log(lower/(1-lower)),
         upper = log(upper/(1-upper))) %>%
  # ggplot(aes(x = beta, y = logit_power, ymin = lower, ymax = upper, group = n, col = n)) +
  ggplot(aes(x = beta, y = logit_power, ymin = lower, ymax = upper, group = interaction(j, n), col = n, lty = j)) +
  # geom_linerange() +
  geom_line() +
  # facet_wrap(~j) +
  theme_bw() +
  NULL
### Evidence for effect modification between n and beta


power_glm <- power_results %>%
  group_by(n, j, beta) %>%
  mutate(reject = score_p < 0.05) %>%
  glm(reject ~ j + beta * n, data=., family = binomial(link = "logit"))

power_glm %>%
  summary

# Call:
#   glm(formula = reject ~ j + beta * n, family = binomial(link = "logit"),
#       data = .)
#
#
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -2.6040977  0.1490436 -17.472  < 2e-16 ***
#   j           -0.0002547  0.0002335  -1.091  0.27542
# beta         0.2209221  0.0771352   2.864  0.00418 **
#   n           -0.0036479  0.0020690  -1.763  0.07788 .
# beta:n       0.0295110  0.0015882  18.581  < 2e-16 ***

# saveRDS(power_glm, "power_glm.RDS")
