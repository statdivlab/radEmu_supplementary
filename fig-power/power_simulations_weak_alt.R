devtools::load_all(path="../radEmu/")

sim_dir_name <- "power_weak"

#parameters determining form of simulations
sample_sizes <- c(250,50,10)
Js <- c(250,50,10)
distns <- c("ZINB","Poisson")

#### change from 10 to 500
nsim <- 10

#parameters determining how to execute simulations
ncores <- 4 #cores to parallelize over
timeout <- 300 # how long (in seconds) to let an individual simulation run
# before quitting (to guard against stalled jobs)

#generate simulation settings and initialize directory to save results in
sim_settings <- expand.grid(sample_sizes,Js,distns)
colnames(sim_settings) <- c("n","J","distn")
sim_settings <- as.data.frame(sim_settings)

sim_results <- vector(nrow(sim_settings),
                      mode = "list")

dir.create(sim_dir_name)

set.seed(384093)
B_list <- lapply(Js,
                 function(j)
                   get_sim_bs(j))
for(i in 1:length(Js)){
  J <- Js[i]
  B_list[[i]]$b1[J/2] <- 1
}

results <- parallel::mclapply(1:nrow(sim_settings),
                              function(sim_no){
                                set.seed(290062 + sim_no)
                                message("Simulation Setting ", sim_no,":")
                                n <- sim_settings[sim_no,"n"]
                                J <- sim_settings[sim_no,"J"]
                                which_J <- which(Js ==J)
                                distn <- sim_settings[sim_no,"distn"]
                                message("Simulation setting ", sim_no,": J is ", J,", n is ",n,
                                        ", and distribution is ",distn,".")
                                b0 <- B_list[[which_J]]$b0
                                b1 <- B_list[[which_J]]$b1
                                X <- cbind(1,rep(0:1,each = n/2))
                                rezzies <- data.frame("wald_p" = numeric(nsim),
                                                      "score_p_null_info" = numeric(nsim),
                                                      "score_p_full_info" = numeric(nsim))
                                sim_name <- paste("n",
                                                  sim_settings[sim_no,"n"],
                                                  "J",
                                                  sim_settings[sim_no,"J"],
                                                  "distn",
                                                  sim_settings[sim_no,"distn"],
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
                                           verbose = TRUE,
                                           use_both_cov = FALSE,
                                           use_fullmodel_info = TRUE,
                                           return_both_score_pvals = TRUE)
                                    # )
                                  }
                                  )

                                  if(!inherits(fitted_model,"try-error")){
                                    rezzies[sim,] <- fitted_model$coef[J/2,c("wald_p",
                                                                             "score_pval_null_info",
                                                                             "score_pval_full_info")]

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
                              },
                              mc.cores = ncores,
                              mc.preschedule = FALSE)
saveRDS(results, "results_sims_weak_alt.RDS")
system("say done")

