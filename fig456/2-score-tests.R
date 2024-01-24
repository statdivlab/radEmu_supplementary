#### Scripts to reproduce the data analyses in the `radEmu` paper
####
#### David Clausen and Amy Willis
#### Dec 2023 -- Jan 2024

#### Part 2 -- run score tests

####################################################
################ score tests
####################################################


## Because they are J = 845, n = 566, p = 12, J*p = 10140
## i.e., a huge number of parameters to estimate,
## and the score test requires refitting under the null,
## this analysis is computationally expensive

## We provide both a parallelized script as well as precomputed
## p-values so that you don't need to run it to reproduce the figures

## To allow parallelization, we call `fit_null` and `get_score_stat`
## outside of the formal function.

## Please let us know if you'd like us to prioritize a wrapper for this functionality

full_model <- readRDS("crc_reanalysis_full_fit_wald.RDS")

X_cup <- X_cup_from_X(X, J = ncol(Y))

B_cup <- B_cup_from_B(full_model$B)

#### compute information for use in fitting null models
I <- f_info(Y = Y,
            B_cup = B_cup,
            B = full_model$B,
            X = X,
            X_cup = X_cup,
            omit = NULL,
            compute_together = FALSE,
            return_separate = FALSE)

#### compute Dy for use in fitting null models
n <- nrow(Y)
J <- ncol(Y)
p <- ncol(X)
j_ref <- which.max(colSums(Y>0)) #should be 161 for bacteroides dorei
indexes_to_remove <- (j_ref - 1)*p + 1:p #indexes for elements of B_cup to be removed under constraint
scores <- vector(n,mode = "list")

B <- full_model$B
for(k in 1:p){
  B[k,] <- B[k,] - B[k,j_ref]
}
z <- update_z(Y = Y,
              X = X,
              B = B)

for(i in 1:n){
  X_cup_i <- X_cup[(i - 1)*J + 1:J,]
  scores[[i]] <- Matrix::crossprod(X_cup_i,Y[i,] - exp(X_cup_i%*%B_cup + z[i]))

}

score_mat <- do.call(cbind,scores)
score_mat <- methods::as(score_mat,"sparseMatrix")
Dy <- Matrix::tcrossprod(score_mat)

dir.create("wirbel_score_tests")

run_jth_score_test <- function(j){

  null_fit <-  fit_null(B = full_model$B, #B (MPLE)
                        Y = full_model$Y_augmented, #Y (with augmentations)
                        X = X, #design matrix
                        X_cup = X_cup,
                        k_constr = 2, #row index of B to constrain
                        j_constr = j, #col index of B to constrain
                        j_ref = j_ref,
                        constraint_fn = (function(x) pseudohuber_center(x,0.1)), #constraint function
                        constraint_grad_fn =(function(x) dpseudohuber_center_dx(x,0.1)), #gradient of constraint fn\
                        rho_init = 1,
                        tau = 5,
                        constraint_tol = 1e-3, ## changed dec 11
                        verbose = TRUE,
                        trackB = FALSE,
                        # I = I,
                        # Dy = Dy
  )

  score_stat <- get_score_stat(Y = full_model$Y_augmented,
                               X_cup = X_cup,
                               X = X,
                               B = null_fit$B,
                               k_constr = 2,
                               j_constr = j,
                               constraint_grad_fn = function(x) dpseudohuber_center_dx(x,0.1),
                               indexes_to_remove = indexes_to_remove,
                               j_ref = j_ref,
                               J = J,
                               n = n,
                               p = p)

  file_loc <- paste("wirbel_score_tests/test_",j,".rds",sep = "",collapse =)
  saveRDS(list(null_fit = null_fit,
               score_stat = score_stat),
          file_loc)


}

##### Time and run a random sample of 16 taxa

# system.time({
#   run_jth_score_test(1) ### 15 mins
# })
#
# js_todo <- sample(1:845, 16)


#### They take about 15 minutes each on my computer, so 211 hours in series or
#### ~2 days in parallel over 4 cores. Be careful before running the next line!


# parallel::mclapply(1:845, FUN=run_jth_score_test, mc.cores=4)



#
# #### Read in the results
# get_score_stat_from_rds <- function(x) {
#   res <- try(readRDS(x)$score_stat$score_stat[1,1], silent=T)
#   if (class(res) == "try-error") {
#     res <- try(readRDS(x)$score_stat)
#   }
#   res
# }
#
# ####
#
# score_stats_tb <- list.files("wirbel_score_tests", full.names=T) %>%
#   sapply(get_score_stat_from_rds) %>%
#   enframe(name="file", value="test_stat") %>%
#   mutate("pval_score" = 1 - pchisq(test_stat, 1)) %>%
#   mutate(category_num = gsub(".*_","",file)) %>%
#   mutate(category_num = str_remove(category_num, ".rds")) %>%
#   mutate(category_num = as.numeric(category_num)) %>%
#   arrange(category_num)
#
# # score_stats_tb %>% inner_join(ps_compare, by = "category_num") %>%
# #   dplyr::select(category_num, test_stat.x, test_stat.y, pval_score.x, pval_score.y)
#
#
# my_ps_compare <- full_model$coef %>%
#   tibble %>%
#   filter(covariate == "groupCRC") %>%
#   inner_join(score_stats_tb, by = "category_num") %>%
#   dplyr::select(-covariate, -score_stat, -pval) ## remove empty columns
#
#
# ps_compare %<>%
#   mutate(score_p_BY = p.adjust(pval_score, method = "BY")) %>%
#   mutate(wald_p_BY = p.adjust(wald_p, method = "BY")) %>%
#   mutate(detection = sapply(category_num, function(x) mean(Y[ , x] > 0)))
#
# saveRDS(object = my_ps_compare, "my_ps_compare.RDS")

#### Alternatively, you can read in the data that we have run for you

ps_compare <- readRDS("our_ps_compare.RDS")

ps_compare %>%
  ggplot(aes(x = wald_p, y = pval_score)) +
  geom_point() +
  xlab("Wald p-values") +
  ylab("Score p-values")

ps_compare %>% pull(pval_score) %>% summary
ps_compare %>% pull(pval_score) %>% sort %>% head(50)
ps_compare %>% pull(score_p_BY) %>% summary
ps_compare %>% pull(score_p_BY) %>% sort %>% head(50)
