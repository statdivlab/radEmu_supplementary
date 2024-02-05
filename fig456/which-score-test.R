devtools::load_all("/Users/adwillis/research/radEmu/")

full_fit_wald <- readRDS("../../research/radEmu/re_wirb_j.RDS")
full_model <- full_fit_wald

X_cup <- X_cup_from_X(X, J = ncol(Y))

B_cup <- B_cup_from_B(full_model$B)

dim(Y)
dim(full_model$B)

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

test1_with_I_Dy <- function(j){

  init_obj <- readRDS(paste("local-results/","test_",j,".rds",sep = "",collapse =))
  # init_obj$null_fit$u

  null_fit <-  fit_null(B = init_obj$null_fit$B,
                        #B (MPLE)
                        Y = full_model$Y_augmented, #Y (with augmentations)
                        X = X, #design matrix
                        X_cup = X_cup,
                        k_constr = 2, #row index of B to constrain
                        j_constr = j, #col index of B to constrain
                        j_ref = j_ref,
                        constraint_fn = (function(x) pseudohuber_center(x,0.1)), #constraint function
                        constraint_grad_fn =(function(x) dpseudohuber_center_dx(x,0.1)), #gradient of constraint fn\
                        rho_init = init_obj$null_fit$rho,
                        tau = 5,
                        kappa = 0.8,
                        score_tol = 1e-2,
                        inner_tol = 1,
                        constraint_tol = 1e-3, ## changed dec 11
                        c1 = 1e-4,
                        maxit = 1000,
                        inner_maxit = 25,
                        n_final_inner_it = 25,
                        verbose = TRUE,
                        trackB = FALSE,
                        I = I,
                        Dy = Dy
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
                               p = p,
                               check_influence = TRUE)

  file_loc <- paste("local-results-postunderflow/test_",j,".rds",sep = "",collapse =)
  saveRDS(list(null_fit = null_fit,
               score_stat = score_stat),
          file_loc)


}



# # run_jth_score_test_from_init(617)
# system.time({run_jth_score_test_from_init(620)})


# run_jth_score_test(844)
#
# lapply(840:700,
#        run_jth_score_test)

# choices <- sample(c(617, 620, 623, 626, 629, 632, 634,
#        635, 637, 638, 640, 641, 643, 644, 646, 647, 649, 650, 652, 653,
#        655, 656, 658, 659, 661, 662, 664, 665, 667, 668, 670, 671, 673,
#        674, 676, 677, 679, 680, 682, 683, 685, 686, 688, 689, 691, 692,
#        694, 695, 697, 698, 700, 701, 703, 704, 706, 707, 709, 710, 712,
#        713, 715, 716, 718, 719, 721, 722, 724, 725, 727, 728, 730, 731,
#        733, 734, 736, 737, 739, 740, 742, 743, 745, 746, 748, 749, 751,
#        752, 754, 755, 757, 758, 760, 761, 763, 764, 766, 767, 769, 770,
#        772, 773, 775, 776, 778, 779, 781, 782, 784, 785, 787, 788, 790,
#        791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803,
#        804, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816,
#        817, 818, 819, 820, 821, 822, 823, 824, 825, 826, 827, 828, 829,
#        830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 844),
#       162,
#       replace=F)
#
# parallel::mclapply(choices, run_jth_score_test_from_init, mc.cores=3)


test2_without_I_Dy <- function(j){

  null_fit <-  fit_null(B = full_model$B,
                        #B (MPLE)
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
                        kappa = 0.8,
                        # score_tol = 1e-2,
                        inner_tol = 1,
                        constraint_tol = 1e-3, ## changed dec 11
                        c1 = 1e-4,
                        maxit = 1000,
                        inner_maxit = 25,
                        n_final_inner_it = 25,
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

  file_loc <- paste("local-results-postunderflow/test_",j,".rds",sep = "",collapse =)
  saveRDS(list(null_fit = null_fit,
               score_stat = score_stat),
          file_loc)


}

# run_jth_score_test(355)
system("say done")
null_fit
score_stat

j=356

# parallel::mclapply(100:616, run_jth_score_test, mc.cores=4)
run_jth_score_test(490)
parallel::mclapply(491:616, run_jth_score_test, mc.cores=2)
