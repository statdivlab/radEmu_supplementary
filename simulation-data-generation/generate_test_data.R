generate_test_data <- function(n, # number of samples
                               J, # number of taxa
                               b0 = NULL, # intercept (optional, required if B is null)
                               b1 = NULL, # beta_1 vector (optional, required if B is null)
                               distn, # distribution: Poisson or ZINB
                               zinb_size = NULL, # size argument for negative binomial distribution for ZINB data
                               zinb_zero_prop = NULL, # proportion zeros in ZINB data
                               mean_count_before_ZI, # expected z parameter for sample with mean expected count of 0
                               X = NULL, # design matrix (optional)
                               B = NULL, # B matrix (optional, required if b0 and b1 are NULL)
                               cluster = NULL ) { # cluster vector (optional)
  
  if (is.null(X)) {
    X <- cbind(1,rep(c(0,1),each = n/2))
  }
  if (!is.null(b0) & !is.null(b1)) {
    B <- rbind(b0,b1)
    if (nrow(B) != ncol(X)) {
      stop("You've input b0 and b1 but your X matrix does not have 2 columns. Please use the B argument when your design matrix does not have 2 columns.")
    }
  } else if (is.null(b0) | is.null(b1)) {
    if (is.null(B)) {
      stop("Please input either parameter vectors b0 and b1, or parameter matrix B.")
    }
  }
  log_means <- do.call(cbind,
                       lapply(1:J,
                              function(j) X%*%B[,j,drop = FALSE]))
  
  row_means <- rowSums(exp(log_means))/J
  
  z <- sapply(row_means,function(x) log(mean_count_before_ZI) - log(x) + stats::rnorm(1))
  Y <- matrix(0, ncol = J, nrow = n)
  
  for(i in 1:n){
    log_means[i,] <- log_means[i,] + z[i]
  }
  
  if (is.null(cluster)) {
    for(i in 1:n){
      accepted <- FALSE
      while(!accepted){
        for(j in 1:J){
          if(distn == "Poisson"){
            Y[i,j] <- stats::rpois(1,lambda = exp(log_means[i,j]))
          }
          if(distn == "ZINB"){
            Y[i,j] <- stats::rnbinom(1,mu = exp(log_means[i,j]), size= zinb_size)*
              (1 - stats::rbinom(1,1,prob = zinb_zero_prop))
          }
          if(sum(Y[i,])>0){
            accepted <- TRUE
          }
        }
      }
    }
  } else {
    cluster_effs <- lapply(1:length(unique(cluster)), function(i) log(matrix(stats::rexp(2*J), nrow = 2)))
    for(i in 1:n){
      accepted <- FALSE
      while(!accepted){
        for(j in 1:J){
          if(distn == "Poisson"){
            Y[i,j] <- stats::rpois(1,lambda = exp(log_means[i,j] + 
                                                    cluster_effs[[cluster[i]]][, j]))
          }
          if(distn == "ZINB"){
            Y[i,j] <- stats::rnbinom(1, 
                                     mu = exp(log_means[i,j] + 
                                                cluster_effs[[cluster[i]]][, j]), 
                                     size= zinb_size)*
              (1 - stats::rbinom(1,1,prob = zinb_zero_prop))
          }
          if(sum(Y[i,])>0){
            accepted <- TRUE
          }
        }
      }
    }
  }
  
  return(Y)
}

generate_test_data_include_delta <- function(n, # number of samples
                                             J, # number of taxa
                                             distn, # distribution: Poisson or ZINB
                                             zinb_size = NULL, # size argument for negative binomial distribution for ZINB data
                                             zinb_zero_prop = NULL, # proportion zeros in ZINB data
                                             mean_count_before_ZI, # expected z parameter for sample with mean expected count of 0
                                             X, # design matrix
                                             B, # B matrix 
                                             cluster = NULL ) { # cluster vector (optional)
  
  
  delta_j <- seq(from = -1, to = 1, length.out = J)
  perm_delta_j <- delta_j
  ind_keep <- c(ceiling(J / 4), J / 2, ceiling(3 * J / 4))
  perm_delta_j[-ind_keep] <- sample(delta_j[-ind_keep], J - 3, replace = FALSE)
  
  log_means <- do.call(cbind,
                       lapply(1:J,
                              function(j) X%*%B[,j,drop = FALSE] + perm_delta_j[j]))
  
  row_means <- rowSums(exp(log_means))/J
  
  z <- sapply(row_means,function(x) log(mean_count_before_ZI) - log(x) + stats::rnorm(1))
  Y <- matrix(0, ncol = J, nrow = n)
  
  for(i in 1:n){
    log_means[i,] <- log_means[i,] + z[i]
  }
  
  if (is.null(cluster)) {
    for(i in 1:n){
      accepted <- FALSE
      while(!accepted){
        for(j in 1:J){
          if(distn == "Poisson"){
            Y[i,j] <- stats::rpois(1,lambda = exp(log_means[i,j]))
          }
          if(distn == "ZINB"){
            Y[i,j] <- stats::rnbinom(1,mu = exp(log_means[i,j]), size= zinb_size)*
              (1 - stats::rbinom(1,1,prob = zinb_zero_prop))
          }
          if(sum(Y[i,])>0){
            accepted <- TRUE
          }
        }
      }
    }
  } else {
    cluster_effs <- lapply(1:length(unique(cluster)), function(i) log(matrix(stats::rexp(2*J), nrow = 2)))
    for(i in 1:n){
      accepted <- FALSE
      while(!accepted){
        for(j in 1:J){
          if(distn == "Poisson"){
            Y[i,j] <- stats::rpois(1,lambda = exp(log_means[i,j] + 
                                                    cluster_effs[[cluster[i]]][, j]))
          }
          if(distn == "ZINB"){
            Y[i,j] <- stats::rnbinom(1, 
                                     mu = exp(log_means[i,j] + 
                                                cluster_effs[[cluster[i]]][, j]), 
                                     size= zinb_size)*
              (1 - stats::rbinom(1,1,prob = zinb_zero_prop))
          }
          if(sum(Y[i,])>0){
            accepted <- TRUE
          }
        }
      }
    }
  }
  
  return(list(Y = Y, z = z))
}