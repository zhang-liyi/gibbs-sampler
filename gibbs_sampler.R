library(tidyverse)
library(parallel)
options(mc.cores = parallel::detectCores())

# =================================================================+
# Get Data & Clean Data
# =================================================================+

data <- read_delim('ap/ap.dat', '\t', col_names=FALSE)
data <- data$X1
N <- length(data)

vocab <- read_delim('ap/vocab.txt', '\t', col_names=FALSE)
vocab <- vocab$X1
V <- length(vocab)

X <- array(0, c(N, V))

for(i in 1:N){
  article_str <- strsplit(data[i], ' ')[[1]]
  for(j in 2:length(article_str)){
    word_count <- as.integer(strsplit(article_str[j], ':')[[1]])
    X[i, word_count[1]+1] <- word_count[2]
  }
}

# Remove words that occur only in at most one article
idx_remove <- c()
for(i in 1:V){
  word_found <- 0
  for(n in 1:N){
    if(X[n, i] > 0){
      word_found <- word_found + 1
    }
    if(word_found > 1){
      break
    }
  }
  if(word_found <= 1){
    idx_remove <- c(idx_remove, -i)
  }
}

if(length(idx_remove)>0){
  X <- X[,idx_remove]
}
vocab_trim <- vocab[idx_remove]

# Divide out held-out data
idx_heldout <- sample(1:N, as.integer(N/10))
X_heldout <- X[idx_heldout,]
X <- X[-idx_heldout,]
N <- dim(X)[1]
V <- dim(X)[2]

X_listrep <- list()
for(n in 1:N){
  X_listrep[[n]] <- X[n,]
}


# Function to save data; data is intended to be the list outputed by gibbs_sampler()
quick_save <- function(data){
  t <- strsplit(as.character(Sys.time()), ' ')[[1]]
  date <- t[1]
  time <- paste(strsplit(t[2], ':')[[1]], collapse='')
  time_str <- paste(c(date, time), collapse='-')
  dir.create(file.path('results', time_str), recursive=TRUE)
  
  names <- list('logjoint.Rda', 'lambda.Rda', 'theta.Rda', 'z.Rda')
  for(i in 1:length(data)){
    data_i <- data[[i]]
    name <- names[[i]]
    saveRDS(data_i, file=file.path('results', time_str, name))
  }
  return(time_str)
}

# =================================================================+
# Model-building and Inference
# =================================================================+

# Computes Dirichlet density; parameter a must be a vector
ddirichlet <- function(x, a){
  dsimplex <- prod(x^(a-1))
  denom <- prod(gamma(a))/gamma(sum(a))
  return(dsimplex/denom)
}

# Helper for compute_log_joint
compute_logjoint_helper_z <- function(z_row, theta){
  return(dmultinom(z_row, prob=theta, log=TRUE))
}

# Helper for compute_log_joint
compute_logjoint_helper_x <- function(x_row, z_row, lambda){
  lambda_k <- lambda[which(z_row %in% 1),]
  return(sum(dpois(x_row, lambda_k, log=TRUE)))
}

# Computes log-joint of this mixture model (parallelized)
compute_logjoint <- function(x, z, lambda, theta, alpha, a, b){
  N <- dim(x)[1]
  logp_theta <- log(ddirichlet(theta, alpha))
  logp_lambda <- sum(dgamma(lambda, a, b, log=TRUE))
  logp_z <- sum(unlist(mclapply(z, compute_logjoint_helper_z, theta, mc.silent=TRUE)))
  logp_x <- sum(unlist(mcmapply(compute_logjoint_helper_x, x, z, MoreArgs=list(lambda), mc.silent=TRUE)))
  return(logp_theta + logp_lambda + logp_z + logp_x)
}

# Computes complete-conditional of z (for parallelization)
comp_cond_z <- function(X_row, lambda, theta, K, V){
  loglik <- dpois(matrix(rep(X_row, K), K, V, byrow=TRUE), 
                  lambda,
                  log=TRUE)
  loglik <- rowSums(loglik)
  log_unnorm_param <- log(theta) + loglik
  log_unnorm_param <- log_unnorm_param - max(log_unnorm_param)
  z_row_new <- t(rmultinom(1, 1, exp(log_unnorm_param)))
  return(z_row_new)
}

# Convert data strucutre
z_to_classmap <- function(z_row){
  class <- which(z_row %in% 1)
  return(class)
}

# Main function
gibbs_sampler <- function(X, K, num_iter=1000){
  N <- dim(X)[1]
  V <- dim(X)[2]
  
  # Set hyperparameters:
  alpha <- rep(1, K)
  a <- 0.1
  b <- 0.1
  
  # Initialize parameters:
  lambda_0_perturb <- matrix(rnorm(K*V, 0, 1), K, V)
  lambda_0 <- matrix(colMeans(X), K, V, byrow=TRUE) # dim KxV; prior Gamma(a, b)
  lambda_0 <- pmax(lambda_0 + lambda_0 * lambda_0_perturb * 0.5, matrix(0.02, K, V))
  theta_0 <- array(1/K, K) # dim 1xK; prior Dirichlet(alpha)
  z_0 <- t(rmultinom(N, 1, theta_0)) # dim NxK; Categorical(theta)
  
  # Prepare data structures for parameters and log-joint:
  log_joint <- array(0, num_iter+1)
  lambda <- array(0, c(num_iter+1, K, V))
  theta <- array(0, c(num_iter+1, K))
  z <- array(0, c(num_iter+1, N, K))
  lambda[1,,] <- lambda_0
  theta[1,] <- theta_0
  z[1,,] <- z_0
  z_0_listrep <- as.list(as.data.frame(t(z_0)))
  log_joint[1] <- compute_logjoint(X_listrep, z_0_listrep, lambda_0, theta_0, alpha, a, b)
  print(paste('Initial log-joint:', as.character(log_joint[1]), collapse=' '))
  
  # Auxilliary data structure for computation:
  class_map <- c()
  for(n in 1:N){
    k <- which(z[1,n,] %in% 1)
    class_map <- c(class_map, k)
  }
  
  # Gibbs Sampling begins:
  for(i in 2:(num_iter+1)){
    # Complete conditional of proportions theta
    param_a <- colSums(z[i-1,,]) + alpha
    samples <- rgamma(K, param_a, 1)
    theta[i,] <- samples/sum(samples)
    
    # Complete conditional of component-param lambda
    for(k in 1:K){
      indices <- which(class_map %in% k)
      X_inclass <- X[indices,]
      lambda[i,k,] <- rgamma(V, a+colSums(X_inclass), b+length(indices))
    }
    
    # Complete conditional of assignments z 
    z_list <- mclapply(X_listrep, comp_cond_z, lambda[i,,], theta[i,], K, V, mc.silent=TRUE)
    z[i,,] <- do.call(rbind, z_list)
    class_map <- unlist(mclapply(z_list, z_to_classmap, mc.silent=TRUE))
    
    log_joint[i] <- compute_logjoint(X_listrep, z_list, lambda[i,,], theta[i,], alpha, a, b)
    print(paste('Iteration:', as.character(i-1), '      ',
                'Log-joint:', as.character(log_joint[i]), 
                collapse=' '))
  }
  
  output <- list('logjoint'=log_joint, 'lambda'=lambda, 'theta'=theta, 'z'=z)
  quick_save(output)
  return(output)
}


# Run the program:
# - It's running two chains. 
# - Just add more X into the list (e.g. list(X,X,X,X)) to run more chains.
K <- 5
output <- mclapply(list(X,X), gibbs_sampler, K=K, num_iter=1000, mc.silent=TRUE)

# =================================================================+
# Model Analysis
# =================================================================+

# Compute autocorrelation (within sequence correlation):
autocor <- function(x, t){
  n <- length(x)
  mu <- mean(x)
  sigma2 <- var(x)
  numer <- sum((x[1:(n-t)] - mu) * (x[(t+1):n] - mu))
  return(numer/((n-t)*sigma2))
}

# Compute number of effective sample-size (Neff):
neff <- function(x){
  n <- length(x)
  autocor_arr <- array(0, n-1)
  for(t in 1:(n-1)){
    autocor_arr[t] <- autocor(x, t)
    if(t > 1){
      if(autocor_arr[t] + autocor_arr[t-1] < 0){
        autocor_arr[t] <- 0
        autocor_arr[t-1] <- 0
        break
      }
    }
  }
  return(n/(1+2*sum(autocor_arr)))
}

lambda_c1 <- readRDS('results/2020-10-12-225004/lambda.Rda')
lambda_c2 <- readRDS('results/2020-10-12-225016/lambda.Rda')
logjoint_c1 <- readRDS('results/2020-10-12-225004/logjoint.Rda')
logjoint_c2 <- readRDS('results/2020-10-12-225016/logjoint.Rda')
X_col_mean <- colMeans(X)
w <- 1 # which word's corresponding param do we want to analyze


# ---------- Within-chain Analysis --------------------------------+

# Plot log-joints of each chain
plot(logjoint_c1, type='l', col='red', xlab='Iteration', ylab='Log-joint')
lines(logjoint_c2, type='l', col='blue')

# Autocorrelation
param_c2 <- apply(lambda_c2[501:1001,,w], c(1), mean)
for(t in 1:10){
  if(t==1){print(paste('Autocorrelation for word index', as.character(w), collapse=' '))}
  print(autocor(param_c2, t))
}

# Neff
idx_samp <- sample(1:V, 10)
for(i in 1:10){
  w <- idx_samp[i]
  if(i==1){print(paste('Neff'))}
  #param_c2 <- apply(lambda_c2[501:1001,,w], c(1), mean)
  param_c2 <- lambda_c2[501:1001,1,w]
  print(neff(param_c2))
}


# ---------- Across-chain Analysis --------------------------------+

w <- 1 # which word's corresponding param do we want to analyze
for(k in 1:K){
  if(k==1){plot(lambda_c1[,k,w], type='l', col='red', 
                ylim=c(0,4*X_col_mean[w]),
                ylab='Sampled parameter',
                xlab='Iteration')}
  else{lines(lambda_c1[,k,w], type='l', col='red')}
  lines(lambda_c2[,k,w], type='l', col='blue')
}


# ---------- Visualize Topics -------------------------------------+

# Computes term-score
term_score <- function(lambda, i1, i2){
  K <- dim(lambda)[2]
  scores <- array(0, c(K,V))
  lambda_hat <- apply(lambda[i1:i2,,], c(2,3), mean)
  for(k in 1:K){
    for(v in 1:V){
      denom <- prod(lambda_hat[,v])**(1/K)
      scores[k,v] <- lambda_hat[k,v] * log(lambda_hat[k,v] / denom)
    }
  }
  return(scores)
}

scores <- term_score(lambda_c2, 501, 1001)

for(k in 1:K){
  scores_k <- scores[k,]
  idx_vocab <- rev(order(scores_k))[1:20]
  print(vocab_trim[idx_vocab])
  print('------------------------------')
}








