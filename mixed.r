###########
## @Description           a simple normal mixture modeling function for clustering to the subfamily EII
## @data                  dataframe of n observations and p variables  
## @k                     number of initial clusters 
## @return                returns a list of cluster weights, means, and covariance matrices
############  

mixed <- function(data,k){
    n <- nrow(data)
    p <- ncol(data)
    weights <- rep(1/k, k)
    means <- matrix(runif(p*k, min = min(data), max = max(data)), nrow = p, ncol = k)
    covar <- rep(diag(p) * sd(data) / sqrt(k), k)
    
  for(i in 1:n){
    for(j in 1:k){
      probs[i,j] <- weights[j] * dmvnorm(data[i,], means[,j], covar[,,j], log = FALSE)
    }
    probs[i,] <- probs[i,] / sum(probs[i,])
    }
  weights.new <- rowSums(probs) / n
  means.new <- t(apply(probs, 2, function(x) t(data) %*% x / sum(x)))
  covar.new <- array(0, dim = c(p, p, k))
  
  for(j in 1:k){
    covar.new[,,j] <- (t(data - means.new[,j]) %*% (probs[,j] * (data - means.new[,j]))) / sum(probs[,j])
    }
  covar.new <- array(covar.new, c(p, p, k))
  weights <- weights.new
  means <- means.new
  covar <- covar.new
  
  return(list(weights = weights, means = means, covariances = covar))
}
