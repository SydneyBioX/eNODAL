#' Calculate LC-Stats
#' @param D Experimental matrix
#' @param M A matrix indicate which axis used in D to calculate LC
#' @param Z A vector of gene of interest.
#' @param lmd Scale factor to calculate weighted matrix in LC-Stats.
#' @importFrom stats var dist
#' @return A number of LC-Stats
lc_stat = function(D, M, Z, lmd = NULL){

  eps = 0.001

  if (is.vector(Z)){
    Z = matrix(Z,ncol = 1)
  }else{
    Z = as.matrix(Z)
  }

  D = scale(D)
  #  Z = scale(Z)

  D = as.matrix(D)


  if(is.null(lmd)){
    adp_weight = 1/sqrt(2)*(apply(Z,2,var) + eps)
    lmd = adp_weight
  }

  #  print(Z)
  D_z  = abs(outer(Z[,1],Z[,1],"-"))


  D_M = exp(as.matrix(-dist(D %*% M)^2))
  # print(D_M)
  # print(D_z)

  R = sum(D_M * abs(D_z))
  return(R)
}

#' Calculate LC-Test
#' @param D Experimental matrix
#' @param Z A vector of gene of interest.
#' @param n Number of permutation
#' @return Pvalue of LC-test
LC_test0 <- function(D, Z, n = 100){


  M <- diag(rep(1,ncol(D)))
  D = scale(D)
  Z = scale(Z)
  lc_0 = lc_stat(D,M,Z)
  lc = c()
  for(i in 1:n){
    Z1 = sample(Z,replace = F)
    lc[i] = lc_stat(D, M, Z1)
  }

  p = sum(lc<lc_0)/n
  return(p)
}

