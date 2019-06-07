#' Auxiliary function that compute the Legendre quadrature of order K
#'
#'
#' Generate nodes and weights for Legendre-Gauss quadrature on [-1,1]. Note that t is a column vector and w is a
#' row vector. Also normalizes and returns the eigenvectors of J so that they are samples of the unit-norm Legendre
#' polynomials
#'
#'@param  K order of the Legendre quadraure
#'
#'@return a list containing, in order:
#'
#' - t : points of the Legendre quadrature
#'
#' -  w : weigths for the Legendre quadrature
#'
#' -  Pbarmat : the eigenvectors of J
#'
#'@examples
#'K=30
#'res2 <-  legendrequad(K)
#'
#'
legendrequad <- function(K){

  u <- sqrt(1/(4-1/seq(1,(K-1))^2))
  n = length(u)+1
  trans = myDiag(matrix(0,n, n),u,1) + myDiag(matrix(0,n, n),u,-1)
  eigen_trans <- eigen(trans)
  V<- eigen_trans$vectors
  Lambda <- eigen_trans$values
  t = sort(Lambda)
  i= sort(seq(1, length(Lambda)), decreasing=TRUE)
  V = V[,i,drop=FALSE]
  Vtop = V[1,,drop=FALSE]
  w = 2*Vtop^2
  Pbarmat = V/repmat(Vtop*sqrt(2),K,1)
  res <- vector("list")
  res[[1]] <- t
  res[[2]] <- w
  res[[3]] <- Pbarmat
  return(res)
}

