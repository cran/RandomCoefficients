#' Auxiliary function that extends the matrix X
#'
#'@param  X A vector of inputs
#'@param  m the first dimension of the desired output matrix
#'@param  n the second dimension of the desired output matrix
#'
#'
#'@return A matrix os size (dim(X)[1]*m,dim(X)[2]*n)
#'
#'@examples
#'library(orthopolynom)
#'library(polynom)
#'library(tmvtnorm)
#'library(ks)
#'library(sfsmisc)
#'library(snowfall)
#'library(fourierin)
#'library(rdetools)
#'library(statmod)
#'library(RCEIM)
#'library(robustbase)
#'library(VGAM)
#'library(RandomCoefficients)
#'K=3
#'u <- sqrt(1/(4-1/seq(1,(K-1))^2))
#'n = length(u)+1
#'trans = myDiag(matrix(0,n, n),u,1) + myDiag(matrix(0,n, n),u,-1)
#'eigen_trans <- eigen(trans)
#'V<- eigen_trans$vectors
#'Lambda <- eigen_trans$values
#'t = sort(Lambda)
#'i= sort(seq(1, length(Lambda)), decreasing=TRUE)
#'V = V[,i,drop=FALSE]
#'Vtop = V[1,,drop=FALSE]
#'w = 2*Vtop^2
#'Pbarmat = V/repmat(Vtop*sqrt(2),K,1)
#'
repmat = function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  return(matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T))
  }
