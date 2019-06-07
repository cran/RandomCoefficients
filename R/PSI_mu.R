#' Ausiliary function  that evaluates the SVD of F_c on a pre-specified grid divided by the singular values to the square.
#'
#' @param  x the pre-specified grid
#' @param N number of singular values to compute
#' @param c parameter indexing the singular values
#' @param K ordre of the Legendre quadrature
#' @param L1 number of Legendre polynomials used in the computation of the coefficients of the singular functions;
#'
#'
#' @return a list containing, in order:
#'
#' - ipeven
#'
#' - ipodd
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
#'#### Bandwidth 1
#'L =15
#'L1 = L+1
#'N2 = max(L,3)
#'twoN = 2*N2
#'#### Bandwidth 1
#'c= 1
#'K1 = max(twoN+2,30)
#'K = K1
#'res2 <-  legendrequad(K)
#'t <- res2[[1]]
#'psi_even <- abs(PSI_mu(t,N2,c,K, L1)[[1]])
#'

PSI_mu <- function(x,N,c,K, L1){
  twoN = 2*N
  ## get the beta's
  res1 <- get_u(c,N)
  ueven <- res1[[1]]
  uodd <- res1[[2]]

  ## get the legendre quadrature
  res2 <-  legendrequad(K)
  t <- res2[[1]]
  w<- res2[[2]]
  Pbarmat<- res2[[3]]
  Pbeven = Pbarmat[seq(1,twoN+1,by=2),]
  Pbodd = Pbarmat[seq(2,twoN+2,by=2),]


  ### compute the scalar product exp(ictx)* psi_n(x)
  psieven = t(ueven[,1:L1])%*%Pbeven
  psiodd = t(uodd[,1:L1])%*%Pbodd
  q = (t(matrix(w, 1, length(w)))%*%matrix(1,1, length(x)))*exp(1i*c*t(matrix(1,1, length(t)))%*%matrix(x,1, length(x))*t)
  ipeven = psieven%*%q
  ipodd = psiodd%*%q


  res <- vector("list")
  res[[1]] <- (ipeven)
  res[[2]] <- (ipodd)
  return(res)
}
