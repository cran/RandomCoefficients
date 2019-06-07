#' Computation of the coefficients of the PSWF on the Legendre polynomials basis of L^2(-1,1)
#'
#'  The matrix ueven has columns u^0,u^2,...,u^{2N}, and the matrix uodd has columns u^1,u^3,...,u^{2N+1}.
#' Each column has length N+1. Below, note that i must have odd length so that ai and ci have even length and can be split into two subvectors containing alternate entries of ai and ci.
#'
#' @param  c parameter indexing the functions of the SVD
#' @param  N maximal index of the functions computed in the SVD
#'
#' @return a list containing, in order:
#'
#' - ueven : the coefficients of the decomposition of the SVD on the Legendre polynomials for the even functions
#'
#' - uodd : the coefficients of the decomposition of the SVD on the Legendre polynomials for the odd functions
#' @examples
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
#'#### Number of Psis
#'L =15
#'L1 = L+1
#'N2 = max(L,3)
#'twoN = 2*N2
#'#### Bandwidth 1
#'c1 = 1
#'K1 = max(twoN+2,30)
#'K = K1
#'### get the coefficients beta (for the computation of Psi) of the PSWF on the
#'### basis of Legendre polynomials.
#'out <- get_u(c1,N2)
#'
get_u <- function(c,N){
  #
  i = seq(0,2*N+2)
  i1 = i +1
  twoi = 2*i
  twoiP3 = twoi+3
  fi = i1/twoiP3
  ai = fi[1:(length(fi)-1)]*fi[2:length(fi)]
  aeven = ai[seq(1,(length(ai)-1), by=2)]
  aodd = ai[seq(2,(length(ai)), by=2)]
  gi = (i-1)/(twoi-3)
  ci = gi[1:(length(gi)-1)]*gi[2:length(gi)]
  ceven = ci[seq(1,length(ci)-1, by=2)]
  codd = ci[seq(2,length(ci), by=2)]
  c2 = c*c
  ii1 = i*i1
  bi = ii1+((2*ii1-1)/((twoi-1)*(twoiP3)))*c2
  be = bi[seq(1,length(bi)-1, by=2)]
  bo = bi[seq(2,length(bi), by=2)]


  e = c2*sqrt(aeven[1:(length(aeven)-1)]*ceven[2:length(ceven)])

  trans = diag(be) + myDiag(matrix(0,nrow(diag(be)), ncol(diag(be))),e,1) + myDiag(matrix(0,nrow(diag(be)), ncol(diag(be))),e,-1)


  eigen_trans <- eigen(trans)
  ueven <- eigen_trans$vectors
  le<- eigen_trans$values

  chi <- sort(le)
  i <- sort(seq(1, length(chi)),decreasing = TRUE)

  ueven = ueven[,i,drop=FALSE]
  e = c2*sqrt(aodd[1:(length(aodd)-1)]*codd[2:length(codd)])

  trans2 = diag(bo) + myDiag(matrix(0,nrow(diag(bo)), ncol(diag(bo))),e,1) + myDiag(matrix(0,nrow(diag(bo)), ncol(diag(bo))),e,-1)
  eigen_trans2 <- eigen(trans2)
  uodd <- eigen_trans2$vectors
  lo<- eigen_trans2$values

  chi <- sort(lo)
  i <- sort(seq(1, length(chi)),decreasing = TRUE)

  uodd = uodd[,i,drop=FALSE]
  res <- vector("list")
  res[[1]] <- ueven
  res[[2]] <- uodd
  return(res)
}
