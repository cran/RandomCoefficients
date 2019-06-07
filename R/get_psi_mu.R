#' Auxiliary function for the computation of the eigenvalues of the SVD of F_c
#'
#' This function compute the  eigenvalues of the SVD of F_c.
#'
#' @param c1 parameter indexing the SVC
#' @param N2 maximal number of elements of the SVD that are computed.
#' @param twoN number of  Legendre polynomials that are loaded
#' @param K1 order of the Legendre quadrature
#' @param L1 number of Legendre polynomials used in the computation
#'
#' @return a list containing, in order:
#'
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
#'### get the beta's (for the computation of Psi)
#'out <- get_psi_mu(c1,N2,twoN,K1, L1)
#'Psi1 <- out[[1]]
#'mu1<- out[[2]]

get_psi_mu <- function(c1,N2,twoN,K1, L1){

  res11 <- get_u(c1,N2)
  ueven1 <- res11[[1]]
  uodd1 <- res11[[2]]

  ###############################################

  leg1 <- legendre.polynomials(twoN+1, normalized=TRUE)
  Psi_even1 <- vector("list")
  Psi_odd1 <- vector("list")
  for(n in seq(1,L1)){
    poly1 <- NULL
    for(i in seq(0,N2)){
      poly1 <- poly1 + ueven1[i+1,n]*leg1[[2*i+1]]
    }
    Psi_even1 [[n]] <- poly1
  }

  for(n in seq(1,L1)){
    poly1 <- NULL
    for(i in seq(1,N2+1)){
      poly1 <- poly1 + uodd1[i,n]*leg1[[2*i]]
    }
    Psi_odd1[[n]] <- poly1
  }

  #### allows to sort psis
  output <- vector("list")
  output[[1]] <- insertEO(Psi_odd1, Psi_even1)
  ##### compute the eigenvalues mu, K must be sufficiently large
  mu1<- MU(N2,c1,K1, L1)
  # mu_even1 <-mu1[[1]]
  # mu_odd1 <- mu1[[2]]
  output[[2]] <- insertEO(mu1[[2]], mu1[[1]])
  return(output)
}
