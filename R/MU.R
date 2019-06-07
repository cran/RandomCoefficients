#' Auxiliary function that computes the singular values of the SVD of the operator F_c in Gaillac and Gautier 2018.
#'
#' @param N number of singular values to compute
#' @param c parameter indexing the singular values
#' @param K ordre of the Legendre quadrature
#' @param L1 number of Legendre polynomials used in the computation of the coefficients of the singular functions;
#'
#' @return a list containing, in order:
#'
#' - mu_even
#'
#' - mu_odd
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
#'c1 = 1
#'K1 = max(twoN+2,30)
#'K = K1
#'mu1<- MU(N2,c1,K1, L1)


MU <- function(N,c,K, L1){
  res2 <-  legendrequad(K)
  t <- res2[[1]]
  w<- res2[[2]]
  psi_even <- abs(PSI_mu(t,N,c,K, L1)[[1]])
  psi_odd<- abs(PSI_mu(t,N,c,K, L1)[[2]])
  mu_even <- NULL
  for(n in seq(1,L1)){
    q = t(w)*psi_even[n,]
    ipeven = psi_even%*%q
    mu_even<- rbind(mu_even, (1i)^(2*(n-1))*sqrt(ipeven[n]))
    #   mu_even<- rbind(mu_even, (1i)^(2*(n-1))*sqrt(abs(ipeven[n])))
  }

  mu_odd <- NULL
  for(n in seq(1,L1)){
    q = t(w)*psi_odd[n,]
    ipodd= psi_odd%*%q
    mu_odd<- rbind(mu_odd, (1i)^(2*(n-1)+1)*sqrt(ipodd[n]))
  }
  res <- vector("list")
  res[[1]] <- (mu_even)
  res[[2]] <- (mu_odd)
  return(res)
}

