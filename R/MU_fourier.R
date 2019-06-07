#' Auxiliary function that computes the singular values of the SVD of the operator F_c in Gaillac and Gautier 2018 using the Fast Fourier transform for the integration.
#'
#'@param  psi Prolate spheroidal wave functions
#'@param xseq grid on which to evaluate them, output of the Legendre quadrature
#'@param splin use interpolation by splines or not (boolean).
#'
#'@return mu
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
#'c = 1
#'b=1
#'bound=1
#'out <- get_psi_mu(c,N2,twoN,K, L1)
#'Psi <- out[[1]]
#'mu<- out[[2]]
#'xseq = seq(-bound,bound, length.out=10)
#'resol=2^7
#'psix <- PSI_mu_fourier(xseq,c,b,Psi,resol)
#'psi <- psix[[1]]
#'xseq <- psix[[2]]
#'xseq2 <- xseq/b
#'splin =FALSE
#'mu_ev<- MU_fourier(psi,xseq,splin)
#'mu <- mu_ev[[1]]



MU_fourier <- function(psi,xseq,splin){
  ### integration over R approximated on large compact support
  mu <- NULL
  for(j in seq(1,dim(psi)[1])){
    ip = integrate.xy(xseq, abs(psi[j,])*abs(psi[j,]),use.spline=splin,xtol=10e-15)
    # ip = sum(abs(psi[j,])*abs(psi[j,])*wseq)
    mu<- rbind(mu,(1i)^(2*(j-1))*sqrt(ip))
    # mu_even<- rbind(mu_even, (1i)^(2*(n-1))*sqrt(abs(ipeven[n])))
  }
  res <- vector("list")
  res[[1]] <- (mu)
  return(res)
}
