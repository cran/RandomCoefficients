#'  Auxiliary function that put together even and odd functions of the SVD in an only one output list.
#'
#'@param Psi_odd odd singular functions Psi
#'@param Psi_even even singular functions Psi
#'
#'@return Psi
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
#'mu1<- MU(N2,c1,K1, L1)
#'output <- insertEO(mu1[[2]], mu1[[1]])
#'


insertEO <- function( Psi_odd, Psi_even){
  Psi <- vector("list")
  for(i in seq(1,length(Psi_even)) ){
    Psi[[2*(i-1)+1]] <- Psi_even[[i]]
    Psi[[2*i]] <- Psi_odd[[i]] }
  return(Psi)
}
