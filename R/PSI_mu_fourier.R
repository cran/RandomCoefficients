
#' Auxiliary funciton for the evaluation of the SVD of F_c on a pre-specified grid divided by the singular values to the square.
#'
#' Computation is done using the Fast Fourier transform for the integration.
#'
#'@param  x a pre-specified grid
#'@param  c the parameter indexing the singular functions
#'@param  b the parameter indexing the smoothness class (see Gaillac and Gautier 2018)
#'@param  Psi the Prolate Spheroidal wave functions
#'@param  res the resolution level for the FFT.
#'
#'
#'@return a list containing, in order:
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



PSI_mu_fourier  <- function(x,c,b, Psi, res){

  ipev <- NULL
  for(i in 1:length(Psi)){

    pst <- function(t){
      out <- matrix(0,1,length(t))
      out[abs(t)>1] <- 0
      out[abs(t)<=1] <- as.function(Psi[[i]])(t)
    }

    pst0 <- fourierin_1d(  pst, -1, 1,lower_eval = min(x), upper_eval = max(x),freq_adj=-c/b,const_adj=1,resolution =res)
    ipev <- rbind(ipev, pst0$values*((pst0$w <= 1)& (pst0$w >=- 1)))
    # ipev <- rbind(ipev, pst0$values)
  }

  res <- vector("list")
  res[[1]] <- (ipev)
  res[[2]] <- pst0$w
  return(res)
}
