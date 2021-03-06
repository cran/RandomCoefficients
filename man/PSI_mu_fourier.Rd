% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PSI_mu_fourier.R
\name{PSI_mu_fourier}
\alias{PSI_mu_fourier}
\title{Auxiliary funciton for the evaluation of the SVD of F_c on a pre-specified grid divided by the singular values to the square.}
\usage{
PSI_mu_fourier(x, c, b, Psi, res)
}
\arguments{
\item{x}{a pre-specified grid}

\item{c}{the parameter indexing the singular functions}

\item{b}{the parameter indexing the smoothness class (see Gaillac and Gautier 2018)}

\item{Psi}{the Prolate Spheroidal wave functions}

\item{res}{the resolution level for the FFT.}
}
\value{
a list containing, in order:

- ipeven

- ipodd
}
\description{
Computation is done using the Fast Fourier transform for the integration.
}
\examples{
library(orthopolynom)
library(polynom)
library(tmvtnorm)
library(ks)
library(sfsmisc)
library(snowfall)
library(fourierin)
library(rdetools)
library(statmod)
library(RCEIM)
library(robustbase)
library(VGAM)
library(RandomCoefficients)
#### Bandwidth 1
L =15
L1 = L+1
N2 = max(L,3)
twoN = 2*N2
#### Bandwidth 1
c1 = 1
K1 = max(twoN+2,30)
K = K1
c = 1
b=1
bound=1
out <- get_psi_mu(c,N2,twoN,K, L1)
Psi <- out[[1]]
mu<- out[[2]]
xseq = seq(-bound,bound, length.out=10)
resol=2^7
psix <- PSI_mu_fourier(xseq,c,b,Psi,resol)
}
