% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot_stat.R
\name{boot_stat}
\alias{boot_stat}
\title{Auxiliary function for parallel implementation of rc_estim}
\usage{
boot_stat(t_ind, T_seq, und_X, twoN, N2, L1, d, f_X, Y, X, N, limit, g1, M,
  n_f_X)
}
\arguments{
\item{t_ind}{index of the parameter t}

\item{T_seq}{discrete grid for the first variable of the Partial Fourier transform}

\item{und_X}{bound on the suport of the regressors X}

\item{twoN}{parameter useful to compute the SVD}

\item{N2}{parameter useful to compute the SVD}

\item{L1}{parameter useful to compute the SVD}

\item{d}{dimension of X}

\item{f_X}{estimated density of X}

\item{Y}{outcome data}

\item{X}{regressors data}

\item{N}{sample size}

\item{limit}{apriori on the support of the density of the random slope}

\item{g1}{evaluation grid for the estimator of the density}

\item{M}{number of points in g1.}

\item{n_f_X}{estimated supnorm of the inverse of the density of the regressors}
}
\value{
a list containing, in order:

-  output

-  res0
}
\description{
Auxiliary function that compute for fixed value of t the partial Fourier transform of the density of the random coefficients. This function is used in the parallel computation of the estimator.
}
\examples{
\dontshow{
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

# beta (output) Grid
M=100
limit =7.5
b_grid  <- seq(-limit ,limit ,length.out=M)
a = limit

up =1.5
down = -up
und_beta <- a
x2 <- b_grid
x.grid <- as.matrix(expand.grid(b_grid ,b_grid ))
# DATA generating process
d = 1
Mean_mu1 = c(-2,- 3)
Mean_mu2= c(3,  0)
Sigma= diag(2, 2)
Sigma[1,2] = 1
Sigma[2,1] = 1
limit2 = 6
set.seed(2019)
M=2
b_grid  <- seq(-limit ,limit ,length.out=M)
N <-5
beta <- runif(N,  -limit2,-limit2)
X <- as.matrix(runif(N,  -up,-up))
X_t <- cbind(matrix(1, N,1),X)
Y <-rowSums(beta*X_t)
und_X= 1
n_f_X = 1/2
f_X <-rep(1/2,length(Y))
#### High level parameters
d = dim(X)[2]
N = length(Y)
M= length(b_grid)
g1 = b_grid
limit = max(abs(b_grid))
L =15
L1 = L+1
N2 = max(L,3)
twoN = 2*N2
T_seq=c(-1.5,1.5)
out <- boot_stat(1, T_seq,und_X ,twoN,N2,L1,d,f_X,Y,X,N,limit,g1,M, n_f_X )
}

}
