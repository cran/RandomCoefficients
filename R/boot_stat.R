#' Auxiliary function for parallel implementation of rc_estim
#'
#' Auxiliary function that compute for fixed value of t the partial Fourier transform of the density of the random coefficients. This function is used in the parallel computation of the estimator.
#'
#' @param t_ind index of the parameter t
#' @param T_seq discrete grid for the first variable of the Partial Fourier transform
#' @param und_X bound on the suport of the regressors X
#' @param twoN parameter useful to compute the SVD
#' @param N2 parameter useful to compute the SVD
#' @param L1 parameter useful to compute the SVD
#' @param d dimension of X
#' @param f_X estimated density of X
#' @param Y outcome data
#' @param X regressors data
#' @param N sample size
#' @param limit apriori on the support of the density of the random slope
#' @param g1 evaluation grid for the estimator of the density
#' @param M number of points in g1.
#' @param n_f_X estimated supnorm of the inverse of the density of the regressors
#'
#' @return a list containing, in order:
#'
#' -  output
#'
#' -  res0
#'@examples
#'\dontshow{
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
#'
#'# beta (output) Grid
#'M=100
#'limit =7.5
#'b_grid  <- seq(-limit ,limit ,length.out=M)
#'a = limit
#'
### Support apriori limits (taken symmetric)
#'up =1.5
#'down = -up
#'und_beta <- a
#'x2 <- b_grid
#'x.grid <- as.matrix(expand.grid(b_grid ,b_grid ))
#'# DATA generating process
#'d = 1
#'Mean_mu1 = c(-2,- 3)
#'Mean_mu2= c(3,  0)
#'Sigma= diag(2, 2)
#'Sigma[1,2] = 1
#'Sigma[2,1] = 1
#'limit2 = 6
#'set.seed(2019)
#'M=2
#'b_grid  <- seq(-limit ,limit ,length.out=M)
#'N <-5
#'beta <- runif(N,  -limit2,-limit2)
#'X <- as.matrix(runif(N,  -up,-up))
#'X_t <- cbind(matrix(1, N,1),X)
#'Y <-rowSums(beta*X_t)
#'und_X= 1
#'n_f_X = 1/2
#'f_X <-rep(1/2,length(Y))
#'#### High level parameters
#'d = dim(X)[2]
#'N = length(Y)
#'M= length(b_grid)
#'g1 = b_grid
#'limit = max(abs(b_grid))
#'L =15
#'L1 = L+1
#'N2 = max(L,3)
#'twoN = 2*N2
#'T_seq=c(-1.5,1.5)
#'out <- boot_stat(1, T_seq,und_X ,twoN,N2,L1,d,f_X,Y,X,N,limit,g1,M, n_f_X )
#'}
#'
boot_stat <- function(t_ind, T_seq,und_X ,twoN,N2,L1,d,f_X,Y,X,N,limit,g1,M, n_f_X ){
  t=T_seq[t_ind]
  #       cat(paste0( "  t   : ", t ," \n"))
  #### define the adapted bandwidth
  c = 1.03*limit* und_X * abs(t)
  K = max(twoN+2,100)

  ####
  out <- get_psi_mu(c,N2,twoN,K, L1)
  Psi <- out[[1]]
  mu<- out[[2]]

  ##########################
  ##########################
  ##########################
  ##### compute the scalar products
  index <- seq(1,twoN)
  psi_c<- matrix(1,N,twoN )
  for( i in 1:twoN){
    p_cur <-as.function(Psi[[i]])
    pp <- und_X^(-d)*p_cur(X/und_X)
    psi_c[,i]<- psi_c[,i]* matrix(pp,N,1)
  }
  inner_r_estimate <- as.matrix(exp(1i*t*Y)/f_X)
  inner_r_estimate <- 1/N*t(psi_c)%*%inner_r_estimate


  x.grid_norm <- g1/limit


  ######
  c = und_X * abs(t)
  b=1/(1.03*limit)
  bound=1
  xseq = seq(-bound,bound, length.out=400)
  resol=2^15
  psix <- PSI_mu_fourier(xseq,c,b,Psi,resol)
  psi <- psix[[1]]
  xseq <- psix[[2]]
  xseq2 <- xseq/b

  splin =FALSE
  mu_ev<- MU_fourier(psi,xseq,splin)
  mu <- mu_ev[[1]]
  mu <- unlist(mu)
  mup <-mu%*%matrix(1,1,dim(psi)[2])
  psi_basis_j <- t(psi/(mup)^2*b)
  # dim( psi_basis_j)
  lambda_1_j <- t(as.matrix(sqrt(b/Re(mu)^2),1,twoN ))

  ##############################
  N_u=min(40,max(10,ceiling(log(N)/(2*d)/max(1,lambertW(7*pi/(c/b*abs(t)), tolerance = 1e-10, maxit = 50)))))
  Q = 10^(-20)

  out_stock <- matrix(0,1,N_u)
  Pen  <-  matrix(0,1,N_u)
  U_part_bias <- matrix(0,1,N_u)

  for (N_ad in 2:(min(dim(inner_r_estimate)[1],N_u)-1)){

    ### computation of B1
    res0 <- sum(Re(abs(lambda_1_j[1,N_ad:min( dim(lambda_1_j)[2],N_ad)]))*Mod(inner_r_estimate[N_ad:min( dim(inner_r_estimate)[1],N_ad),1])^2)
    Pen0 = Q*20*(1+2*log(N))*n_f_X/N*(und_X*abs(t))^d*(2*max(1,N_ad)/max(limit*und_X*abs(t),1) )^d*exp(2*d*(N_ad)*log(max(7*pi*(N_ad+1)/c,1)))
    res <- max( res0  - Pen0,0 )
    for( n_kprime in N_ad:min(dim(inner_r_estimate)[1],N_u)){
      res0 <- sum(Re(abs(lambda_1_j[1,N_ad:min( dim(lambda_1_j)[2],max(N_ad,n_kprime))]))*Mod(inner_r_estimate[N_ad:min( dim(inner_r_estimate)[1],max(N_ad,n_kprime)),1])^2)
      Pen0 = Q*20*(1+2*log(N))*n_f_X/N*(und_X*abs(t))^d*(2*max(1,n_kprime)/max(limit*und_X*abs(t),1) )^d*exp(2*d*(n_kprime)*log(max(7*pi*(n_kprime+1)/c,1)))
      res0 <- max( res0  - Pen0,0 )
      if(res0 > res){
        N_adp = n_kprime
        res = res0
      }

    }

    Pen[1,N_ad] = Q*20*(1+2*log(N))*n_f_X/N*(und_X*abs(t))^d*(2*max(1,N_ad)/max(limit*und_X*abs(t),1) )^d*exp(2*d*(N_ad)*log(max(7*pi*(N_ad+1)/c,1)))
    out_stock[1,N_ad]=  res +  Pen[1,N_ad]
    # out_stock[1,N_ad]=  res0
  }
  # out_stock <- out_stock + abs(min( U_part_bias ))
  out_stock

  N_adp = max(1,which.min(abs(out_stock[2:(min(dim(inner_r_estimate)[1],N_u)-1)])) ) +1
  norm_stock <- matrix(0,1,N_u)
  mat <- matrix(0,1,length(g1))
  # beta_model_cur <- beta_model
  resol=2^13
  # output_N <-matrix(F_N, M,1)
  # T_1_opt[t_ind,] =    output_N


  # #       ###### select Psis according to the regularisation rule
  sel <- NULL
  for(i in 1:twoN){
    if(i <= N_adp )
      sel = c(sel,i)
  }

  F_N <- psi_basis_j[,sel]%*%inner_r_estimate[sel,]
  xseq2 <- xseq/b
  Fourier_funRe <- approxfun(xseq2, y = Re(F_N),       method = "linear")
  Fourier_funIm <- approxfun(xseq2, y = Im(F_N),       method = "linear")
  Fourier_fun0 <- function(x){  return(Fourier_funRe(x) + 1i*Fourier_funIm(x))}
  output <-matrix( Fourier_fun0(g1), M,1)
  res0 <- sum(Re(abs(lambda_1_j[,sel]))%*%abs(inner_r_estimate[sel,])^2)
  # norm_stock[1,N_opt]=  sqrt(sum(t(Mod(mat - t(output))^2)))
  # T_1_mat[t_ind,] =   output

  out <- vector("list")
  out[[1]] <-  output
  out[[2]] <-  res0
  return(out)
}



