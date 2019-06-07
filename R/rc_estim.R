
#' Adaptive estimation of the joint density of random coefficients model in a linear random coefficient model
#'
#' This function implements the adaptive estimation of the joint density of random coefficients model as in  Gaillac and Gautier (2019).
#' It takes as inputs data (Y,X) then estimates the density and return its evaluation on a grid b_grid times a_grid.
#' By setting nbCores greater than 1 computations are done in parallel.
#'
#' @param  X  Vector of size $N$, $N$ being the number of observation and the number of regressors limited to 1 in this version of the package.
#' @param  Y Outcome vector of size $N	$.
#' @param  b_grid vector grid on which the estimator of the density of the random slope will we evaluated.
#' @param  a_grid Vector grid on which the estimator of the density of the random intercept will we evaluated.
#' @param  nbCores number of cores for the parallel implementation. Default is 1, no parallel computation.
#' @param  M_T number of discretisation points for the estimated partial Fourier transform. Default is 60.
#' @param  N_u aximal number of singular functions used in the estimation. Default is the maximum of 10 and $N_{max}$.
#' @param  epsilon parameter for the interpolation. Default is (log(N)/log(log(N)))^(-4) as is (T5.1) in Gaillac and Gautier (2019).
#' @param  n_0  Parameter for the sample splitting. If n_0 = 0 then no sample splitting is done and we use the same sample of size $N$ to perform the estimation of the truncated density. If n_0 >0 , then this is the size of the sample used to perform the estimation of the truncated density. Default is $n_0 =0$.\\
#' @param  trunc Dummy for the truncation of the density of the regressors to an hypercube [x_0,x_0]^p. If trunc=1, then truncation is performed and [x_0,x_0]^p is defined using the argmin of the ratio of the estimated constant c_X  over the percentage of observation in  [x_0,x_0]^p. Default is 0, no truncation.
#' @param  center Dummy to trigger the use of X -x_bar instead of $X$ as regressor. If center=1, then use X - x_bar where x_bar is the vector of the medians coordinates by  coordinates for $X$. Default is center=0, where regressors are left unchanged.
#' @return a list containing, in order:
#'
#'  -  outcome : the matrix of size length(b_grid)x length(a_grid) which contains the evaluation of the estimator of the density on the grid b_grid times a_grid
#'
#'  -  b_grid : vector grid on which the estimator of the density is evaluated for the random slope
#'
#'  -  a_grid : vector grid on which the estimator of the density is evaluated for the random intercept
#'
#'  -  x_bar : vector used to center the regressors, if center=1
#'
#'  -  n_f_X : estimated supnorm of the inverse of the density of the regressors
#'
#'  -  und_X : parameter x_0 used for the truncation (truncate =1) of the density of the regressors.
#'
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
#'
##### Estimation
#'\dontshow{
#'set.seed(2019)
#'M=2
#'b_grid  <- seq(-limit ,limit ,length.out=M)
#'N <-5
#'beta <- runif(N,  -limit2,-limit2)
#'X <- as.matrix(runif(N,  -up,-up))
#'X_t <- cbind(matrix(1, N,1),X)
#'Y <-rowSums(beta*X_t)
#'out <- rc_estim(X,Y,b_grid,b_grid,nbCores = 1, M_T = 2,N_u=2)
#'
#'M=100
#'b_grid  <- seq(-limit ,limit ,length.out=M)
#'}
#'\donttest{
#'N <- 1000
#'xi1 <- rtmvnorm(N, mean = Mean_mu1, sigma=Sigma, lower=c( -limit2,-limit2), upper=c(limit2,limit2))
#'xi2 <- rtmvnorm(N, mean = Mean_mu2, sigma=Sigma, lower=c( -limit2,-limit2), upper=c(limit2,limit2))
#'theta = runif(N, -1 , 1)
#'beta <- 1*(theta >=0) * xi1  + 1*(theta <0) * xi2
#'X <- rtmvnorm(N, mean = c(0), sigma=2.5, lower=c( down), upper=c(up))
#'X_t <- cbind(matrix(1, N,1),X)
#'Y <-rowSums(beta*X_t)
#'out <- rc_estim( X,Y,b_grid,b_grid,nbCores = 1, M_T = 60)
#'}
#'
#'

rc_estim <- function( X,Y,b_grid,a_grid  = b_grid,nbCores = 1, M_T = 60 , N_u=10, epsilon = (log(N)/log(log(N)))^(-2), n_0 = 0, trunc=0, center=0){

  X <- as.matrix(X)
  ## sample splitting option
  if( n_0 != 0){
    select0 = sample(1:length(Y), n_0 )
    X_sample <- X[select0]
    X <- as.matrix(X[-select0,])
    Y <- Y[-select0]
  }else{
    X_sample <- X
    n_0 = length(Y)
  }

  ### estimation of f_X
  if(center !=0){
    x_bar = colMedians(X)
    X <- X - x_bar
    X_sample <- X - x_bar
  }else{
    x_bar =  rep(0,d)
  }


  ## truncation based on the density of X
  if( trunc != 0){
    xq_seq <- seq(0.50,1.0, by = 0.1)
    qq_x = NULL
    for(kk in 1:d){
      qq_x <- rbind(qq_x,c(quantile( abs(X),xq_seq), 1.05*max(abs(abs(X)))))
    }
    xq_seq <- c(xq_seq,1.05)

    out_cx <- NULL
    for (kk in 1:length(qq_x)){
      und_X_1= qq_x[kk]
      est <- kde( X_sample[1:min(1000,n_0)], gridsize=c(n_0), xmin=rep(-1,d),xmax=rep(1,d) )
      po <- ceimOpt(OptimFunction=function(x){-1/predict(est,x=x)*(abs(x)<und_X_1)}, maxIter=30, epsilon=0.3,nParams=d)
      n_f_X = -po$BestMember[2]
      out_cx <- rbind(out_cx, c(  n_f_X  ,mean((abs( X_sample)<und_X_1))))
    }

    indic_truncate  = which.min(out_cx[,1]/out_cx[,2])
    n_f_X  = out_cx[indic_truncate,1]
    und_X= qq_x[indic_truncate]
    select0 = abs(X)>=und_X
    select1 = abs(X_sample)>=und_X
    X_sample <-  X_sample[-select1]
    n_0= length(X_sample)
    X <- as.matrix(X[-select0,])
    Y <- Y[-select0]
    est <- kde( X_sample[1:min(1000,n_0)], gridsize=c(n_0), xmin=rep(-1,d),xmax=rep(1,d) )
    f_X <- predict(est,x=X)
    f_X[f_X==0] <- min(f_X[f_X>0])
  }else{
    und_X= 1.05*max(abs(abs(X)))
    if(length(Y) <= 10){
      n_f_X = 1/(2*und_X)
      f_X <-rep(1/(2*und_X),length(Y))
      f_X[f_X==0] <- min(f_X[f_X>0])
    }else{
      est <- kde( X_sample[1:min(1000,n_0)], gridsize=c(n_0), xmin=rep(-1,d),xmax=rep(1,d) )
      po <- ceimOpt(OptimFunction=function(x){-1/predict(est,x=x)*(abs(x)<und_X)}, maxIter=30, epsilon=0.3,nParams=d)
      n_f_X = -po$BestMember[2]
      f_X <- predict(est,x=X)
      f_X[f_X==0] <- min(f_X[f_X>0])
    }
  }


  #### High level parameters
  a = max(abs(a_grid))
  und_X = 1.05*max(abs(X))
  d = dim(X)[2]
  N = length(Y)
  M= length(b_grid)
  limit = max(abs(b_grid))
  L =15


  #################################################################################################################################
  #   Interpolation initialisation
  #################################################################################################################################
  #### Number of Psis
  L1 = L+1
  N2 = max(L,3)
  twoN = 2*N2


  #### Bandwidth 1
  c1 = 1
  K1 = max(twoN+2,30)
  K = K1

  #### get the beta's (for the computation of Psi)

  out <- get_psi_mu(c1,N2,twoN,K1, L1)
  Psi1 <- out[[1]]
  mu1<- out[[2]]

  if( M_T <= 2){
    outcome<- array(0, c(M,M))
  }else{
  #### Tuning parameters
  T_lim = epsilon +  log(N)
  T_seq  <- sort(c(seq(-T_lim,-epsilon,length.out=M_T/2),seq(epsilon,T_lim,length.out=M_T/2)))
  M_T  =length(T_seq)
  delta_T <- 2*(T_lim -epsilon)/M_T

  M_eps =max(floor(2*epsilon/ delta_T), 2)

  ######
  T_1_outcome<- array(0, c(M_T,rep(M,d)))
  T_1_opt<- array(0, c(M_T,rep(M,d)))
  #### Fourier transform of the Marginal f_a of alpha
  Tf_a <- matrix(0,1,M_T)
  #### Fourier transform of the Marginal f_a of alpha
  Tf_a2 <- matrix(0,1,M_T)
  ##### Loop on the values of the parameter t
  stock_T <- matrix(0,1,M_T)
  ##

  #### if one core, then direct work
  if(nbCores==1){
    ######## Bootstrap test  b_rep=1
    for (t_ind in 1:M_T){
      out <-boot_stat(t_ind ,T_seq,und_X ,twoN,N2,L1,d,f_X,Y,X,N,limit,b_grid,M,  n_f_X )
      T_1_outcome[t_ind,] <- out[[1]]
      stock_T[1,t_ind] <- out[[2]]
    }
  }else{
    # nbCores=4
    sfInit(parallel=TRUE, cpus=nbCores, type="SOCK")
    # sfExportAll( except=c() )
    sfExportAll()
    #sfExport("c_cube")
    sfLibrary("orthopolynom",character.only=TRUE)
    sfLibrary(fourierin)
    sfLibrary("sfsmisc",character.only=TRUE)

    out <-sfLapply( 1:M_T, boot_stat,T_seq,und_X ,twoN,N2,L1,d,f_X,Y,X,N,limit,b_grid,M,  n_f_X)

    sfStop()
    # T_reps2 <- unlist(  T_reps2)
    for (t_ind in 1:M_T){
      out0 <-  out[[t_ind]]
      T_1_outcome[t_ind,] <- out0[[1]]
      stock_T[1,t_ind] <- out0[[2]]
    }
  }

  #   Adaptive choice of T
  T_scale <- sort(T_seq)
  T_min=1
  T_scale2 = T_scale[  T_scale >T_min]
  pas=1
  T_scale2 = T_scale2[seq(1,length(T_scale2),by=pas)]
  SUM2 = NULL
  SUM3 <-  matrix(0,1,length(T_scale2))
  SUM3pp <-  matrix(0,1,length(T_scale2))
  SUM3res <-  matrix(0,1,length(T_scale2))
  T_table <- matrix(0,1,length(T_scale2))
  Q = 10^(-5)
  ind=1
    for(ind in 1:length(T_scale2) ){
      t= T_scale2[ind]

      res = 0
      t_adp =0
      ## computation of B2
      for( t_kprime in ind:(length(T_scale2)-1)){
        select <- stock_T[1,(T_scale2[t_kprime] >= abs(T_seq))& ( abs(T_seq)>= t)]
        Pen1 <- function(t){
          N_u = min(40,max(10,ceiling(log(N)/(2*d)/max(1,lambertW(7*pi/(limit*und_X*abs(t)), tolerance = 1e-10, maxit = 50)))))

          return(Q*20*(1+2*log(N))*n_f_X/N*(und_X*abs(t))^d*(2*max(1,N_u)/max(limit*und_X*abs(t),1) )^d*exp(2*d*N_u*log(max(7*pi*(N_u+1)/(limit*und_X*abs(t)),1))))
        }
        res0 <- sum(pmax(sum(select) -  Pen1(T_scale2[1:t_kprime]),0     )*(2*(T_scale2[t_kprime]-t)/(max(1,length(select)-1))))
        if(res0 > res){
          t_adp = t_kprime
          res = res0
        }

      }

      PP <- sum(Pen1(T_scale2[1:ind]))*(2*t/(max(1,length(1:ind)-1)))
      SUM3[1,ind]=  res +   PP
      SUM3pp[1,ind]=   PP
      SUM3res[1,ind]=  res
    }

    num <-min(1+ which.min( SUM3[2:length(T_scale2)]),length(T_table[1,]))
  T_adpt <- T_scale2[num ]

  M_T_adpt <-T_table[1,num]
  T_1_adpt<- array(0, c(M_T_adpt,rep(M,d)))
  t_vect <- NULL
  for(i in 1:length(T_seq)){
    if(abs(T_seq[i])<= T_adpt){
      t_vect <- c(t_vect,i)
    }
  }
    T_seq_a <- T_seq[t_vect]
  T_1_adpt<- T_1_outcome[t_vect,1:M]
  M_T_adpt <- dim( T_1_adpt)[1]



  #############################################################################
  ################################################################################
  ##### compute the scalar products
  sizeT_a <- length(T_seq_a)
  lim <-3/4
  N_I=20

  T_1_adp_eps = matrix(0,M_eps,M)
  min_T <- min(abs(T_seq_a ))
  max_T <- max(abs(T_seq_a ))
  dif_T <-( max_T - min_T)/sizeT_a
  m_ind=1
  for(m_ind in 1:M){


    m=b_grid[m_ind]

    index1 <- seq(1,twoN)
    psi_c1<- matrix(1,sizeT_a,twoN )

    for( i in 1:twoN){
      p_cur1 <-as.function(Psi1[[i]])
      pp1 <- p_cur1( T_seq_a )
      psi_c1[,i]<- psi_c1[,i]* matrix(pp1,sizeT_a,1)
    }

    inner_r_estimate1 <- t(psi_c1)%*%T_1_adpt[,m_ind ]*dif_T
    grid_eps <- seq(-epsilon, epsilon, length.out=M_eps)

    #### form  the psi basis
    psi_basis_j1<- matrix(1,length(grid_eps),twoN )
    lambda_1_j1<- matrix(1,1,twoN )
    for( i in 1:twoN){
      p_cur1 <-as.function(Psi1[[i]])
      pp1 <- p_cur1(grid_eps)
      lambda_1_j1[,i]<-min(lim , abs(mu1[[i]])^2*(c1/2*pi))
      psi_basis_j1[,i]<- psi_basis_j1[,i]*lambda_1_j1[,i]/(1-lambda_1_j1[,i])* matrix(pp1,length(grid_eps),1)
    }
    sel <- NULL
    for(i in 1:twoN){
      if( i<= N_I)
        sel = c(sel,i)
    }
    F_N1 <- psi_basis_j1[,sel]%*%inner_r_estimate1[sel,] /(2*pi)
    output_N1 <-matrix(F_N1,M_eps ,1)
    T_1_adp_eps[,m_ind] =    output_N1

  }
  T_seq_a2<- c(T_seq_a[1:(M_T_adpt /2)],grid_eps, T_seq_a[(M_T_adpt /2+1):dim(T_1_adpt)[1]] )
  ##

  TOUT_a <- rbind(T_1_adpt[1:(M_T_adpt/2),], T_1_adp_eps)
  TOUT_a <-rbind(TOUT_a, T_1_adpt[(M_T_adpt /2+1):dim(T_1_adpt)[1],])
  M_T_adpt= dim(TOUT_a)[1]

  ##########################################################################
  ##########################################################################

  ### computation of f(t)
  #### integrate the partial fourier transform estimation over t (fourier variable related to the intercept)
  M_a =M
  ######
  T_1_cur1<- array(0, c(dim(TOUT_a)[1]/2,M))
  outcome<- array(0, c(M_a,M))
  maxT <- max(T_seq_a2)
  for( a_ind in 1:M_a){
    a_grid0 <- a_grid[a_ind]
    for( i in 1:(dim(TOUT_a)[1]/2)){
      #     T_1_cur[i,  ]<-  exp(-T_seq[i]*1i*a_grid) * T_1_adpt[i,  ]
      T_1_cur1[i,  ]<-  exp(-T_seq_a2[i+dim(TOUT_a)[1]/2]*1i*a_grid0) * TOUT_a[i+dim(TOUT_a)[1]/2,  ]
    }
    for( i in 1:M){
      outcome[a_ind, i ]<- sum(T_1_cur1[,i ])* 2* maxT/dim(TOUT_a)[1]

    }
  }

  res = 2^13
  outcome<- array(0, c(M,M))
  dim(TOUT_a)
  length(b_grid)
  for( i in 1:M){
    pst <- function(t){
      outRe <- approxfun(T_seq_a2, Re(TOUT_a[,i  ]))
      outIm <- approxfun(T_seq_a2, Im(TOUT_a[,i  ]))
      return(outRe(t) + 1i*outIm(t))
    }
    pst0 <- fourierin_1d(  pst,0, max(T_seq_a2),lower_eval = min(b_grid), upper_eval = max(b_grid),freq_adj=-1,const_adj=1,resolution =res)
    out1 <- approxfun( pst0$w, Re(pst0$values))
    outcome[,i]<- out1(b_grid)
  }

  }
  outcome <-Re(outcome)/2
  outcome[outcome<0]<- 0
  outcome[is.na(outcome)] <- 0

  out_export <-  vector("list")
  out_export[[1]] <- outcome
  out_export[[2]] <- b_grid
  out_export[[3]] <- a_grid
  out_export[[4]] <- x_bar
  out_export[[5]] <- n_f_X
  out_export[[6]] <- und_X

  return(out_export)

}
