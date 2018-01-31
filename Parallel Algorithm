library(doParallel)
registerDoParallel(cores=2) #Run this for parallelization on two cores
registerDoParallel(cores=4) #Run this for parallelization on four cores

#************************************
Y.mat <- as.matrix(predicted)
  Yc <- scale(Y.mat)
  Yc <- Yc[,!colSums(!is.finite(Yc))]

X.mat <- as.matrix(predictor)
  Xcr <- scale(X.mat)
  Xcr <- Xcr[,!colSums(!is.finite(Xcr))]

  #variable length of X and Y
  p <- ncol(Xcr)
  q <- ncol(Yc)

  CRTs          <- c()
  sum_abs_Betas <- c()
  Nr_iterations = 0         #iteration counter

  WeContinnue = TRUE        #T/F value for stop criterium based on CRT values
                            #in last two iteration
  CRT = 1                   #convergance measure between alpha and beta

#****************
#Start algorithm
#****************

ptime2 <- system.time({
  BETA = matrix( rep(1, q), nrow = q, ncol = 1, byrow = TRUE)
  ALPHA = matrix( rep(1, p), nrow = p, ncol = 1, byrow = TRUE)
r <- foreach(i=1) %dopar% {

  while(CRT > tolerance && WeContinnue && Nr_iterations < max_iterations) {

    ETA = Yc %*% BETA
    ETA = scale(ETA)

    XI = Xcr %*% ALPHA
    XI = scale(XI)

    ALPH_0 <- solve(t(Xcr) %*% Xcr) %*% t(Xcr) %*% ETA

    XI = Xcr %*% ALPH_0
    XI = scale(XI)

    BETA_0 = solve(t(XI)%*%XI) %*% t(XI) %*% Yc
    BETA_0 = t(as.matrix(BETA_0))

    ETA = Yc %*% BETA_0
    ETA = scale(ETA)

    CRT = sum((ALPHA - ALPH_0)^2, (BETA - BETA_0)^2);

    ALPHA=            ALPH_0
    BETA =            BETA_0

    Nr_iterations                   =   Nr_iterations + 1
    CRTs[[Nr_iterations]]           =   CRT
    sum_abs_Betas[[Nr_iterations]]  =   sum(abs(BETA))

    if (Nr_iterations>1){

      stop_condition <- abs(CRTs[[Nr_iterations]] - CRTs[[Nr_iterations-1]])
      stop_criterium <- 1 * 10^-6

      if (stop_condition < stop_criterium){
        WeContinnue <- FALSE
      }

    }#END Check if last two iterations CR converges*************************#
}# End of main loop
}
})[3]
ptime
