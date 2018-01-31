library(gpuR)

Y.matgpu <- gpuMatrix(predicted)
  #scale (SD = 1) and center (mean = 0) Y
  Ycgpu <- scale(Y.matgpu)
  Ycgpu <- Ycgpu[,!colSums(!is.finite(Ycgpu))]
  Ycgpu <- gpuMatrix(Ycgpu)

X.matgpu <- gpuMatrix(predictor)
  #scale (SD = 1) and center (mean = 0) X
  Xcrgpu <- scale(X.matgpu)
  Xcrgpu <- Xcrgpu[,!colSums(!is.finite(Xcrgpu))]
  Xcrgpu <- gpuMatrix(Xcrgpu)

  #variable length of X and Y
  p <- ncol(Xcrgpu)
  q <- ncol(Ycgpu)

 CRTs          <- c()
  sum_abs_Betas <- c()
  Nr_iterations = 0         #iteration counter

  WeContinnue = TRUE        #T/F value for stop criterium based on CRT values
                            #in last two iteration
  CRT = 1                   #convergance measure between alpha and beta

#****************
#Start algorithm
#****************

stime_gpu       <-    system.time({

BETAgpu <- gpuMatrix(rep(1, q), nrow = q, ncol = 1)
ALPHAgpu <- gpuMatrix( rep(1, p), nrow = p, ncol = 1)

  while(CRT > tolerance && WeContinnue && Nr_iterations < max_iterations) {

    ETAgpu = Ycgpu %*% BETAgpu
    ETAgpu = scale(ETAgpu)

    XIgpu = Xcrgpu %*% ALPHAgpu
    XIgpu = scale(XIgpu)

    ALPHgpu_0 <- base::solve(t(Xcrgpu) %*% Xcrgpu) %*% (t(Xcrgpu) %*% ETAgpu)

    XIgpu = Xcrgpu %*% ALPHgpu_0
    XIgpu = scale(XIgpu)

    BETAgpu_0 = base::solve(t(XIgpu)%*%XIgpu) %*% (t(XIgpu) %*% Ycgpu)
    BETAgpu_0 = t(BETAgpu_0)

    ETAgpu = Ycgpu %*% BETAgpu_0
    ETAgpu = scale(ETAgpu)

    CRT = sum((ALPHAgpu - ALPHgpu_0)^2, (BETAgpu - BETAgpu_0)^2);

    ALPHAgpu=            ALPHgpu_0
    BETAgpu =            BETAgpu_0

    Nr_iterations                   =   Nr_iterations + 1
    CRTs[[Nr_iterations]]           =   CRT
    sum_abs_Betas[[Nr_iterations]]  =   sum(abs(BETAgpu))

    if (Nr_iterations>1){

      stop_condition <- abs(CRTs[[Nr_iterations]] - CRTs[[Nr_iterations-1]])
      stop_criterium <- 1 * 10^-6

      if (stop_condition < stop_criterium){
        WeContinnue <- FALSE
      }

    }#END Check if last two iterations CR converges*************************#
}# End of main loop
})[3]#end of measure time
stime_gpu
