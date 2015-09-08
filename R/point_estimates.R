###################NOTE#############################
#For simplicity we use three different functions for
#the three different cases
#1. Binary simple estimate
#2. Binary treatment for ATT
#3. Multiple treatment arms
###################NOTE##############################
#####################################################

#CASE 1: Binary treatment. Simple case, i.e. not ATT
get.est.simple<- function(ini1,ini2, X, Y, Ti, theta, max.iter = 100,
                          tol = 1e-3,  bt.a = 0.3, bt.b = 0.5,
                          verbose = TRUE){

  #initialize some objects
  N<- length(Y)
  umat<- cbind(1,X)
  ubar <- colMeans(umat)

  #Perform some simple checks
  if(length(ini1) != length( umat[1,] ) ){
    stop("Incorrect length of initial vector")
  }
  if(length(ini2) != length( umat[1,] ) ){
    stop("Incorrect length of initial vector")
  }

  if(verbose){
    cat("Running BFGS algorithm for estimating Weights p: \n")
  }

  #Implement the BFGS algorithm to get the optimal lambda_p
  #and correspoding weights \hat{p}_k
  p.hat<- cpp_quasi_newt(ini1,umat, ubar, Ti, theta, bt.a,bt.b,
                 max.iter , tol)

  if(verbose){
    cat("\nRunning BFGS algorithm for estimating Weights q: \n")
  }

  #Implement the BFGS algorithm to get the optimal lambda_q
  #and correspoding weights \hat{q}_k
  q.hat<- cpp_quasi_newt(ini2,umat, ubar, 1 - Ti, theta, bt.a,bt.b,
                         max.iter, tol)

  #Obtain estimates for E[Y(1)], E[Y(0)] and tau = E[Y(1)] - E[Y(0)]
  Y_one <-  sum((p.hat$weights*Y)[Ti==1])
  Y_zero<- sum((q.hat$weights*Y)[Ti==0])
  tau.hat<-  Y_one - Y_zero

  #Obtain weights and lambda's
  w.p<- p.hat$weights
  w.q<- q.hat$weights

  lam.p<- p.hat$res
  lam.q<- q.hat$res

  #Send warning if algorithm did not converge
  conv = TRUE
  if(!p.hat$Conv | !q.hat$Conv){
    warning("BFGS Algorithm did not converge for atleast one objective function.")
    conv<- FALSE
  }

  #Return list of objects
  res.l<- list("lam.p" = lam.p, "lam.q" = lam.q, "weights.p" = w.p,
               "weights.q" = w.q, "Y1" = Y_one, "Y0"= Y_zero,"tau" = tau.hat,
               "conv" = conv)
  return(res.l)
}

###############################################################################
#Main function to obtain the point estimate for ATT
get.est.ATT<- function(ini2, X, Y, Ti, theta, max.iter = 100,
                       tol = 1e-3,  bt.a = 0.3, bt.b = 0.5,
                       verbose = TRUE){

  #Initialize
  N<- length(Y)
  umat<- cbind(1,X)

  #Main difference is that now u_bar chnages
  umat2<- umat[Ti==1,]
  ubar <- colMeans(umat2)

  if(length(ini2) != length( umat[1,] ) ){
    stop("Incorrect length of initial vector")
  }

  if(verbose){
    cat("\nRunning BFGS algorithm Raphson for estimating Weights q: \n")
  }

  #Obtain lambda'_q and \hat{q}_K for the case of ATT
  q.hat<- cpp_quasi_newt(ini2,umat, ubar, 1 - Ti, theta, bt.a,bt.b,
                         max.iter, tol)

  #Note that treatment effect on treated is simple to estimate
  Y_one <- mean(Y[Ti==1])
  #The calibration estimator for E[Y(0)|T = 1]
  Y_zero<- sum((q.hat$weights*Y)[Ti==0])
  tau.hat<-  Y_one - Y_zero

  #Weights and lambda vectors
  w.q<- q.hat$weights
  lam.q<- q.hat$res

  #Warning message if the algorithm did not converge
  conv = TRUE
  if(!q.hat$Conv){
    warning("\nBFGS algorithm did not converge for the objective function.")
    conv<- FALSE
  }

  #Return list
  res.l<- list("lam.q" = lam.q, "weights.q" = w.q,
               "Y1" = Y_one, "Y0"= Y_zero,"tau" = tau.hat,
               "conv" = conv)
  return(res.l)
}

###############################################################################
#Main function to obtain the point estimate for multiple treatment effect
#Much of the syntax is the same as for the binary case
get.est.MT<- function(ini.mat, X, Y, Ti, theta, max.iter = 100,
                      tol = 1e-3,  bt.a = 0.3, bt.b = 0.5,
                      verbose = TRUE){

  #initialize some objects
  N<- length(Y)
  J<-length(unique(Ti))
  umat<- cbind(1,X)
  ubar <- colMeans(umat)

  if(ncol(ini.mat) != length( umat[1,] ) ){
    stop("Incorrect length of initial vector")
  }
  if(verbose){
    cat("\nRunning BFGS for estimating Weights: \n")
  }

  lam.mat<- matrix(0, ncol = ncol(ini.mat), nrow = J  )
  weights.mat<- matrix(0, ncol = N, nrow = J  )
  Yj.mat<- numeric( J )

  #Loop through the treatment arms
  for(j in 0:(J-1 ) ){
    #Inidicator of treatment arm being treamtent j
    temp.Ti<- 1*(Ti==j)

    #Implement BFGS algorithm
    pj.hat<- cpp_quasi_newt(ini.mat[j+1,],umat, ubar, temp.Ti, theta, bt.a,bt.b,
                           max.iter, tol)

    #Find estimate for  E[Y(j)]
    Yj.hat <-  sum((pj.hat$weights*Y)[temp.Ti==1])

    Yj.mat[j+1]<- Yj.hat
    lam.mat[j+1,]<- pj.hat$res
    weights.mat[j+1,] <- pj.hat$weights

    #Warning for non-convergence
    conv = TRUE
    if(!pj.hat$Conv){
      warning("BFGS algorithm did not converge for atleast one objective function.")
      conv<- FALSE
    }
  }
  #Rerutn result.
  res.l<- list("lam.mat" = lam.mat, "weights.mat" = weights.mat,
               "Yj.hat" = Yj.mat, "conv" = conv)
  return(res.l)
}

###############################################################################
