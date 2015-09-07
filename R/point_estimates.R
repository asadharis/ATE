#We note that we combine the newton raphson for both the simple estimate
#and the treatment effect of the treated since there is only a small modification
#needed.


#########################################################################
#Main function to obtain the simple point estimate
get.est.simple<- function(ini1,ini2, X, Y, Ti, theta, max.iter = 100,
                          tol = 1e-3,  bt.a = 0.3, bt.b = 0.5,
                          verbose = TRUE){

  N<- length(Y)

  umat<- cbind(1,X)
  ubar <- colMeans(umat)

  if(length(ini1) != length( umat[1,] ) ){
    stop("Incorrect length of initial vector")
  }
  if(length(ini2) != length( umat[1,] ) ){
    stop("Incorrect length of initial vector")
  }

  if(verbose){
    cat("Running BFGS algorithm for estimating Weights p: \n")
  }


  p.hat<- cpp_quasi_newt(ini1,umat, ubar, Ti, theta, bt.a,bt.b,
                 max.iter , tol)

  if(verbose){
    cat("\nRunning BFGS algorithm for estimating Weights q: \n")
  }

  q.hat<- cpp_quasi_newt(ini2,umat, ubar, 1 - Ti, theta, bt.a,bt.b,
                         max.iter, tol)

  Y_one <-  sum((p.hat$weights*Y)[Ti==1])
  Y_zero<- sum((q.hat$weights*Y)[Ti==0])
  tau.hat<-  Y_one - Y_zero

  w.p<- p.hat$weights
  w.q<- q.hat$weights

  lam.p<- p.hat$res
  lam.q<- q.hat$res

  conv = TRUE
  if(!p.hat$Conv | !q.hat$Conv){
    warning("BFGS Algorithm did not converge for atleast one objective function.")
    conv<- FALSE
  }

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

  N<- length(Y)
  umat<- cbind(1,X)

  umat2<- umat[Ti==1,]
  ubar <- colMeans(umat2)

  if(length(ini2) != length( umat[1,] ) ){
    stop("Incorrect length of initial vector")
  }


  if(verbose){
    cat("\nRunning BFGS algorithm Raphson for estimating Weights q: \n")
  }
  q.hat<- cpp_quasi_newt(ini2,umat, ubar, 1 - Ti, theta, bt.a,bt.b,
                         max.iter, tol)

  Y_one <- mean(Y[Ti==1])
  Y_zero<- sum((q.hat$weights*Y)[Ti==0])
  tau.hat<-  Y_one - Y_zero

  w.q<- q.hat$weights
  lam.q<- q.hat$res

  conv = TRUE
  if(!q.hat$Conv){
    warning("\nBFGS algorithm did not converge for the objective function.")
    conv<- FALSE
  }

  res.l<- list("lam.q" = lam.q, "weights.q" = w.q,
               "Y1" = Y_one, "Y0"= Y_zero,"tau" = tau.hat,
               "conv" = conv)
  return(res.l)
}

###############################################################################
#Main function to obtain the point estimate for multiple treatment effect
get.est.MT<- function(ini.mat, X, Y, Ti, theta, max.iter = 100,
                      tol = 1e-3,  bt.a = 0.3, bt.b = 0.5,
                      verbose = TRUE){

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

  for(j in 0:(J-1 ) ){
    temp.Ti<- 1*(Ti==j)

    pj.hat<- cpp_quasi_newt(ini.mat[j+1,],umat, ubar, temp.Ti, theta, bt.a,bt.b,
                           max.iter, tol)

    Yj.hat <-  sum((pj.hat$weights*Y)[temp.Ti==1])

    Yj.mat[j+1]<- Yj.hat
    lam.mat[j+1,]<- pj.hat$res
    weights.mat[j+1,] <- pj.hat$weights

    conv = TRUE
    if(!pj.hat$Conv){
      warning("BFGS algorithm did not converge for atleast one objective function.")
      conv<- FALSE
    }
  }

  res.l<- list("lam.mat" = lam.mat, "weights.mat" = weights.mat,
               "Yj.hat" = Yj.mat, "conv" = conv)
  return(res.l)
}

###############################################################################
