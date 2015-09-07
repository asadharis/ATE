

#########################################################
#########################################################
#We begin with the estimating equations for the simple estimate
get.gk.simple<- function(u,Y,Ti,lam,beta,tau1,tau0,theta,ulam,ubeta){

  #N<- length(Y)
  #ulam<- u%*%lam
  #ubeta<- u%*%beta

  temp1<- Ti*cpp_dcr_rho(ulam,theta)
  gk1<- apply(u,2,"*", temp1)-u

  temp2<- (1-Ti)*cpp_dcr_rho(ubeta,theta)
  gk2<- apply(u,2,"*", temp2)-u

  gk3<- (Ti*Y)*cpp_dcr_rho(ulam,theta) - tau1
  gk4<- ((1-Ti)*Y)*cpp_dcr_rho(ubeta,theta) - tau0

  cbind(gk1,gk2,gk3,gk4)
}


#This function obtains the big covariance matrix for all coefficients
#THe bottom right element of this will be the variance estimate for tau-simple
get.cov.simple<- function(u, Y, Ti, obj,theta){
  N<- length(Y)
  n1<- sum(Ti)
  lam<- (obj$lam.p)
  beta<- (obj$lam.q)
  ulam<- u%*%lam
  ubeta<- u%*%beta

  tau1<- obj$Y1
  tau0<- obj$Y0
  K<- length(lam)

  gk<- get.gk.simple(u,Y,Ti,lam,beta,tau1,tau0,theta,ulam,ubeta)

  meat<- crossprod(gk)/N

  temp1<- Ti*cpp_ddcr_rho(ulam, theta)
  temp2 <- (1-Ti)*cpp_ddcr_rho(ubeta, theta)
  #We define the matrix A
  tempA1<- apply(u,2,"*",temp1)
  A1<- crossprod(tempA1,u)/N
  tempA2<- apply(u,2,"*", temp2)
  A2<- crossprod(tempA2,u)/N
  A<- bdiag(A1,A2)

  C1<- crossprod(u,temp1*Y)
  C2<- crossprod(u,temp2*Y)
  C<- t(bdiag(C1,C2))/N

  Ainv<- bdiag(solve(A1),solve(A2))
  tempMat<- Matrix(0,ncol = 2, nrow = 2*K)
  bread<- cbind( rbind(Ainv,C%*%Ainv), rbind(tempMat, -diag(rep(1,2)))  )

  (bread%*% tcrossprod(meat,bread))/N
}


###############################################################################

get.gk.ATT<- function(u,Y,Ti,beta,tau1,tau0,theta,ubeta){

  n1<- sum(Ti)
  N<- length(Y)
  frac<-n1/N
  #ubeta<- u%*%beta

  temp1<-  (1-Ti)*cpp_dcr_rho(ubeta, theta)
  gk1<- apply(u,2,"*", temp1) - apply(u,2,"*", Ti)/frac
  #gk1<- (1-Ti)*rho1(crossprod(be, uk),...)*uk - 1/frac*Ti.int*uk
  gk2<- Ti*(Y -  tau1)/frac
  #gk2<- N/n1*(Ti.int*(Y.scaler- tau1))
  gk3<- temp1*Y - tau0
  #gk3<- (1-Ti.int)*rho1(crossprod(be, uk),...)*Y.scaler - tau0
  gk4<- Ti-frac

  cbind(gk1,gk4,gk2,gk3)
}


#This function obtains the big covariance matrix for all coefficients for ATT
#THe bottom right element of this will be the variance estimate for tau-ATT
get.cov.ATT<- function(u, Y, Ti, obj,theta){
  N<- length(Y)
  n1<- sum(Ti)
  frac<-n1/N

  beta<- (obj$lam.q)
  ubeta<- u%*%beta

  tau1<- obj$Y1
  tau0<- obj$Y0
  K<- length(beta)

  gk<- get.gk.ATT(u,Y,Ti,beta,tau1,tau0,theta,ubeta)
  meat<- crossprod(gk)/N

  temp1 <- (1-Ti)*cpp_ddcr_rho(ubeta, theta)
  tempA<- apply(u,2,"*",temp1)
  A<- cbind(crossprod(tempA,u) , crossprod(u,Ti)/(frac^2))/N
  A<- rbind(A,c(rep(0,K),-1) )

  C<- Matrix(0, ncol = K, nrow = 2)
  C[2,]<- crossprod(u,temp1*Y)
  C<- cbind(C,c(0,0))/N

  Ainv<- solve(A)
  tempMat<- Matrix(0,ncol = 2, nrow = K+1)
  bread<- cbind( rbind(Ainv,C%*%Ainv), rbind(tempMat, -diag(rep(1,2))) )

  (bread%*%tcrossprod(meat, bread))/N
}


###############################################################################


get.gk.MT<- function(u,Y,Ti,lams, taus,theta, ulams ){
  J<- nrow(lams)
  gk1list<- vector("list", J)
  gk2list<- vector("list", J)

  for(j in 1:J){
    ulam<- ulams[,j]
    temp1<- 1*(Ti== j-1)*cpp_dcr_rho(ulam,theta)
    gk1list[[j]] <- apply(u,2,"*", temp1)-u
    gk2list[[j]]<- 1*(Ti== j-1)*Y*cpp_dcr_rho(ulam,theta) - taus[j]
  }
  cbind(do.call(cbind, gk1list), do.call(cbind, gk2list) )

}

#This function obtains the BIG covariance matrix for all coefficients for MT
#This time there is no bottom right element.
get.cov.MT<- function(u, Y, Ti, obj,theta){
  N<- length(Y)
  lams<- obj$lam.mat
  K<- ncol(lams)
  J<- length(unique(Ti))
  taus<- obj$Yj.hat
  ulams<- tcrossprod(u,lams)

  gk<- get.gk.MT(u, Y, Ti, lams, taus, theta, ulams)

  meat<- crossprod(gk)/N

  Alist<- vector("list", J)
  AinvList<- vector("list", J)

  Clist<- vector("list", J)
  for(j in 1:J){
    ulam<- ulams[,j]
    temp1<- 1*(Ti== j-1)*cpp_ddcr_rho(ulam, theta)
    #We define the matrix A
    tempA<- apply(u,2,"*",temp1)
    Alist[[j]] <- crossprod(tempA,u)/N
    AinvList[[j]]<- solve(Alist[[j]])
    Clist[[j]] <- t(crossprod(u,temp1*Y)/N)

  }
  Ainv<- bdiag(AinvList)
  C<- bdiag(Clist)
  tempMat<- Matrix(0,ncol = J, nrow = J*K)

  bread<- cbind( rbind( Ainv ,C%*%Ainv ), rbind(tempMat, -diag(rep(1,J))))
  (bread%*%tcrossprod(meat, bread))/N
}

###############################################################################


estimate_variance<- function(object){
  K<-  object$K
  u<- cbind(1,object$X)
  if(object$gp == "simple"){
    cov.mat<- get.cov.simple(u, object$Y, object$Ti, object, object$theta)
    covEY<- cov.mat[-(1:(2*K) ),-(1:(2*K))]
    fin<- covEY
    #fin[3,3]<- covEY[1,1]+covEY[2,2] - 2*covEY[1,2]

  }else if(object$gp == "ATT"){
    cov.mat<-  get.cov.ATT(u, object$Y, object$Ti,  object, object$theta)

    covEY<- cov.mat[-(1:K),-(1:K)]
    covEY<- covEY[-1,-1]
    fin<- covEY
  }else{
    cov.mat<-  get.cov.MT(u, object$Y, object$Ti, object, object$theta)
    covEY<- cov.mat[-(1:(object$J*K)),-(1:(object$J*K))]
    fin<- covEY
  }
  return(fin)
}
