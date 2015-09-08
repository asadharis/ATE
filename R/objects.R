#This is the main function available to the user
#This creates the object with some routine checks
#This also gives us the point estimates and can be used
#with some summary and plot functions
ATE<- function(Y, Ti, X, theta=0,
               ATT = FALSE, verbose = FALSE,
               max.iter = 100, tol = 1e-10, initial.values = NULL,
               backtrack.alpha = 0.3, backtrack.beta = 0.5){
  #The alpha and beta values for backtracking line search
  bt.a<- backtrack.alpha
  bt.b<- backtrack.beta

  #J is the number treatment arms
  J<- length(unique(Ti))
  #The totoal length of the parameter vector
  #Recall: Using the notation of the Package Vignette we
  #have u_K(X) = cbind(1,X)
  K<- ncol(X)+1

  #In some cases if we have a data.frame we need to
  #convert it into a numeric matrix
  if(class(X)== "data.frame"){
    X<- as.matrix(X)
  }

  #Take care of the case of a single covariate
  if(is.vector(X)){
    X<- matrix(X, ncol = 1, nrow = length(X))
    warning("Data matrix 'X' is a vector, will be treated as n x 1 matrix")
  }

  #Some simple checks before running the main functions
  if(nrow(X)!=length(Y)){
    stop("Dimensions of covariates and response do not match")
  }
  if(J==1){
    stop("There must be atleast two treatment arms")
  }
  if(!all(0:(J-1) == sort(unique(Ti)))){
    stop("The treatment levels must be labelled 0,1,2,...")
  }
  if(!is.numeric(theta)){
    stop("theta must be a real number")
  }
  if(J>2 & ATT){
    stop("Treatment effect on the treated cannot be calculated for multiple treatment groups.")
  }

  #If the inital estimates is not provided we generate
  #simple initial parameters
  if(is.null(initial.values)){
    if(ATT){
      initial.values<- numeric(K)
    }else{
      initial.values<- matrix(0, ncol = K, nrow = J )
    }
  }
  #For ATT there is only one BFGS applied
  if(ATT & !is.vector(initial.values)){
    stop("For ATT we only need one vector of initial values for BFGS")
  }
  #In all other cases this has to be a matrix.
  if(!ATT){
    if( !is.matrix(initial.values) ){
      stop("Initial values must be a matrix")
    }
    if(any(dim(initial.values) !=  c(J,K)) ){
      stop("Matrix of initial values must have dimensions J x K.")
    }
  }

  #Now we specify the category of the problem
  #The simple case with binary treatment
  gp<- "simple"

  #The case of average treatment effect on the treated
  if(ATT) gp<- "ATT"

  #The case of Multiple treatment arms
  if(J>2) gp<- "MT"

  #Obtain the estimates for whatever the case may be
  if(gp == "simple"){
    ini1<- initial.values[1,]
    ini2<- initial.values[2,]
    est<- get.est.simple(ini1,ini2, X, Y, Ti, theta, max.iter ,
                              tol,  bt.a , bt.b ,verbose)

  }else if(gp == "ATT"){
    ini2<- initial.values
    est<- get.est.ATT(ini2, X, Y, Ti, theta, max.iter,
                           tol,  bt.a , bt.b , verbose )
  }else if(gp == "MT"){
    est<- get.est.MT(initial.values, X, Y, Ti, theta, max.iter,
                          tol,  bt.a , bt.b ,
                          verbose)
  }

  #Begin building the "ATE" object which is a list
  res<- est
  res$X<- X
  res$Y<- Y
  res$Ti<- Ti
  res$theta<- theta
  res$gp<- gp
  res$J<- J
  res$K<- K

  if(verbose){
    cat("\nEstimating Variance")
  }

  #Estimate the variance covariance matrix of
  #E[Y(1)] and E[Y(0)]
  res$vcov<- estimate_variance(res)
  res$call<- match.call()

  #Rename some of the elements of the list and organize the output object
  if(gp=="simple"){
    est<- c(res$Y1,res$Y0,res$tau)
    names(est)<- c("E[Y(1)]", "E[Y(0)]", "ATE")
    res$est<- est
    res$Y1<- NULL
    res$Y0<- NULL
    res$tau<- NULL
  }else if(gp == "ATT"){
    est<- c(res$Y1,res$Y0,res$tau)
    names(est)<- c("E[Y(1)|T=1]", "E[Y(0)|T=1]", "ATT")
    res$est<- est
    res$Y1<- NULL
    res$Y0<- NULL
    res$tau<- NULL
  }else{
    est<- res$Yj.hat
    names(est)<- paste("E[Y(",0:(J-1),")]",sep = "")
    res$est<- est
    res$Yj.hat<- NULL
  }

  #Define the class ATE
  class(res)<- "ATE"
  return(res)
}

#S3 print method for class ATE
#The function prints the point estimates and
print.ATE<- function(x, ...){
  object<- x
  if(object$gp == "simple"){
    cat("Call:\n")
    print(object$call)
    cat("\nThe analysis was completed for a simple study design with binary treatment.\n")
    cat("\nPoint Estimates:\n")
    print(object$est)
  }else if(object$gp == "ATT"){
    cat("Call:\n")
    print(object$call)
    cat("\nThe analysis was completed for a binary treatment for estimating treatment effect on the treated.\n")
    cat("\nPoint Estimates:\n")
    print(object$est)
  }else{
    cat("Call:\n")
    print(object$call)
    cat("\nThe analysis was completed for a study design with multiple treatment arms.\n")
    cat("\nPoint Estimates:\n")
    print(object$est)
  }
}


#S3 summary method for ATE
#THis function calculates the SE, Z-statistic and
#confidence intervals and P values
summary.ATE<- function(object, ...){

  #For binary treatment we calculate the variance of
  # EY(1) - EY(0). Using the formula var(a-b) = var(a)+var(b)-2cov(a,b)
  if(object$gp== "simple" || object$gp== "ATT"){
    var.tau<- object$vcov[1,1]+object$vcov[2,2]-
      2*object$vcov[1,2]
    se<- c(sqrt(diag(object$vcov)), sqrt(var.tau))
  }else{#The case of multiple treatments
    #Now we do not have a specific variable like EY(1)-EY(0)
    #So we just obtain the SE for EY(0), EY(1), EY(2),...
    se<- sqrt(diag(object$vcov))
  }
  #Obtain confidence intervals
  Ci.l<- object$est+se*qnorm(0.025)
  Ci.u<- object$est+se*qnorm(0.975)

  #Evaluate the Z-statistic and P-value
  z.stat<- object$est/se
  p.values<- 2*pnorm(-abs(z.stat))
  #Create a coefficient matrix which we can print later
  #to make it look like the out put of summary for lm/glm etc.
  coef<- cbind(Estimate = object$est,
               StdErr = se,
               "95%.Lower" = Ci.l,
               "95%.Upper" = Ci.u,
               Z.value = z.stat,
               p.value = p.values)

  if(object$gp == "simple"){
    res<- list(call = object$call, Estimate = coef, vcov = object$vcov,
               Conv = object$conv, Weights.p= object$weights.p,
               weights.q = object$weights.q)

  }else if(object$gp == "ATT"){
    res<- list(call = object$call, Estimate = coef, vcov = object$cov ,
               Conv = object$conv, weights.q = object$weights.q)

  }else{
    res<- list(call = object$call, Estimate = coef, vcov = object$cov ,
               Conv = object$conv, weights = object$weights.mat)
  }
  class(res)<- "summary.ATE"
  return(res)

}

#A Print method for "summary.ATE"
#uses the coefficient matrix to give us print statements like
#the summary function of lm/glm.
print.summary.ATE <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$Estimate, P.values = TRUE, has.Pvalue=TRUE)
}


#A simple function to evaluate the emiprical CDF
#This function allows us to calculate a weighted emiprical CDF too.
my.ecdf<- function(t,x, weights = NULL){
  if(is.null(weights)){
    sum(1*(x <= t))/length(x)
  }else{
    sum(1*(x<=t)*weights)
  }
}

#The S3 plot function for ATE
#This function plots the empirical CDF for
#all the covariates for the different tretment groups
#And the weighted eCDF after covariate balancing
plot.ATE<- function(x, ...){
  object<- x
  Ti<- object$Ti

  ##################Case 1: Simple binary treatment###################
  if(object$gp == "simple"){
    #Obtain weights
    w.p<- object$weights.p[Ti==1]
    w.q<- object$weights.q[Ti==0]

        #Check if covariates are named
    names<- colnames(object$X)
    if(is.null(names)){
      p<- ncol(object$X)
      names<- paste("X",1:p,sep = "")
    }
    #Obtain the subsample of subjects for each treatment arm
    x1<- as.matrix(object$X[Ti==1,])
    x0<- as.matrix(object$X[Ti==0,])
    for(i in 1:ncol(x1)){
      if(i==2) par(ask = TRUE)

      if(length(unique(object$X[,i])) == 2){
        Treatment<- x1[,i]
        Placebo<- x0[,i]
        plot(c(0.5,1,2,2.5), c(2,mean(Treatment), mean(Placebo),2), pch = 16, cex = 1.5,
             ylim = c(0,1), ylab = "Mean of group", xlab = "",col = c("blue","red"),
             main = "Unweighted", xaxt = "n")
        axis(side = 1, at = c(1,2), labels = c("Treatment", "Control") )
        abline(h = mean(object$X[,i]), lty = 2)
        new_treat<- sum(w.p*Treatment)
        new_control<- sum(w.q*Placebo)
        plot(c(0.5,1,2,2.5), c(2, new_treat, new_control,2), pch = 16, cex = 1.5,
             ylim = c(0,1), ylab = "Mean of group", xlab = "",col = c("blue","red"),
             main = "Weighted", xaxt = "n")
        axis(side = 1, at = c(1,2), labels = c("Treatment", "Control") )
        abline(h = mean(object$X[,i]), lty = 2)

      }else{
        rng<- range(c(x1[,i],x0[,i]))
        my.seq<- seq(rng[1],rng[2],length = 100)
        temp1<- sapply(my.seq, my.ecdf,x = x1[,i])
        temp0<- sapply(my.seq, my.ecdf,x = x0[,i])
        par(mfrow = c(1,2))
        plot(my.seq,temp1,cex = 0.4,pch = 16,type = "l",lty = 1,col = "red",
             xlab = names[i],ylab = "empirical CDF",main = "Unweighted empirical CDF")
        lines(my.seq, temp0, cex = 0.4, pch = 16, lty = 2,col = "blue")
        legend("bottomright", c("Treatment", "Control"), lty = c(1,2), col = c("red", "blue"))

        temp1<- sapply(my.seq, my.ecdf,x = x1[,i],weights = w.p)
        temp0<- sapply(my.seq, my.ecdf,x = x0[,i], weights = w.q)
        plot(my.seq,temp1,cex = 0.4,pch = 16,type = "l",lty = 1,col = "red",
             xlab = names[i],ylab = "emperical CDF",main = "Weighted emperical CDF")
        lines(my.seq, temp0, cex = 0.4, pch = 16, lty = 2,col = "blue")
      }
    }
    par(ask = FALSE)
    #Second, for average treatment effect on the treated
    #We only balance on group in this case
  }else if(object$gp == "ATT"){
    #w.p<- object$weights.p[Ti==1]
    w.q<- object$weights.q[Ti==0]

    names<- colnames(object$X)
    if(is.null(names)){
      p<- ncol(object$X)
      names<- paste("X",1:p,sep = "")
    }
    x1<- object$X[Ti==1,]
    x0<- object$X[Ti==0,]
    for(i in 1:ncol(x1)){
      if(i==2) par(ask = TRUE)

      if(length(unique(object$X[,i])) == 2){
        Treatment<- x1[,i]
        Placebo<- x0[,i]
        plot(c(0.5,1,2,2.5), c(2,mean(Treatment), mean(Placebo),2), pch = 16, cex = 1.5,
             ylim = c(0,1), ylab = "Mean of group", xlab = "",col = c("blue","red"),
             main = "Unweighted", xaxt = "n")
        axis(side = 1, at = c(1,2), labels = c("Treatment", "Control") )
        abline(h = mean(x1[,i]), lty = 2)

        new_control<- sum(w.q*Placebo)
        plot(c(0.5,1,2,2.5), c(2, mean(Treatment), new_control,2), pch = 16, cex = 1.5,
             ylim = c(0,1), ylab = "Mean of group", xlab = "",col = c("blue","red"),
             main = "Weighted", xaxt = "n")
        axis(side = 1, at = c(1,2), labels = c("Treatment", "Control") )
        abline(h = mean(x1[,i]), lty = 2)

      }else{
        rng<- range(c(x1[,i],x0[,i]))
        my.seq<- seq(rng[1],rng[2],length = 100)
        temp1<- sapply(my.seq, my.ecdf,x = x1[,i])
        temp0<- sapply(my.seq, my.ecdf,x = x0[,i])
        par(mfrow = c(1,2))
        plot(my.seq,temp1,cex = 0.4,pch = 16,type = "l",lty = 1,col = "red",
             xlab = names[i],ylab = "emperical CDF",main = "Unweighted emperical CDF")
        lines(my.seq, temp0, cex = 0.4, pch = 16, lty = 2,col = "blue")
        legend("bottomright", c("Treatment", "Control"), lty = c(1,2), col = c("red", "blue"))

        temp1<- sapply(my.seq, my.ecdf,x = x1[,i])
        temp0<- sapply(my.seq, my.ecdf,x = x0[,i], weights = w.q)
        plot(my.seq,temp1,cex = 0.4,pch = 16,type = "l",lty = 1,col = "red",
             xlab = names[i],ylab = "emperical CDF",main = "Weighted emperical CDF")
        lines(my.seq, temp0, cex = 0.4, pch = 16, lty = 2,col = "blue")
      }

    }
    par(ask = FALSE)
  }else{

    wgt <- object$weights.mat
    names<- colnames(object$X)
    if(is.null(names)){
      p<- ncol(object$X)
      names<- paste("X",1:p,sep = "")
    }
    p<- ncol(object$X)
    for(i in 1:p){
      if(i==2) par(ask = TRUE)

      if(length(unique(object$X[,i])) == 2){
        J<- object$J

        plot(c(0:(J-1)), c(mean(object$X[Ti==0, i]), rep(2,J-1)), pch = 16, cex = 1.5,
             ylim = c(0,1), ylab = "Mean of group", xlab = "Treatment group",col = 1,
             main = "Unweighted", xaxt = "n", xlim = c(-0.5,J-1+0.5))
        axis(side = 1, at = 0:(J-1), labels = paste(0:(J-1)) )
        abline(h = mean(object$X[,i]), lty = 2)
        for(j in 1:(object$J-1)){
          points(j, mean(object$X[Ti==j,i]), pch = 16, cex = 1.5, col = j+1)
        }

        plot(c(0:(J-1)), c(sum(object$X[Ti==0, i]*wgt[1,Ti==0]), rep(2,J-1)), pch = 16, cex = 1.5,
             ylim = c(0,1), ylab = "Mean of group", xlab = "Treatment group",col = 1,
             main = "Weighted", xaxt = "n",xlim = c(-0.5,J-1+0.5))
        axis(side = 1, at = 0:(J-1), labels = paste(0:(J-1)) )
        abline(h = mean(object$X[,i]), lty = 2)
        for(j in 1:(object$J-1)){
          points(j, sum(object$X[Ti==j,i]*wgt[j+1,Ti==j]), pch = 16, cex = 1.5, col = j+1)
        }

      }else{
        rng<- range(object$X[,i])
        my.seq<- seq(rng[1],rng[2],length = 100)
        temp0<- sapply(my.seq, my.ecdf,x = object$X[object$Ti==0,i])
        par(mfrow = c(1,2))
        plot(my.seq,temp0,cex = 0.4,pch = 16,type = "l",lty = 1,col = 1,
             xlab = names[i],ylab = "emperical CDF",main = "Unweighted emperical CDF")
        for(j in 1:(object$J-1)){
          temp<- sapply(my.seq, my.ecdf,x = object$X[object$Ti==j,i])
          lines(my.seq, temp, cex = 0.4, pch = 16, lty = j+1 , col = j+1)
        }
        legend("bottomright", paste("gp",0:(object$J-1)), lty = 1:object$J, col = 1:object$J )

        temp0<- sapply(my.seq, my.ecdf,x = object$X[object$Ti==0,i], weights = wgt[1,Ti==0])
        plot(my.seq,temp0,cex = 0.4,pch = 16,type = "l",lty = 1,col = 1,
             xlab = names[i],ylab = "emperical CDF",main = "Weighted emperical CDF")
        for(j in 1:(object$J-1)){
          temp<- sapply(my.seq, my.ecdf,x = object$X[object$Ti==j,i], weights = wgt[j+1,Ti==j])
          lines(my.seq, temp, cex = 0.4, pch = 16, lty = j+1 , col = j+1)
        }

      }

    }
    par(ask = FALSE)

  }
}

