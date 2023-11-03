#######################################
### copyright: Hui Zou and Runze Li ###     
#######################################

########################################
### The linear regression model case ###
########################################

scad.onestep<-function(beta0,y,x,lambda){
  #the main function, compute the one-step SCAD estimator
  #x is the predictor matrix
  #y is repsonse vector
  #lambda is the SCAD penalization parameter
  #beta0 is an initial estimator, for example, we used LS-estimator in the simulation. 
  call <- match.call()
  if (missing(beta0)){
    fit0<-lm(y~x-1)
    beta0<-fit0$coef
  }
  w<-dp(abs(beta0),lambda=lambda)*2*length(y)
  if(sum(w!=0)==length(w)){
    object<-mylasso(x=x,y=y,xw=w)
    beta<-predict.mylasso(object,newx=x,s=lambda,mode="penalty")$coef
  }
  else{
    if(sum(w==0)==length(w)){
      beta<-beta0
    }
    else{
      ids<-(1:length(w))
      id0<-ids[w==0]
      id1<-ids[w!=0]
      x1<-x[,id0,drop=FALSE]
      H1<-x1%*%solve(t(x1)%*%x1)%*%t(x1)
      y1<-y-drop(H1%*%y)
      x2<-x[,id1,drop=FALSE]-H1%*%x[,id1,drop=FALSE]
      if(length(id1)>1){
        object<-mylasso(x=x2,y=y1,xw=w[id1])
        beta2<-predict.mylasso(object,newx=x2,s=lambda,mode="penalty")$coef
      }
      else{
        aa<-t(x2)%*%y1/(t(x2)%*%x2)
        bb<-w[id1]/(2*t(x2)%*%x2)
        if (abs(aa)>bb){beta2<-(abs(aa)-bb)*sign(aa)}
        else{beta2<-0}
      }
      beta1<-drop(solve(t(x1)%*%x1)%*%t(x1)%*%(y-x[,id1,drop=FALSE]%*%beta2))
      beta<-beta0
      beta[id0]<-beta1
      beta[id1]<-beta2
    }
  }
  beta
}



cv.scad.onestep<-function(x,y,K=5,s){
  #compute cross-validation error 
  #x represents the predictor matrix
  #y is response vector
  #K is the number of folds in CV
  #s is the tuning parameter in one-step SCAD
  
  all.folds <- cv.folds(length(y), K)
  residmat <- matrix(0, length(s), K)
  for (i in seq(K)) {
    omit <- all.folds[[i]]
    for (j in 1:length(s)){
      
      beta.fit <- scad.onestep(x=x[-omit, ], y=y[-omit], lambda=s[j])
      fit <- drop(x[omit, , drop = FALSE]%*%beta.fit)
      residmat[j, i] <- mean((y[omit] - fit)^2)
      #cat(i,j,"\n")
    }
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object <- list(s = s, cv = cv, cv.error = cv.error)
  invisible(object)
}


