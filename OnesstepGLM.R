#########################################
### The generalized linear model case ###
#########################################

scad.onestep.glm<-function(beta0,y,x,family=c("binomial","poisson"),lambda){
  call <- match.call()
  family <- match.arg(family)
  if (missing(beta0)){
    fit0<-glm(y~x-1,family=family)
    beta0<-fit0$coef
  }
  Sigma<-d2loglik(beta=beta0,y=y,x=x,family=family)
  V<-chol(Sigma) 
  u<-V%*%beta0
  w<-dp(abs(beta0),lambda=lambda)*2*length(y)/lambda
  if(sum(w!=0)==length(w)){
    object<-mylasso(x=V,y=u,xw=w)
    beta<-predict.mylasso(object,newx=V,s=lambda,mode="penalty")$coef
  }
  else{
    if(sum(w==0)==length(w)){
      beta<-beta0
    }  
    else{
      ids<-(1:length(w))
      id0<-ids[w==0]
      id1<-ids[w!=0] 
      x1<-V[,id0,drop=FALSE]
      x2<-V[,id1,drop=FALSE]
      H1<-x1%*%solve(t(x1)%*%x1)%*%t(x1)
      y1<-u-drop(H1%*%u)
      x2<-V[,id1,drop=FALSE]-H1%*%V[,id1,drop=FALSE]
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
      beta1<-drop(solve(t(x1)%*%x1)%*%t(x1)%*%(u-V[,id1,drop=FALSE]%*%beta2))
      beta<-beta0
      beta[id0]<-beta1
      beta[id1]<-beta2
    }
  }  
  beta
}

cv.scad.onestep.glm<-function(x,y,K=5,s,family=c("binomial","poisson")){
  family <- match.arg(family)
  all.folds <- cv.folds(length(y), K)
  residmat <- matrix(0, length(s), K)
  for (i in seq(K)) {
    omit <- all.folds[[i]]
    for (j in 1:length(s)){
      beta.fit <- scad.onestep.glm(x=x[-omit, ], y=y[-omit], lambda=s[j],family=family)
      fit <- drop(x[omit, , drop = FALSE]%*%beta.fit)
      residmat[j, i] <- mean((devresd(y[omit],fit,family=family))^2)#sum of square of deviance residuals
    }
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object <- list(s = s, cv = cv, cv.error = cv.error)
  invisible(object)
}
