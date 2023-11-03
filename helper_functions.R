###########################################################
### the followings are functions called in onestep SCAD ###
### don't need to directly use them in the computation  ###
###########################################################

dp<-function(theta,lambda,a=3.7){
  p<-length(theta)
  b1<-rep(0,p)
  b1[theta>lambda]<-1
  b2<-rep(0,p)
  b2[theta<(lambda*a)]<-1
  lambda*(1-b1)+((lambda*a)-theta)*b2/(a-1)*b1
}

library("lars")

mylasso<-function (x, y, xw,type = c("lasso", "lar", "forward.stagewise"),
                   trace = FALSE, Gram, eps = .Machine$double.eps, max.steps,
                   use.Gram = TRUE)
{   
  call <- match.call()
  type <- match.arg(type)
  TYPE <- switch(type, lasso = "LASSO", lar = "LAR", forward.stagewise = "Forward Stagewise")
  if (trace)
    cat(paste(TYPE, "sequence\n"))
  nm <- dim(x)
  n <- nm[1]
  m <- nm[2]
  im <- inactive <- seq(m)
  one <- rep(1, n)
  vn <- dimnames(x)[[2]]
  meanx<-rep(0,m)
  x <- scale(x, meanx, FALSE)
  normx<-xw
  nosignal <- normx/sqrt(n) < eps
  if (any(nosignal)) {
    ignores <- im[nosignal]
    inactive <- im[-ignores]
    normx[nosignal] <- eps * sqrt(n)
    if (trace)
      cat("LARS Step 0 :  ", sum(nosignal), "Variables with Variance < eps; dropped for good\n")
  }
  else ignores <- NULL
  names(normx) <- NULL
  x <- scale(x, FALSE, normx)
  
  if (use.Gram & missing(Gram)) {
    if (m > 500 && n < m)
      cat("There are more than 500 variables and n<m;\nYou may wish to restart and set use.Gram=FALSE\n")
    if (trace)
      cat("Computing X'X .....\n")
    Gram <- t(x) %*% x
  }
  mu <- 0
  y <- drop(y - mu)
  Cvec <- drop(t(y) %*% x)
  ssy <- sum(y^2)
  residuals <- y
  if (missing(max.steps))
    max.steps <- 50 * min(m, n - 1)
  beta <- matrix(0, max.steps + 1, m)
  penalty<-max(abs(Cvec))
  Gamrat <- NULL
  arc.length <- NULL
  R2 <- 1
  RSS <- ssy
  first.in <- integer(m)
  active <- NULL
  actions <- as.list(seq(max.steps))
  drops <- FALSE
  Sign <- NULL
  R <- NULL
  k <- 0
  maxvar<-m
  while ((k < max.steps) & (length(active) < maxvar)) {
    action <- NULL
    k <- k + 1
    C <- Cvec[inactive]
    Cmax <- max(abs(C))
    if (!any(drops)) {
      new <- abs(C) >= Cmax - eps
      C <- C[!new]
      new <- inactive[new]
      for (inew in new) {
        if (use.Gram) {
          R <- updateR(Gram[inew, inew], R, drop(Gram[inew,
                                                      active]), Gram = TRUE, eps = eps)
        }
        else {
          R <- updateR(x[, inew], R, x[, active], Gram = FALSE,
                       eps = eps)
        }
        if (attr(R, "rank") == length(active)) {
          nR <- seq(length(active))
          R <- R[nR, nR, drop = FALSE]
          attr(R, "rank") <- length(active)
          ignores <- c(ignores, inew)
          action <- c(action, -inew)
          if (trace)
            cat("LARS Step", k, ":       Variable", inew,
                " collinear; dropped for good\n")
        }
        else {
          if (first.in[inew] == 0)
            first.in[inew] <- k
          active <- c(active, inew)
          Sign <- c(Sign, sign(Cvec[inew]))
          action <- c(action, inew)
          if (trace)
            cat("LARS Step", k, ":       Variable", inew,
                " added\n")
        }
      }
    }
    else action <- -dropid
    Gi1 <- backsolve(R, backsolvet(R, Sign))
    dropouts <- NULL
    if (type == "forward.stagewise") {
      directions <- Gi1 * Sign
      if (!all(directions > 0)) {
        if (use.Gram) {
          nnls.object <- nnls.lars(active, Sign, R, directions,
                                   Gram[active, active], trace = trace, use.Gram = TRUE,
                                   eps = eps)
        }
        else {
          nnls.object <- nnls.lars(active, Sign, R, directions,
                                   x[, active], trace = trace, use.Gram = FALSE,
                                   eps = eps)
        }
        positive <- nnls.object$positive
        dropouts <- active[-positive]
        action <- c(action, -dropouts)
        active <- nnls.object$active
        Sign <- Sign[positive]
        Gi1 <- nnls.object$beta[positive] * Sign
        R <- nnls.object$R
        C <- Cvec[-c(active, ignores)]
      }
    }
    A <- 1/sqrt(sum(Gi1 * Sign))
    w <- A * Gi1
    if (!use.Gram)
      u <- drop(x[, active, drop = FALSE] %*% w)
    if (length(active) >= maxvar) {
      gamhat <- Cmax/A
    }
    else {
      if (use.Gram) {
        a <- drop(w %*% Gram[active, -c(active, ignores),
                             drop = FALSE])
      }
      else {
        a <- drop(u %*% x[, -c(active, ignores), drop = FALSE])
      }
      gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
      gamhat <- min(gam[gam > eps], Cmax/A)
    }
    if (type == "lasso") {
      dropid <- NULL
      b1 <- beta[k, active]
      z1 <- -b1/w
      zmin <- min(z1[z1 > eps], gamhat)
      if (zmin < gamhat) {
        gamhat <- zmin
        drops <- z1 == zmin
      }
      else drops <- FALSE
    }
    beta[k + 1, ] <- beta[k, ]
    beta[k + 1, active] <- beta[k + 1, active] + gamhat *
      w
    if (use.Gram) {
      Cvec <- Cvec - gamhat * Gram[, active, drop = FALSE] %*%
        w
    }
    else {
      residuals <- residuals - gamhat * u
      Cvec <- drop(t(residuals) %*% x)
    }
    
    #print(penalty[k]-abs(gamhat*A)-max(abs(Cvec)))
    
    penalty <- c(penalty,max(abs(Cvec)))
    
    #max(Cvec)
    
    Gamrat <- c(Gamrat, gamhat/(Cmax/A))
    arc.length <- c(arc.length, gamhat)
    
    
    if (type == "lasso" && any(drops)) {
      dropid <- seq(drops)[drops]
      for (id in rev(dropid)) {
        if (trace)
          cat("Lasso Step", k + 1, ":    Variable", active[id],
              "   dropped\n")
        R <- downdateR(R, id)
      }
      dropid <- active[drops]
      beta[k + 1, dropid] <- 0
      active <- active[!drops]
      Sign <- Sign[!drops]
    }
    if (!is.null(vn))
      names(action) <- vn[abs(action)]
    actions[[k]] <- action
    inactive <- im[-c(active, ignores)]
  }
  beta <- beta[seq(k + 1), ]
  dimnames(beta) <- list(paste(0:k), vn)
  if (trace)
    cat("Computing residuals, RSS etc .....\n")
  residuals <- y - x %*% t(beta)
  beta <- scale(beta, FALSE, normx)
  RSS <- apply(residuals^2, 2, sum)
  R2 <- 1 - RSS/RSS[1]
  Cp <- ((n - k - 1) * RSS)/rev(RSS)[1] - n + 2 * seq(k + 1)
  object <- list(call = call, type = TYPE, R2 = R2, RSS = RSS,
                 Cp = Cp, actions = actions[seq(k)], entry = first.in,
                 Gamrat = Gamrat, arc.length = arc.length, Gram = if (use.Gram) Gram else NULL,
                 beta = beta, mu = mu, normx = normx, meanx = meanx,penalty=penalty)
  class(object) <- "lars"
  object
}


predict.mylasso<-function (object, newx, s, mode = c("penalty", "fraction", "step"))
{
  mode <- match.arg(mode)
  betas <- object$beta
  sbetas <- scale(betas, FALSE, 1/object$normx)
  kp <- dim(betas)
  k <- kp[1]
  p <- kp[2]
  steps <- seq(k)
  
  sbeta <- switch(mode, fraction = {
    if (any(s > 1) | any(s < 0))
      stop("Argument s out of range")
    nbeta <- drop(abs(sbetas) %*% rep(1, p))
    nbeta/nbeta[k]
  }, penalty = {
    pen <- object$penalty
    #if (any(s > pen[1]) | any(s < 0))
    #    stop("Argument s out of range")
    s[s>pen[1]]<-pen[1]
    s[s<0]<-0
    if (any(s > pen[1]) | any(s < 0))
      stop("Argument s out of range"
      )
    pen
  }, step = {
    if (any(s < 0) | any(s > k))
      stop("Argument s out of range")
    steps
  }
  )
  sfrac <- (s - sbeta[1])/(sbeta[k] - sbeta[1])
  sbeta <- (sbeta - sbeta[1])/(sbeta[k] - sbeta[1])
  usbeta <- unique(sbeta)
  useq <- match(usbeta, sbeta)
  sbeta <- sbeta[useq]
  betas <- betas[useq, ]
  coord <- approx(sbeta, seq(sbeta), sfrac)$y
  left <- floor(coord)
  right <- ceiling(coord)
  newbetas <- ((sbeta[right] - sfrac) * betas[left, , drop = FALSE] +
                 (sfrac - sbeta[left]) * betas[right, , drop = FALSE])/(sbeta[right] -
                                                                          sbeta[left])
  newbetas[left == right, ] <- betas[left[left == right], ]
  
  y<-object$y
  d<-object$d
  x<-object$x
  
  
  robject <- list(s = s,  mode = mode, coefficients = drop(newbetas),fit = drop(newx%*%t(newbetas)))
  robject
}


d1loglik<-function(beta,y,x,family=c("binomial","poisson")){
  family <- match.arg(family)
  n<-length(y)
  p<-dim(x)[2]
  u<-drop(x%*%beta)
  r <- switch(family,binomial = {
    u1<-exp(u)/(1+exp(u)) 
    d1<-t(x)%*%(u1-y)
    d1
  }, poisson = {
    u2<-exp(u)
    d1<-t(x)%*%(u2-y)
    d1
  })
  r
}

d2loglik<-function(beta,y,x,family=c("binomial","poisson")){
  family <- match.arg(family)
  n<-length(y)
  p<-dim(x)[2]
  u<-drop(x%*%beta)
  r <- switch(family, binomial = {
    w1<-exp(u)/((1+exp(u))^2)
    #cat("d2loglink",w1,"\n")
    if(n>1){w<-diag(w1)}
    else{w<-matrix(w1,1,1)}
    d2<-t(x)%*%w%*%x
    d2   
  }, poisson = {
    w2<-exp(u)
    if(n>1){w<-diag(w2)}
    else{w<-matrix(w2,2,1)}
    d2<-t(x)%*%w%*%x
    d2 
  })
  r
}



devresd<-function(y,fit,family=c("binomial","poisson")){
  family <- match.arg(family)
  z<-exp(fit)
  r <- switch(family, poisson = {
    u<-z
    r<-y
    id1<-(1:length(y))[y>0]
    id0<-(1:length(y))[y==0]
    if (length(id1)==length(y)){
      r<-sqrt(2*(y*log(y)-y*fit-y+u))*sign(y-u)
    }
    else{
      if(length(id0)==length(y)){ r<-sqrt(2*(u))*sign(-u) }
      else{ 
        r[id1]<-sqrt(2*(y[id1]*log(y[id1])-y[id1]*fit[id1]-y[id1]+u[id1]))*sign(y[id1]-u[id1])
        r[id0]<-sqrt(2*(u[id0]))*sign(-u[id0])
      }
    }
  }, binomial = {
    u<-z/(1+z)
    r<-sqrt(-y*fit+log(1+z))*sign(y-u)
  })
  r
}