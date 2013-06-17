glmIAL <-
function(family, y,  X, delta=seq(0, 5.0, by=0.5),  
         tau = c(0.01*c(10,8,6,4,2,1), 0.0001*c(50,20,10,5,2,1)),
         pMax=20, offset=NULL, naPercent=0.4, 
         nTotal=NULL, maxit=50, maxitIAL=20, nReEstimate=pMax, 
         conv=1e-5, fitted=NULL, trace=1, prior=NULL, prop=NULL, 
protect=0)
{
  if(!is.numeric(y)){
    stop("y must be a numeric vector\n")
  }
  
  if(!is.matrix(X)){
    stop("X must be a vector\n")
  }

  if(is.null(prior)) prior=rep(1,length(y))
  
  if(is.null(prop)) prop=1

  M = ncol(X)
  N = length(y)
  
  if(nrow(X) != N){
    stop("the dimension of X and y do not match\n")
  }
  
  if(!is.null(offset)){
    useOffset = 1
    if((! is.numeric(offset)) || length(offset) != N){
      stop("offset must be a numeric vector of the same length as y\n")
    }
  }else{
    useOffset = 0
    offset    = 0.0 
  }
  
  if(family == "binomial"){
    if((! is.numeric(nTotal)) || length(nTotal) != N){
      strWarn = "For binomial family, nTotal must be a numeric vector"
      stop(strWarn, "of the same length as y\n")
    }
  }
  
  if(is.null(fitted)){
    init = 0
    fitted = rep(0, N)
  }else{
    init = 1
    if((! is.numeric(fitted)) || length(fitted) != N){
      stop("fitted must be a numeric vector of the same length as y\n")
    }
  }
  
  #define BINOMIAL  1
  #define POISSON   2
  #define GAMMA     4
  
  #define LOGIT     1
  #define LOG       2
  #define INVERSE   4
  
  if(family=="binomial"){
    familyR = 1
    linkR   = 1
  }else if(family=="poisson"){
    familyR = 2
    linkR   = 2
  }else if(family=="gamma"){
    familyR = 4
    linkR   = 4
  }else if(family=="gaussian"){
    familyR = 3
    linkR   = 3
  }else{
    stop("invalid family\n")
  }
  
  isNA = apply(X, 1, function(v){ any(is.na(v)) })
  isNA = is.na(y) | isNA
  
  if(length(which(isNA))/N > naPercent){
    stop("percent of missing data is too high\n")
  }
  
  w2kp  = which(!isNA)
  yR    = y[w2kp]
  XR    = X[w2kp,]
  wXR   = XR
  nXR   = matrix(0, nrow=N, ncol=max(nReEstimate+1,1))
  
  if(useOffset){
    offset = offset[w2kp]
  }
  if(family == "binomial"){
    nTotal = nTotal[w2kp]
  }

  N  = length(yR)
  
  protect=protect-1

  dims    = numeric(11)
  dims[1] = N
  dims[2] = M
  dims[3] = maxit
  dims[4] = maxitIAL
  dims[5] = init
  dims[6] = useOffset
  dims[7] = length(delta)
  dims[8] = length(tau)
  dims[9] = pMax
  dims[10] = nReEstimate
  dims[11]= protect
  
  resid   = weights = numeric(N)
  fitted  = fitted[w2kp]
  nIter   = 0
  phi     = 0.0
  
  w2use   = rep(-9, pMax)
  b2use   = numeric(pMax)
  score   = matrix(0, nrow=length(delta), ncol=length(tau))
  b0      = matrix(0, nrow=length(delta), ncol=length(tau))
  b02use  = 0.0
  convg   = 0
  n2use   = 0
  score2use = 0.0
  delta2use = 0.0
  tau2use   = 0.0
  wss  = numeric(maxit)
  wss0 = 0
  
  Z = .C("glmIAL", as.integer(familyR), as.integer(linkR), 
         as.integer(dims), nIter=as.integer(nIter), 
         as.double(yR), as.double(offset), as.double(XR), 
         as.double(wXR), as.double(nXR), as.double(nTotal), 
         as.double(conv), fitted=as.double(fitted), as.double(resid), 
         as.double(weights), as.double(phi), as.integer(trace), 
         n2use=as.integer(n2use), w2use=as.integer(w2use), 
         b2use=as.double(b2use), b02use=as.double(b02use),
         b0=as.double(b0), score=as.double(score), 
         score2use=as.double(score2use), as.double(delta), 
         delta2use=as.double(delta2use), as.double(tau), 
         tau2use=as.double(tau2use), wss=as.double(wss),
         wss0=as.double(wss0), as.double(prior), as.double(prop),PACKAGE="hmmcov")
  
  n2use=Z$n2use
  
  if(n2use > 0){
    w2use = Z$w2use[1:n2use]+1
    b2use = Z$b2use[1:n2use]
  }else{
    w2use = b2use = NULL
  }
  
  if (Z$nIter >= maxit){
    warning(sprintf("reach mamximum %d iterations\n", maxit))  
  }
  
  score  = matrix(Z$score, nrow=length(delta), ncol=length(tau), byrow=TRUE)
  b0     = matrix(Z$b0,    nrow=length(delta), ncol=length(tau), byrow=TRUE)
  wss    = Z$wss[1:Z$nIter]

  result = list(n2use=n2use, w2use=w2use, b2use=b2use, b0=b0, 
                score2use=Z$score2use, delta2use=Z$delta2use, 
                tau2use=Z$tau2use, score=score, nIter=Z$nIter, 
                wss=wss, wss0=Z$wss0, fitted=Z$fitted)
  
  result
}
