


hmmcov=function(de1=NULL, ta1=NULL, de2=NULL, ta2=NULL, n, y, X,prop1, pi, Pi, m10=NULL, m20=NULL, 
glmtype="pois", maxitIAL=100, thresh=0, XE = NULL, maxitEM=50, 
EMconv=10^-6, glmconv=10^-5){

	library(MASS)

	if(is.null(XE)){
		XE = X
	}


	#for tuning parameter selection, one component is held at the fit of the full model, passed to function
	#if not doing tuning parameter selection, below options are ignored and parameters for both components are estimated
	#if(!is.null(de1) & !is.null(ta1) & !is.null(m20)) m2=m20	
	#if(!is.null(de2) & !is.null(ta2) & !is.null(m10)) m1=m1
	if(!is.null(m10) & !is.null(m20)){
		m1 = m10
		m2 = m20
	}else{
		probi1=(y<=quantile(y, prop1))^2
		pi1=mean(probi1)
	}

	#set number of states and dimension of X
	K=length(pi)
	p = ncol(X)
	p2 = ncol(XE)

	#create initial partition of data into backround and enriched based on count
	coef1=rep(0,p+1)
	coef2=rep(0,p2+1)

	#estimate inital beta for states that are not held fixed (only during tuning parameter selection)
	#if(!is.null(de1) & !is.null(ta1)){	
	#only fit on first full model fit, use m10 and m20 as starting always
	if(!any(is.null(c(de1, de2, ta1, ta2))) & is.null(m10)){	
		if(glmtype=="pois"){
			m1=glmIAL(y=y, X=scale(X), prior=probi1, family="poisson", prop=pi1, pMax=dim(X)[2], delta=de1, tau=ta1, nReEstimate=0,maxitIAL=maxitIAL, maxit=25, conv=glmconv)
		}else{
			m1=glmNB.IAL(y=y, X=scale(X), prior=probi1,   prop=pi1, pMax=dim(X)[2], delta=de1, tau=ta1, nReEstimate=0,maxitIAL=100, maxit=25, conv=glmconv)
		}
		vars1 =as.numeric(sqrt(apply(X, 2, var)))
		m1$b2use = m1$b2use/vars1[m1$w2use]
		m1$b0=log(m1$fitted/exp(X[,m1$w2use]%*%matrix(m1$b2use, length(m1$b2use), 1)))[1]	
		if(length(m1$w2use)>0) coef1[m1$w2use+1]=m1$b2use
		coef1[1]=m1$b0			
	}

	#if(!is.null(de2) & !is.null(ta2)){
	#only fit on full model fit, use m10 and m20 as starting always.  Assume Pi and pi are passed from first fit
	if(!any(is.null(c(de1, de2, ta1, ta2))) & is.null(m20)){	
		if(glmtype=="pois"){
			m2=glmIAL(y=y, X=scale(XE), prior=1-probi1, family="poisson", prop=1-pi1, pMax=dim(XE)[2], delta=de1, tau=ta1, nReEstimate=0,maxitIAL=maxitIAL, maxit=25, conv=glmconv)
		}else{
			m2=glmNB.IAL(y=y, X=scale(XE), prior=1-probi1,   prop=1-pi1, pMax=dim(XE)[2], delta=de2, tau=ta2, nReEstimate=0,maxitIAL=100, maxit=25, conv=glmconv)
		}
		vars2 =as.numeric(sqrt(apply(XE, 2, var)))
		m2$b2use = m2$b2use/vars2[m2$w2use]
		m2$b0=log(m2$fitted/exp(XE[,m2$w2use]%*%matrix(m2$b2use, length(m2$b2use), 1)))[1]
		if(length(m2$w2use)>0) coef2[m2$w2use+1]=m2$b2use
		coef2[1]=m2$b0
	}


	prob=matrix(0, n, K)
	if(glmtype=="pois"){
		prob[,1]=dpois(y, lambda=m1$fitted, log=T)
		prob[,2]=dpois(y, lambda=m2$fitted, log=T)
	}else{
		prob[,1]=dnbinom(y, mu=m1$fitted, size=1/m1$phi, log=T)
		prob[,2]=dnbinom(y, mu=m2$fitted, size=1/m2$phi, log=T)
	}

	resu=.C("forwardback", as.double(log(Pi)), as.double(log(pi)), as.double(prob), as.integer(nrow(prob)), as.integer(K), logalpha = double(nrow(prob)*K),logbeta = double(nrow(prob)*K),LL = double(1), package="hmmcov")
	forward = matrix(resu$logalpha, ncol=2) 
	backward = matrix(resu$logbeta, ncol=2) 
	forwardbackward = exp(forward + backward - resu$LL)	
		
	chsi2=matrix(0, n, K^2)			
	chsi2[-1,1]=forward[-n,1]+log(Pi[1,1])+prob[-1,1]+backward[-1,1]-resu$LL
	chsi2[-1,2]=forward[-n,2]+log(Pi[2,1])+prob[-1,1]+backward[-1,1]-resu$LL
	chsi2[-1,3]=forward[-n,1]+log(Pi[1,2])+prob[-1,2]+backward[-1,2]-resu$LL
	chsi2[-1,4]=forward[-n,2]+log(Pi[2,2])+prob[-1,2]+backward[-1,2]-resu$LL
	chsi2=exp(chsi2)

	c0 = rep(thresh, maxitEM)
	ll=rep(0, maxitEM)
	unstcoef1=unstcoef2=0
	phi1 = m1$phi
	phi2 = m2$phi
	a=c(Pi[1,], Pi[2,])

#begin ECM loop
	ptm <- proc.time()
	for(i in 1:maxitEM){
	  cat(".")
		#if c0>0, then utilizing rejection controlled ECM 
		if(c0[i]>0){
			#save weights
			c1 = forwardbackward[,1]
			c2 = forwardbackward[,2]

			#threshold given rejection procedure
			c1[c1<c0[i]] = rbinom(sum(c1<c0[i]), 1, c1[c1<c0[i]]/c0[i])*c0[i]  
			c2[c2<c0[i]] = rbinom(sum(c2<c0[i]), 1, c2[c2<c0[i]]/c0[i])*c0[i] 

			#normalized thresholded weights, save to f1, f2 
			f1 = c1/(c1+c2)
			f2 = c2/(c1+c2)

			#if both c1 and c2 are zero, set both to 0
			f1[is.na(f1)] = 0
			f2[is.na(f2)] = 0

			#remove observation whose new weights are 0
			XS1 = matrix(X[f1>0,],sum(f1>0), dim(X)[2]) 		
			y1 = y[f1>0]
			which1 = which(f1>0)
			XS2 = matrix(XE[f2>0,],sum(f2>0), dim(XE)[2])
			y2 = y[f2>0]
			which2 = which(f2>0)

			#if non meet threshold, set back to original data
			if(sum(forwardbackward[,1]<c0[i])==0){
				y1 =y
				XS1 = X
				which1 = 1:length(y)
			}	

			if(sum(forwardbackward[,2]<c0[i])==0){
				y2 = y
				XS2 = XE
				which2 = 1:length(y)
			}
		}else{
			f1 = forwardbackward[,1]
			y1 =y
			XS1 = X
			f2 = forwardbackward[,2]
			y2 = y
			XS2 = XE
			which1=which2 = 1:length(y)
		}

		#CM1 for beta using scaled covariate matrix
		if(!is.null(de1) & !is.null(ta1)){	
			if(glmtype=="pois"){
				m1=glmIAL(y=y1, X=scale(XS1), prior=f1[which1],   prop=mean(f1), pMax=dim(XS1)[2], delta=de1, tau=ta1, nReEstimate=0,maxitIAL=maxitIAL, maxit=25, conv=glmconv, fitted=m1$fitted[which1],family="poisson")
			}else{
				m1=glmNB.IAL(y=y1, X=scale(XS1), prior=f1[which1],   prop=mean(f1), pMax=dim(XS1)[2], delta=de1, tau=ta1,  nReEstimate=0, maxitIAL=maxitIAL, maxit=25, conv=glmconv, fitted=m1$fitted[which1], phi=phi1)
			}
			unstcoef1=m1$b2use	
			vars1 =as.numeric(sqrt(apply(XS1, 2, var)))
			m1$b2use = m1$b2use/vars1[m1$w2use]
			m1$b0=log(m1$fitted/exp(XS1[,m1$w2use]%*%matrix(m1$b2use, length(m1$b2use), 1)))[1]
			coef1=rep(0,p+1)
			if(length(m1$w2use)>0) coef1[m1$w2use+1]=m1$b2use
			coef1[1]=m1$b0	
	
			#CM2, dispersion parameters for the NB distribution
			if(glmtype == "nb"){
				#phi1 = 1/theta.ml(mu=m1$fitted, n=sum(f1[which1]), weights=f1[which1], y=y1, 100)	
				if(phi1 == 10^-5){
					phi1 = 1/theta.mm(dfr=sum(f1[which1])-p-2, mu=m1$fitted, weights=f1[which1], y=y1)
					if(phi1 < 10^-5) phi1 = .1
				}	
				phi1 = .C("phi_mlR", as.double(y1), as.double(m1$fitted), as.integer(length(y1)), as.integer(20), as.double(10^-6), as.double(phi1), as.integer(1), as.integer(0), as.double(f1[which1]), package="hmmcov")[[6]]
			}

			if(c0[i]>0){
				if(length(m1$w2use)>0){
					m1$fitted  = exp(as.numeric(X[,m1$w2use]%*%matrix(m1$b2use, m1$n2use, 1))+m1$b0)
				}else{
					m1$fitted  = rep(exp(m1$b0), length(y))
				}
		 }
		}

		if(!is.null(de2) & !is.null(ta2)){
			if(glmtype=="pois"){
				m2=glmIAL(y=y2, X=scale(XS2), prior=f2[which2],   prop=mean(f2), pMax=dim(XS2)[2], delta=de2, tau=ta2, nReEstimate=0,maxitIAL=maxitIAL, maxit=25, conv=glmconv, fitted=m2$fitted[which2],family="poisson")
			}else{
				m2=glmNB.IAL(y=y2, X=scale(XS2), prior=f2[which2],   prop=mean(f2), pMax=dim(XS2)[2], delta=de2, tau=ta2, nReEstimate=0,maxitIAL=maxitIAL, maxit=25, conv=glmconv, fitted=m2$fitted[which2], phi=phi2)
			}
			unstcoef2=m2$b2use
			vars2 =as.numeric(sqrt(apply(XS2, 2, var)))
			m2$b2use = m2$b2use/vars2[m2$w2use]
			m2$b0=log(m2$fitted/exp(XS2[,m2$w2use]%*%matrix(m2$b2use, length(m2$b2use), 1)))[1]
			coef2=rep(0,p2+1)
			if(length(m2$w2use)>0) coef2[m2$w2use+1]=m2$b2use
			coef2[1]=m2$b0	
			#CM2, dispersion parameters for the NB distribution
			if(glmtype == "nb"){
			#phi2 = 1/theta.ml(mu=m2$fitted, n=sum(f1[which2]), weights=f1[which2], y=y2, 100)
			if(phi2 == 10^-5){
				phi2 = 1/theta.mm(dfr=sum(f2[which2])-p-2, mu=m2$fitted, weights=f2[which2], y=y2)
				if(phi2 < 10^-5) phi2 = .1
			}		
				phi2 = .C("phi_mlR", as.double(y2), as.double(m2$fitted), as.integer(length(y2)), as.integer(20), as.double(10^-6), as.double(phi2), as.integer(1), as.integer(0), as.double(f2[which2]), package="hmmcov")[[6]]
			}

			if(c0[i]>0){
				if(length(m2$w2use)>0){
					m2$fitted  = exp(as.numeric(XE[,m2$w2use]%*%matrix(m2$b2use, m2$n2use, 1))+m2$b0)
				}else{
					m2$fitted  = rep(exp(m2$b0), length(y))
				}
			}
		}
		
		#E step
		#calculate emission probabilities for forward-backward algorithm on log scale

		prob=matrix(0, n, K)
		if(glmtype=="pois"){
			prob[,1]=dpois(y, lambda=m1$fitted, log=T)
			prob[,2]=dpois(y, lambda=m2$fitted, log=T)
		}else{
			prob[,1]=dnbinom(y, mu=m1$fitted, size=1/phi1, log=T)
			prob[,2]=dnbinom(y, mu=m2$fitted, size=1/phi2, log=T)
		}

		resu=.C("forwardback", as.double(log(Pi)), as.double(log(pi)), as.double(prob), as.integer(nrow(prob)), as.integer(K), logalpha = double(nrow(prob)*K),logbeta = double(nrow(prob)*K),LL = double(1), package="hmmcov")
		forward = matrix(resu$logalpha, ncol=2) 
		backward = matrix(resu$logbeta, ncol=2) 
		forwardbackward = exp(forward + backward - resu$LL)	
		ll[i] = resu$LL
		chsi2=matrix(0, n, K^2)			
		chsi2[-1,1]=forward[-n,1]+log(Pi[1,1])+prob[-1,1]+backward[-1,1]-resu$LL
		chsi2[-1,2]=forward[-n,2]+log(Pi[2,1])+prob[-1,1]+backward[-1,1]-resu$LL
		chsi2[-1,3]=forward[-n,1]+log(Pi[1,2])+prob[-1,2]+backward[-1,2]-resu$LL
		chsi2[-1,4]=forward[-n,2]+log(Pi[2,2])+prob[-1,2]+backward[-1,2]-resu$LL
		chsi2=exp(chsi2)

		#end of E-step, check stopping criterion
		if(i>1) if(abs((ll[i] - ll[i-1])/ll[i-1])<EMconv) break
		#end of E-step		 

		# calculate transition probability	
		for(i in 1:K){
			for(j in 1:K){
				Pi[i,j] = sum(chsi2[-1,(i-1)*K+j])/sum(forwardbackward[-1,i])
			}
		}
	 	#calculate state probabilities for first state
		pi=forwardbackward[1,] 
		a=c(Pi[1,], Pi[2,])
	}
  cat("\n")
	time = (proc.time()-ptm)

	m1final=rep(0,p)
	m1final[m1$w2use]=m1$b2use
	m2final=rep(0,p2)
	m2final[m2$w2use]=m2$b2use

	finalcoef=matrix(c(m1$b0, m1final, m2$b0,m2final, a, pi), 1, length(c(m1$b0, m1final, m2$b0,m2final, a, pi))) 
	if(glmtype=="nb") finalcoef=matrix(c(m1$b0, m1final, 1/phi1, m2$b0,m2final, 1/phi2, a, pi), 1, length(c(m1$b0, m1final, m2$b0,m2final, a, pi, 1/phi1, 1/phi2))) 

	m11=list()
	m11$fitted=m1$fitted
	m11$b0=m1$b0
	m11$b2use=m1$b2use
	m11$phi = phi1

  m22=list()
  m22$fitted=m2$fitted
  m22$b0=m2$b0
  m22$b2use=m2$b2use
  m22$phi	= phi2

	if((sum(c(is.null(de1), is.null(de2), is.null(ta1), is.null(ta2)))>0)){
		finalresults=list(y=0, X=0, forwardback=0, forward=0, resid=0, chsi2=0, Pi=Pi, pi=pi, coef=finalcoef, unstcoef1=unstcoef1,probi2=1, 
time=time, unstcoef2=unstcoef2, ll=ll[ll<0][sum(ll<0)], monotone=0, monotone2=0, BIC=-2*ll[ll<0][sum(ll<0)]+sum(finalcoef!=0)*log(n), de1=de1, de2 = de2 , ta1=ta1, ta2=ta2, m1=0, m2=0)		
	}else{
		finalresults=list(y=0, X=0, forwardback=forwardbackward, forward=0, resid=0, chsi2=0, Pi=Pi, pi=pi, coef=finalcoef, unstcoef1=unstcoef1,probi2=1, 
time=time, unstcoef2=unstcoef2, ll=max(ll[ll<0]), monotone=0, monotone2=0, BIC=-2*ll[ll<0][sum(ll<0)]+sum(finalcoef!=0)*log(n), de1=de1, de2 = de2 , ta1=ta1, ta2=ta2, m1=m11, m2=m22)		
	}
	return(finalresults)
}


hmm = function(y, X, prop1, maxitIAL=100,thresh=0.05, XE=NULL,maxitEM=50, glmtype="pois", EMconv=10^-6, glmconv=10^-5){
	
	if(length(y)!=nrow(X)) stop("length of y does not match rows of X, or X is not a matrix")

	n = length(y)
	t = quantile(y,prop1)
	Pi = matrix(0, 2, 2)
	Pi[1,1] = sum(y[-1] <= t & y[-n] <= t)/sum(y[-1]<=t)
	Pi[1,2] = sum(y[-1] <= t & y[-n] >  t)/sum(y[-1]<=t)
	Pi[2,1] = sum(y[-1] > t & y[-n] <= t)/sum(y[-1]>t)
	Pi[2,2] = sum(y[-1] > t & y[-n] > t)/sum(y[-1]>t)
	
	pi = abs(c(y[1]<=t, y[1]>t) - 10^-100)

	res = hmmcov(de1=0, ta1=-1, de2=0, ta2=-1, n=n, y=y, X=X, prop1=prop1, pi=pi, Pi=Pi,m10=NULL, m20=NULL, 
		glmtype=glmtype, maxitIAL=maxitIAL, thresh=thresh, XE = XE, maxitEM=maxitEM, EMconv=EMconv, glmconv=glmconv)

	coef  = matrix(res$coef[1:(length(res$coef)-6)], ncol=2)
	rownames(coef) = rep("", nrow(coef))

	if(glmtype == "nb") rownames(coef)[nrow(coef)] = "phi" 
	final = list(forwardbackward = res$forwardback, coef = coef, ll=res$ll, Pi = res$Pi, pi=res$pi)
	return(final)
}
