
augvec=function(y, K=2){
	arp=1
	m <- matrix(1:(K*length(y)), ncol = K, byrow = TRUE)
	yaug = matrix(y, length(y)*K, arp, byrow=F)
	yaug <-yaug[order(m)]
	return(yaug)
}

arhmmcov=function(de1=NULL, ta1=NULL, de2=NULL, ta2=NULL, n, y, X,prop1, pi, Pi, m10=NULL, m20=NULL, 
glmtype="pois", maxitIAL=100, thresh=0, XE = NULL, maxitEM=50, 
EMconv=10^-6, glmconv=10^-5){

	#print(c(de1, de2, thresh, maxitIAL))

	if(is.null(XE)){
		XE = X
	}

	m <- matrix(1:(nrow(X)*2), ncol = 2, byrow = TRUE)
	Xl <- rbind(X, X)[order(m), ]
	Xaug=Xl

	m <- matrix(1:(nrow(XE)*2), ncol = 2, byrow = TRUE)
	Xl <- rbind(XE, XE)[order(m), ]
	XEaug=Xl

	yaug = augvec(y)

	#for tuning parameter selection, one component is held at the fit of the full model, passed to function
	#if not doing tuning parameter selection, below options are ignored and parameters for both components are estimated
	#if(!is.null(de1) & !is.null(ta1) & !is.null(m20)) m2=m20	
	#if(!is.null(de2) & !is.null(ta2) & !is.null(m10)) m1=m1
	if(!is.null(m10) & !is.null(m20)){
		m1 = m10
		m2 = m20
		m1$w2use = 1:length(m1$b2use)
		m2$w2use = 1:length(m2$b2use)
	}else{
		probi1=(yaug<=quantile(yaug, prop1))^2
		pi1=mean(probi1)
	}

	#set number of states and dimension of X
	K=length(pi)
	p = ncol(Xaug)
	p2 = ncol(XEaug)

	coef1=rep(0,p+1)
	coef2=rep(0,p2+1)

	#estimate inital beta for states that are not held fixed (only during tuning parameter selection)
	if(!is.null(de1) & !is.null(ta1) & is.null(m10)){	
		if(glmtype=="pois"){
			m1=glmIAL(y=yaug, X=scale(Xaug), prior=probi1, family="poisson", prop=pi1, pMax=dim(Xaug)[2], delta=de1, tau=ta1, nReEstimate=0,maxitIAL=maxitIAL, maxit=25, conv=glmconv)
		}else{
			m1=glmNB.IAL(y=yaug, X=scale(Xaug), prior=probi1,   prop=pi1, pMax=dim(Xaug)[2], delta=de1, tau=ta1, nReEstimate=0,maxitIAL=100, maxit=25, conv=glmconv)
		}
		vars1 =as.numeric(sqrt(apply(Xaug, 2, var)))
		m1$b2use = m1$b2use/vars1[m1$w2use]
		m1$b0=log(m1$fitted/exp(Xaug[,m1$w2use]%*%matrix(m1$b2use, length(m1$b2use), 1)))[1]	
		if(length(m1$w2use)>0) coef1[m1$w2use+1]=m1$b2use
		coef1[1]=m1$b0			
	}

	if(!is.null(de2) & !is.null(ta2) & is.null(m20)){
		if(glmtype=="pois"){
			m2=glmIAL(y=yaug, X=scale(XEaug), prior=1-probi1, family="poisson", prop=1-pi1, pMax=dim(XEaug)[2], delta=de1, tau=ta1, nReEstimate=0,maxitIAL=maxitIAL, maxit=25, conv=glmconv)
		}else{
			m2=glmNB.IAL(y=yaug, X=scale(XEaug), prior=1-probi1,   prop=1-pi1, pMax=dim(XEaug)[2], delta=de2, tau=ta2, nReEstimate=0,maxitIAL=100, maxit=25, conv=glmconv)
		}
		vars2 =as.numeric(sqrt(apply(XEaug, 2, var)))
		m2$b2use = m2$b2use/vars2[m2$w2use]
		m2$b0=log(m2$fitted/exp(XEaug[,m2$w2use]%*%matrix(m2$b2use, length(m2$b2use), 1)))[1]
		if(length(m2$w2use)>0) coef2[m2$w2use+1]=m2$b2use
		coef2[1]=m2$b0
	}

	prob=matrix(0, n, K^2)
	if(glmtype!="nb"){
		prob[,1]=dpois(y, lambda=m1$fitted[seq(1, 2*n, by=2)], log=T)
		prob[,2]=dpois(y, lambda=m1$fitted[seq(2, 2*n, by=2)], log=T)
		prob[,3]=dpois(y, lambda=m2$fitted[seq(1, 2*n, by=2)], log=T)
		prob[,4]=dpois(y, lambda=m2$fitted[seq(2, 2*n, by=2)], log=T)
	}else{
		prob[,1]=dnbinom(y, mu=m1$fitted[seq(1, 2*n, by=2)], size=1/m1$phi, log=T)
		prob[,2]=dnbinom(y, mu=m1$fitted[seq(2, 2*n, by=2)], size=1/m1$phi, log=T)
		prob[,3]=dnbinom(y, mu=m2$fitted[seq(1, 2*n, by=2)], size=1/m2$phi, log=T)
		prob[,4]=dnbinom(y, mu=m2$fitted[seq(2, 2*n, by=2)], size=1/m2$phi, log=T)
	}
	resu=.C("forwardbackaug", as.double(log(Pi)), as.double(log(pi)), as.double(prob), as.integer(nrow(prob)), as.integer(K), logalpha = double(nrow(prob)*K),logbeta = double(nrow(prob)*K),LL = double(1), package="hmmcov")
	forward = matrix(resu$logalpha, ncol=2) 
	backward = matrix(resu$logbeta, ncol=2) 
	forwardbackward = exp(forward + backward - resu$LL)	
		
	chsi2=matrix(0, n, K^2)			
	chsi2[-1,1]=forward[-n,1]+log(Pi[1,1])+prob[-1,1]+backward[-1,1]-resu$LL
	chsi2[-1,2]=forward[-n,2]+log(Pi[2,1])+prob[-1,2]+backward[-1,1]-resu$LL
	chsi2[-1,3]=forward[-n,1]+log(Pi[1,2])+prob[-1,3]+backward[-1,2]-resu$LL
	chsi2[-1,4]=forward[-n,2]+log(Pi[2,2])+prob[-1,4]+backward[-1,2]-resu$LL
	chsi2=exp(chsi2)

	chsiaug=matrix(0, n*K,K)	
	chsiaug[seq(1, 2*n, by=2),1] = chsi2[,1]
	chsiaug[seq(2, 2*n, by=2),1] = chsi2[,2]
	chsiaug[seq(1, 2*n, by=2),2] = chsi2[,3]
	chsiaug[seq(2, 2*n, by=2),2] = chsi2[,4]

	chsiaug[1,] = c(0,0)
	chsiaug[2,] = forwardbackward[1,]


	res=res0=rep(0,2*n)
	#res0[seq(1, 2*n, by=2)]=(log(y + 1) - (X%*%matrix(m1$b2use, m1$n2use, 1)+m1$b0))
  #res0[seq(2, 2*n, by=2)]=(log(y + 1) - (XE%*%matrix(m2$b2use,m2$n2use, 1)+m2$b0))
	if(length(m1$b2use)>0){
		res0[seq(1, 2*n, by=2)]=(log(y + 1) - log(exp(X%*%matrix(m1$b2use[1:p], p, 1)+m1$b0)+1))
	}else{
		res0[seq(1, 2*n, by=2)]=(log(y + 1) - log(exp(m1$b0)+1))		
	}

	if(length(m2$b2use)>0){
	  res0[seq(1, 2*n, by=2)]=(log(y + 1) - log(exp(XE%*%matrix(m2$b2use[1:p2], p2, 1)+m2$b0)+1))
	}else{
		res0[seq(1, 2*n, by=2)]=(log(y + 1) - log(exp(m2$b0)+1))
	}

	res=c(0,0,res0[-((2*n-1):(2*n))])
	Xaug=cbind(Xaug, res)
	XEaug=cbind(XEaug, res)
	p=p+1
	p2=p2+1

	a=c(Pi[1,], Pi[2,])

	c0 = rep(thresh, maxitEM)
	ll=rep(0, maxitEM)
	unstcoef1=unstcoef2=0
	phi1 = m1$phi
	phi2 = m2$phi

	library(MASS)
	#begin ECM loop
	ptm <- proc.time()
	for(i in 1:maxitEM){
		#if(i==5) break
		#if c0>0, then utilizing rejection controlled ECM 
		if(c0[i]>0){
			#save weights
			chsi2[1,] = c(0, forwardbackward[1,1], 0 , forwardbackward[1,2])
	
			c1 = chsi2[,1]
			c2 = chsi2[,2]
			c3 = chsi2[,3]
			c4 = chsi2[,4]

			#threshold given rejection procedure
			c1[c1<c0[i]] = rbinom(sum(c1<c0[i]), 1, c1[c1<c0[i]]/c0[i])*c0[i]  
			c2[c2<c0[i]] = rbinom(sum(c2<c0[i]), 1, c2[c2<c0[i]]/c0[i])*c0[i] 
			c3[c3<c0[i]] = rbinom(sum(c3<c0[i]), 1, c3[c3<c0[i]]/c0[i])*c0[i]  
			c4[c4<c0[i]] = rbinom(sum(c4<c0[i]), 1, c4[c4<c0[i]]/c0[i])*c0[i] 

			#normalized thresholded weights, save to f1, f2 
			f1 = c1/(c1+c2+c3+c4)
			f2 = c2/(c1+c2+c3+c4)
			f3 = c3/(c1+c2+c3+c4)
			f4 = c4/(c1+c2+c3+c4)

			#if c1-c4 are all zero, set f1-f4 to 0
			f1[is.na(f1)] = 0
			f2[is.na(f2)] = 0
			f3[is.na(f3)] = 0
			f4[is.na(f4)] = 0

			ff=matrix(0, n*K,K)	
			ff[seq(1, 2*n, by=2),1] = f1
			ff[seq(2, 2*n, by=2),1] = f2
			ff[seq(1, 2*n, by=2),2] = f3
			ff[seq(2, 2*n, by=2),2] = f4

			f1 = ff[,1]
			f2 = ff[,2]

			#remove observation whose new weights are 0
			XS1 = matrix(Xaug[f1>0,],sum(f1>0), dim(Xaug)[2]) 		
			y1 = yaug[f1>0]
			which1 = which(f1>0)
			XS2 = matrix(XEaug[f2>0,],sum(f2>0), dim(XEaug)[2])
			y2 = yaug[f2>0]
			which2 = which(f2>0)

			#if non meet threshold, set back to original data
			if(sum(chsiaug[,1]<c0[i])==0){
				y1 =yaug
				XS1 = Xaug
				which1 = 1:length(yaug)
			}	

			if(sum(chsiaug[,2]<c0[i])==0){
				y2 = yaug
				XS2 = XEaug
				which2 = 1:length(yaug)
			}
		}else{
			f1 = chsiaug[,1]
			y1 =yaug
			XS1 = Xaug
			f2 = chsiaug[,2]
			y2 = yaug
			XS2 = XEaug
			which1=which2 = 1:length(yaug)
		}


		#CM1 for beta using scaled covariate matrix
		if(!is.null(de1) & !is.null(ta1)){	
			if(glmtype=="pois"){
				m1=glmIAL(y=y1, X=scale(XS1), prior=f1[which1],   prop=mean(f1), pMax=dim(XS1)[2], delta=de1, tau=ta1, nReEstimate=0,maxitIAL=maxitIAL, maxit=25, conv=glmconv, fitted=m1$fitted[which1],family="poisson", protect=dim(XS1)[2])
			}else{
				m1=glmNB.IAL(y=y1, X=scale(XS1), prior=f1[which1],   prop=mean(f1), pMax=dim(XS1)[2], delta=de1, tau=ta1,  nReEstimate=0, maxitIAL=maxitIAL, maxit=25, conv=glmconv, fitted=m1$fitted[which1], phi=phi1, protect=dim(XS1)[2])
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
				#if previous phi on boundary, get new initial estimate for phi so doesn't forever stay on boundary
				if(phi1 == 10^-5){
					phi1 = 1/theta.mm(dfr=sum(f1[which1])-p-2, mu=m1$fitted, weights=f1[which1], y=y1)
					if(phi1 < 10^-5) phi1 = .1
				}
				phi1 = .C("phi_mlR", as.double(y1), as.double(m1$fitted), as.integer(length(y1)), as.integer(20), as.double(10^-6), as.double(phi1), as.integer(1), as.integer(0), as.double(f1[which1]), package="hmmcov")[[6]]
			}

			if(c0[i]>0){
				if(length(m1$w2use)> 0){
					m1$fitted  = exp(as.numeric(Xaug[,m1$w2use]%*%matrix(m1$b2use, m1$n2use, 1))+m1$b0)
				}else{
					m1$fitted  = rep(exp(m1$b0), length(y))
				}
			}	
		}

		if(!is.null(de2) & !is.null(ta2)){
			if(glmtype=="pois"){
				m2=glmIAL(y=y2, X=scale(XS2), prior=f2[which2],   prop=mean(f2), pMax=dim(XS2)[2], delta=de2, tau=ta2, nReEstimate=0,maxitIAL=maxitIAL, maxit=25, conv=glmconv, fitted=m2$fitted[which2],family="poisson", protect=dim(XS2)[2])
			}else{
				m2=glmNB.IAL(y=y2, X=scale(XS2), prior=f2[which2],   prop=mean(f2), pMax=dim(XS2)[2], delta=de2, tau=ta2, nReEstimate=0,maxitIAL=maxitIAL, maxit=25, conv=glmconv, fitted=m2$fitted[which2], phi=phi2, protect=dim(XS2)[2])
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
				#if previous phi on boundary, get new initial estimate for phi so doesn't forever stay on boundary
				if(phi2 == 10^-5){
					phi2 = 1/theta.mm(dfr=sum(f2[which2])-p-2, mu=m2$fitted, weights=f2[which2], y=y2)
					if(phi2 < 10^-5) phi2 = .1
				}
				phi2 = .C("phi_mlR", as.double(y2), as.double(m2$fitted), as.integer(length(y2)), as.integer(20), as.double(10^-6), as.double(phi2), as.integer(1), as.integer(0), as.double(f2[which2]), package="hmmcov")[[6]]
			}

			if(c0[i]>0){
				if(length(m2$w2use)> 0){
					m2$fitted  = exp(as.numeric(XEaug[,m2$w2use]%*%matrix(m2$b2use, m2$n2use, 1))+m2$b0)
				}else{
					m2$fitted  = rep(exp(m2$b0), length(y))
				}
			}
		}
		
		#E step
		#calculate emission probabilities for forward-backward algorithm on log scale

		prob=matrix(0, n, K^2)
		if(glmtype!="nb"){
			prob[,1]=dpois(y, lambda=m1$fitted[seq(1, 2*n, by=2)], log=T)
			prob[,2]=dpois(y, lambda=m1$fitted[seq(2, 2*n, by=2)], log=T)
			prob[,3]=dpois(y, lambda=m2$fitted[seq(1, 2*n, by=2)], log=T)
			prob[,4]=dpois(y, lambda=m2$fitted[seq(2, 2*n, by=2)], log=T)
		}else{
			prob[,1]=dnbinom(y, mu=m1$fitted[seq(1, 2*n, by=2)], size=1/phi1, log=T)
			prob[,2]=dnbinom(y, mu=m1$fitted[seq(2, 2*n, by=2)], size=1/phi1, log=T)
			prob[,3]=dnbinom(y, mu=m2$fitted[seq(1, 2*n, by=2)], size=1/phi2, log=T)
			prob[,4]=dnbinom(y, mu=m2$fitted[seq(2, 2*n, by=2)], size=1/phi2, log=T)
		}

		if(any(pi == 0)) pi = abs(pi - 10^-100)
		resu=.C("forwardbackaug", as.double(log(Pi)), as.double(log(pi)), as.double(prob), as.integer(nrow(prob)), as.integer(K), logalpha = double(nrow(prob)*K),logbeta = double(nrow(prob)*K),LL = double(1), package="hmmcov")
		forward = matrix(resu$logalpha, ncol=2) 
		backward = matrix(resu$logbeta, ncol=2) 
		forwardbackward = exp(forward + backward - resu$LL)	
		ll[i] = resu$LL

		chsi2=matrix(0, n, K^2)			
		chsi2[-1,1]=forward[-n,1]+log(Pi[1,1])+prob[-1,1]+backward[-1,1]-resu$LL
		chsi2[-1,2]=forward[-n,2]+log(Pi[2,1])+prob[-1,2]+backward[-1,1]-resu$LL
		chsi2[-1,3]=forward[-n,1]+log(Pi[1,2])+prob[-1,3]+backward[-1,2]-resu$LL
		chsi2[-1,4]=forward[-n,2]+log(Pi[2,2])+prob[-1,4]+backward[-1,2]-resu$LL
		chsi2=exp(chsi2)

		chsiaug=matrix(0, n*K,K)	
		chsiaug[seq(1, 2*n, by=2),1] = chsi2[,1]
		chsiaug[seq(2, 2*n, by=2),1] = chsi2[,2]
		chsiaug[seq(1, 2*n, by=2),2] = chsi2[,3]
		chsiaug[seq(2, 2*n, by=2),2] = chsi2[,4]

		chsiaug[1,] = c(0,0)
		chsiaug[2,] = forwardbackward[1,]


		#end of E-step, check stopping criterion
		if(i>1) if(abs((ll[i] - ll[i-1])/ll[i-1])<EMconv) break
		#end of E-step		 

		#recalculate AR covariate
		
		res=res0=rep(0,2*n)
			if(length(m1$b2use)>1){
			res0[seq(1, 2*n, by=2)]=(log(y + 1) - log(exp(X[,m1$w2use[-length(m1$w2use)]]%*%matrix(m1$b2use[-length(m1$b2use)], length(m1$b2use)-1, 1)+m1$b0)+1))
		}else{
			res0[seq(1, 2*n, by=2)]=(log(y + 1) - log(exp(m1$b0)+1))
		}

		if(length(m2$b2use)>1){
			res0[seq(2, 2*n, by=2)]=(log(y + 1) - log(exp(XE[,m2$w2use[-length(m2$w2use)]]%*%matrix(m2$b2use[-length(m2$b2use)],length(m2$b2use)-1, 1)+m2$b0)+1))
		}else{
			res0[seq(2, 2*n, by=2)]=(log(y + 1) - log(exp(m2$b0)+1))
		}

		res=c(0,0,res0[-((2*n-1):(2*n))])
		Xaug[,p]=res
		XEaug[,p2]=res

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
	time=(ptm <- proc.time())

	m1final=rep(0,p)
	m1final[m1$w2use]=m1$b2use
	m2final=rep(0,p2)
	m2final[m2$w2use]=m2$b2use

	finalcoef=matrix(c(m1$b0, m1final, m2$b0,m2final, a, pi), 1, length(c(m1$b0, m1final, m2$b0,m2final, a, pi))) 
	if(glmtype=="nb") finalcoef=matrix(c(m1$b0, m1final, 1/phi1, m2$b0,m2final,1/phi2, a, pi), 1, length(c(m1$b0, m1final, m2$b0,m2final, a, pi, 1/phi1, 1/phi2))) 

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
time=time, unstcoef2=unstcoef2, ll=max(ll[ll<0]), monotone=0, monotone2=0, BIC=-2*ll[ll<0][sum(ll<0)]+sum(finalcoef!=0)*log(n), de1=de1, de2 = de2 , ta1=ta1, ta2=ta2, m1=0, m2=0)		
	}else{
		finalresults=list(y=0, X=0, forwardback=forwardbackward, forward=0, resid=0, chsi2=0, Pi=Pi, pi=pi, coef=finalcoef, unstcoef1=unstcoef1,probi2=1, 
time=time, unstcoef2=unstcoef2, ll=ll[ll<0][sum(ll<0)], monotone=0, monotone2=0, BIC=-2*ll[ll<0][sum(ll<0)]+sum(finalcoef!=0)*log(n), de1=de1, de2 = de2 , ta1=ta1, ta2=ta2, m1=m11, m2=m22)		
	}

	return(finalresults)
}

#wrapper for hmm and arhmm, not to be used with variable selection
#X must not have an intercept column
arhmm = function(y, X, prop1,maxitIAL=1000,thresh=0.05, XE=NULL,maxitEM=50, glmtype="pois", EMconv=10^-6, glmconv=10^-5){
	
	if(length(y)!=nrow(X)) stop("length of y does not match rows of X, or X is not a matrix")

	n = length(y)
	t = quantile(y,prop1)
	Pi = matrix(0, 2, 2)
	Pi[1,1] = sum(y[-1] <= t & y[-n] <= t)/sum(y[-1]<=t)
	Pi[1,2] = sum(y[-1] <= t & y[-n] >  t)/sum(y[-1]<=t)
	Pi[2,1] = sum(y[-1] > t & y[-n] <= t)/sum(y[-1]>t)
	Pi[2,2] = sum(y[-1] > t & y[-n] > t)/sum(y[-1]>t)
	
	pi = abs(c(y[1]<=t, y[1]>t) - 10^-100)

	res = arhmmcov(de1=0, ta1=-1, de2=0, ta2=-1, n=n, y=y, X=X, prop1=prop1, pi=pi, Pi=Pi,m10=NULL, m20=NULL, 
		glmtype=glmtype, maxitIAL=maxitIAL, thresh=thresh, XE = XE, maxitEM=maxitEM, EMconv=EMconv, glmconv=glmconv)

	coef  = matrix(res$coef[1:(length(res$coef)-6)], ncol=2)
	rownames(coef) = rep("", nrow(coef))

	if(glmtype == "nb"){
		rownames(coef)[nrow(coef)] = "phi"
		rownames(coef)[nrow(coef)-1] = "nu"
	}else if(glmtype=="pois"){
		rownames(coef)[nrow(coef)] = "nu"
	}
	final = list(forwardbackward = res$forwardback, coef = coef, ll=res$ll, Pi = res$Pi, pi=res$pi)	
	return(final)
}



arhmmvarsel = function(y, X, prop1,maxitIAL=1000,thresh=0.05, XE=NULL,maxitEM=50, glmtype="pois", EMconv=10^-6, glmconv=10^-5, pen=NULL, de1=NULL, ta1=NULL, de2=NULL, ta2=NULL, m1 = NULL, m2=NULL, Pi=NULL, pi=NULL){
	
	if(length(y)!=nrow(X)) stop("length of y does not match rows of X, or X is not a matrix")


	n = length(y)

	if(is.null(de1) & is.null(ta1) & is.null(de2) & is.null(ta2) & is.null(m1) & is.null(m2)){
		#initial full fit 
		#do all initializations
		t = quantile(y,prop1)
		Pi = matrix(0, 2, 2)
		Pi[1,1] = sum(y[-1] <= t & y[-n] <= t)/sum(y[-1]<=t)
		Pi[1,2] = sum(y[-1] <= t & y[-n] >  t)/sum(y[-1]<=t)
		Pi[2,1] = sum(y[-1] > t & y[-n] <= t)/sum(y[-1]>t)
		Pi[2,2] = sum(y[-1] > t & y[-n] > t)/sum(y[-1]>t)
	
		pi = abs(c(y[1]<=t, y[1]>t) - 10^-100)

		res = arhmmcov(de1=0, ta1=-1, de2=0, ta2=-1, n=n, y=y, X=X, prop1=prop1, pi=pi, Pi=Pi,m10=NULL, m20=NULL, 
			glmtype=glmtype, maxitIAL=maxitIAL, thresh=thresh, XE = XE, maxitEM=maxitEM, EMconv=EMconv, glmconv=glmconv)
			final = list(forwardbackward= res$forwardback, m1=res$m1, m2 = res$m2, Pi = res$Pi, pi=res$pi, BIC=res$BIC)	
	}else{
		if(is.null(Pi)) stop("Pi needs to be specified")
		if(is.null(pi)) stop("pi needs to be specified")
		if(pen == "scad"){
			ta1 = ta2 = -2
			if(is.null(de1) & is.null(de2)) stop("de1 or de2 at least needs to be specified")
		}else if(pen == "lasso"){
			ta1 = ta2 = -1
			if(is.null(de1) & is.null(de2)) stop("de1 or de2 at least needs to be specified")
		}else if(pen == "log"){
			if(is.null(de1) & is.null(de2) & is.null(ta1) & is.null(ta2)) stop("de1, de2, ta1,and ta2 all need to be specified")
		}
		
		res = arhmmcov(de1=de1, ta1=ta1, de2=de2, ta2=ta2, n=n, y=y, X=X, prop1=prop1, pi=pi, Pi=Pi,m10=m1, m20=m2, 
			glmtype=glmtype, maxitIAL=maxitIAL, thresh=thresh, XE = XE, maxitEM=maxitEM, EMconv=EMconv, glmconv=glmconv)
			final = list(forwardbackward= res$forwardback, Pi = res$Pi, pi=res$pi, BIC=res$BIC, de1=res$de1, de2=res$de2, ta1=res$ta1, ta2=res$ta2)	
	}

	coef  = matrix(res$coef[1:(length(res$coef)-6)], ncol=2)
	rownames(coef) = rep("", nrow(coef))

	if(glmtype == "nb"){
		rownames(coef)[nrow(coef)] = "phi"
		rownames(coef)[nrow(coef)-1] = "nu"
	}else if(glmtype=="pois"){
		rownames(coef)[nrow(coef)] = "nu"
	}
	final$coef = coef	
	return(final)
}

