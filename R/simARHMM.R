simARHMM = function(type, n=10000, beta1 = NULL, beta2 = NULL, nu=c(0,0), glmtype = "pois",phi=4){
  
  library(HiddenMarkov)
  
	K=2
	if(!is.null(beta1) & !is.null(beta2)){
			if(length(beta1)!=length(beta2)) stop("for this simulation, number of covariates need to be the same in each state")
	}

	#histone parameters
	if(type=="histone"){
		prop1=.9
    pi=c(.9, 1-.9)
		a=c(prop1, 1-prop1, 1-prop1, prop1)
		if(is.null(beta1)) beta1 = c(0, 1, 1)
		if(is.null(beta2)) beta2 = c(.5, 2, 2)
		p=length(beta1)-1
	}else if(type=="ctcf"){
		prop1=.9
    pi=c(.9, 1-.9)
		a=c(prop1, 1-prop1, prop1, 1-prop1)
		if(is.null(beta1)) beta1 = c(0, 1, 1)
		if(is.null(beta2)) beta2 = c(1.5, 2, 2)
		p=length(beta1)-1
	}

	Pi <- matrix(a,byrow=TRUE, nrow=K)
	#simulate arbitrary hmm object using Pi and pi, only using to pull simulated states from object 
	x <- dthmm(NULL, Pi, pi, "pois",
	           list(lambda=c(4, 10)), discrete = TRUE)
	x <- simulate(x, nsim=n)
	#pull states
	states=x$y
	#make n x p matrix	
	X=matrix(rbeta(n*p, 1, 1), n, p)
	X1=cbind(rep(1, n), X)
	beta=list(beta1, beta2)

	#begin recursive simulation
	y=rep(0, n)
	if(glmtype=="pois"){		
		y[1]=rpois(1, lambda=exp((X1%*%beta[[states[1]]])[1]))
		for(i in 2:n){
			Ri=log(y[i-1]+1 ) - log(exp((X1%*%beta[[states[i-1]]])[i-1])+1)
			y[i]=rpois(1, lambda=exp((X1%*%beta[[states[i]]])[i]+ nu[states[i]]*(Ri)))
		}
	}else{
		#same phi in each component
		y[1]=rnbinom(1, mu=exp((X1%*%beta[[states[1]]])[1]), size=phi)
		for(i in 2:n){
			Ri=log(y[i-1]+1) - log(exp((X1%*%beta[[states[i-1]]])[i-1])+1)
			y[i]=rnbinom(1, mu=exp((X1%*%beta[[states[i]]])[i]+ nu[states[i]]*(Ri)), size=phi)
		}
	}
	return(list(X=X, y=y, states=states))
}


