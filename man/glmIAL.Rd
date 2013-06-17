\name{glmIAL}
\alias{glmIAL}
\title{
Fit a generalized linear model using interative adaptive Lasso penalty
}
\description{
Fit a generalized linear model, with varaible selection by
penalized likelihood method. Now only support faimliy 
"binomial", "poisson", and "gamma".
}
\usage{
glmIAL(family, y, X, delta=seq(0, 5.0, by=0.5),  
       tau = c(0.01*c(10,8,6,4,2,1), 0.0001*c(50,20,10,5,2,1)),
       pMax=20, offset=NULL, naPercent=0.4, 
       nTotal=NULL, maxit=50, maxitIAL=20, nReEstimate=pMax, 
       conv=1e-5, fitted=NULL, trace=1)
}
\arguments{
  \item{family}{
	a character string specifying the distribution family, the allowed values 
	inlcude "binomial", "poisson", and "gamma". Canonical link 
	function are assumed for them, which is "logit", "log", and "inverse" for 
	"binomial", "poisson", and "gamma"}
	
  \item{y}{
	a vector of the response variable. }
	
  \item{X}{
	a matrix of covaraite data, each column is one covariate and 
	each row is one observation. Note, do not include intercept in 
	this covariate matrix. }
	
  \item{delta}{
	one of the two penalization parameters for iterative adpative Lasso }
	
  \item{tau}{
	one of the two penalization parameters for iterative adpative Lasso }

  \item{pMax}{
	the maiximum number of covariates to be kept in the model }
	
  \item{offset}{
	a vector specifying offset in the linear predictor scale, eta = Xbeta + offset.}
	
  \item{naPercent}{
	the maximum percent of missing that are allowed in the data. }
	
  \item{nTotal}{
	if faimly is "binomail", nTotal specify the total number of trails 
	in each subject}
	
  \item{maxit}{
	the maximum number of iteration of Iterated Reweighted Least Square }
	 
  \item{maxitIAL}{
	the maximum number of iteration of for IAL }

  \item{nReEstimate}{
	If there is less than nReEstimate non-zero coefficients, re-estimate these 
  coefficients by least squares in IAL function. }

  \item{conv}{
	convergence tolerence. Let wss be the weighted sum of perason residual square, 
	the algorithm converges if |wss-wss_last|/wss_last < conv }
	
  \item{fitted}{
	a vector of the fitted value, should be the same length as y}
	
  \item{trace}{ Print out tracking information if trace > 1.}
}

\value{
A list. 
}
\author{
Wei Sun weisun@email.unc.edu
}
\note{
Part of codes is adopted from CNVtools/glm_test.cpp
}

\keyword{ methods }
