#
#  repeated : A Library of Repeated Measurements Models
#  Copyright (C) 2000, 2001 J.K. Lindsey
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public Licence as published by
#  the Free Software Foundation; either version 2 of the Licence, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public Licence for more details.
#
#  You should have received a copy of the GNU General Public Licence
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#  SYNOPSIS
#
#     gausscop(response=NULL, distribution="gamma", mu=NULL, shape=NULL,
#	autocorr="exponential", pmu=NULL, pshape=NULL,
#	par=NULL, pre=NULL, delta=NULL, shfn=FALSE, common=FALSE,
#	envir=parent.frame(), print.level=0, ndigit=10,
#	gradtol=0.00001, steptol=0.00001, iterlim=100, fscale=1,
#	stepmax=10*sqrt(theta%*%theta), typsize=abs(c(theta)))
#
#  DESCRIPTION
#
#    Function to fit the multivariate Gaussian copulas with various
# marginal distributions, various autocorrelation functions, one or
# two levels of random effects, and nonlinear regression.

gausscop <- function(response=NULL, distribution="gamma", mu=NULL, shape=NULL,
	autocorr="exponential", pmu=NULL, pshape=NULL, par=NULL,
	pre=NULL, delta=NULL, shfn=FALSE, common=FALSE, envir=parent.frame(),
	print.level=0, ndigit=10, gradtol=0.00001, steptol=0.00001,
	iterlim=100, fscale=1, stepmax=10*sqrt(theta%*%theta),
	typsize=abs(c(theta))){
#
# distribution functions and densities
#
pinvgauss <- function(q, m, s){
	t <- q/m
	v <- sqrt(q*s)
	pnorm((t-1)/v)+exp(2/(m*s))*pnorm(-(t+1)/v)}
dinvgauss <- function(y, m, s)
	exp(-(y-m)^2/(2*y*s*m^2))/sqrt(2*pi*s*y^3)
plaplace <- function(q, m=0, s=1){
	u <- (q-m)/s
	t <- exp(-abs(u))/2
	ifelse(u<0,t,1-t)}
dlaplace <- function(y, m=0, s=1)
	exp(-abs(y-m)/s)/(2*s)
plevy <- function(q, m=0, s=1){
	if(any(q<m))stop("some y <= m")
	len <- length(q)
	z <- .C("plevy",
		as.double(q),
		as.double(m),
		as.double(s),
		as.double(1),
		len=as.integer(len),
		eps=as.double(1.0e-6),
		pts=as.integer(5),
		max=as.integer(16),
		err=integer(1),
		res=double(len),
		DUP=FALSE,
		PACKAGE="rmutil")
	z$res}
dlevy <- function(y, m=0, s=1){
	if(any(y<=m))stop("some y <= m")
	sqrt(s/(2*pi*(y-m)^3))*exp(-s/(2*(y-m)))}
ppareto <- function(q, m, s){
	if(any(s<=1))stop("s must be > 1")
	1-(1+q/(m*(s-1)))^-s}
dpareto <- function(y, m, s){
	if(any(s<=1))stop("s must be > 1")
	m <- m*(s-1)
	s*(1+y/m)^(-s-1)/m}
#
# likelihood function
#
gcopula <- function(theta){
	mu <- mu1(theta)
	if(dst>1)shape <- exp(shape1(theta))
	if(!indep){
		corr <- 1/(1+exp(-theta[npd1:np]))
		switch(dst,
			"1"=u <- pexp(y,exp(mu)),
			"2"=u <- pgamma(y,shape,scale=exp(mu)/shape),
			"3"=u <- pweibull(y,shape,exp(mu)),
			"4"=u <- ppareto(y,exp(mu),shape),
			"5"=u <- pinvgauss(y,exp(mu),shape),
			"6"=u <- plogis(y,mu,shape),
			"7"=u <- pcauchy(y,mu,shape),
			"8"=u <- plaplace(y,mu,shape),
			"9"=u <- plevy(y,mu,shape))
		if(any(u==1,na.rm=TRUE))u[u==1] <- 1-1e-15
		if(any(u==0,na.rm=TRUE))u[u==0] <- 1e-15
		q <- qnorm(u)
		z <- .Fortran("gcopula",
			theta=as.double(corr),
			like=double(1),
			x=as.double(response$response$times),
			y=as.double(q),
			nobs=as.integer(nobs),
			nest=as.integer(response$response$nest),
			lnest=as.integer(lnest),
			nind=as.integer(nind),
			nld=as.integer(nld),
			npre=as.integer(npre),
			npar=as.integer(npar),
			ar=as.integer(ar),
			v=double(nld*nld),
			tmp1=double(nld),
			tmp2=double(nld),
			tmp3=double(nld),
			warn=logical(1),
			DUP=FALSE,
			PACKAGE="repeated")
		if(z$warn)stop("Correlation matrix for some individuals not positive definite")}
	switch(dst,
		"1"=u <- sum(dexp(y,exp(mu),TRUE)),
		"2"=u <- sum(dgamma(y,shape,scale=exp(mu)/shape,log=TRUE)),
		"3"=u <- sum(dweibull(y,shape,exp(mu),TRUE)),
		"4"=u <- sum(log(dpareto(y,exp(mu),shape))),
		"5"=u <- sum(log(dinvgauss(y,exp(mu),shape))),
		"6"=u <- sum(dlogis(y,mu,shape,TRUE)),
		"7"=u <- sum(dcauchy(y,mu,shape,TRUE)),
		"8"=u <- sum(log(dlaplace(y,mu,shape))),
		"9"=u <- sum(log(dlevy(y,mu,shape))))
	if(indep) -u else z$like-u}
call <- sys.call()
#
# check type of model to be fitted
#
tmp <- c("exponential","gamma","Weibull","Pareto",
	"inverse Gauss","logistic","Cauchy","Laplace","Levy")
dst <- match(distribution <- match.arg(distribution,tmp),tmp)
tmp <- c("exponential","gaussian","cauchy","spherical")
ar <- match(autocorr <- match.arg(autocorr,tmp),tmp)
npr <- length(pmu)
npr1 <- npr+1
npshape <- length(pshape)
npd1 <- npr+npshape+1
if(dst==1){
	npshape <- 0
	pshape <- NULL}
npar <- length(par)
npre <- length(pre)
indep <- (npar+npre)==0
#
# check if location and shape can have common parameters
#
if(common){
	if(!is.function(mu)&&!inherits(mu,"formula"))
		stop("with common parameters, mu must be a function or formula")
	if(!is.function(shape)&&!inherits(shape,"formula"))
		stop("with common parameters, shape must be a function or formula")
	pshape <- NULL}
#
# check if a data object is being supplied
#
type <- "unknown"
respenv <- exists(deparse(substitute(response)),envir=parent.frame())&&
	inherits(response,"repeated")&&!inherits(envir,"repeated")
if(respenv){
	if(dim(response$response$y)[2]>1)
		stop("gausscop only handles univariate responses")
	if(!is.null(response$NAs)&&any(response$NAs))
		stop("gausscop does not handle data with NAs")
	type <- response$response$type}
envname <- if(respenv)paste(deparse(substitute(response)))
	else if(!is.null(class(envir)))deparse(substitute(envir))
	else NULL
#
# if envir, remove extra (multivariate) responses
#
if(!respenv&&inherits(envir,"repeated")){
	if(!is.null(envir$NAs)&&any(envir$NAs))
		stop("gausscop does not handle data with NAs")
	cn <- deparse(substitute(response))
	if(length(grep("\"",cn))>0)cn <- response
	if(length(cn)>1)stop("only one response variable allowed")
	response <- envir
	col <- match(cn,colnames(response$response$y))
	if(is.na(col))stop(paste("response variable",cn,"not found"))
	type <- response$response$type[col]
	if(dim(response$response$y)[2]>1){
		response$response$y <- response$response$y[,col,drop=FALSE]
		if(!is.null(response$response$delta)){
			response$response$delta <- response$response$delta[,col,drop=FALSE]
			if(all(response$response$delta==1)||all(is.na(response$response$delta)))response$response$delta <- NULL}}}
if(inherits(response,"response")){
	type <- response$type
	response <- rmna(response)}
if((inherits(envir,"repeated")&&(length(nobs(response))!=length(nobs(envir))||
	any(nobs(response)!=nobs(envir))))||(inherits(envir,"tvcov")&&
	(length(nobs(response))!=length(envir$tvcov$nobs)||
	any(nobs(response)!=envir$tvcov$nobs))))
	stop("response and envir objects are incompatible")
if((distribution=="exponential"||distribution=="gamma"||distribution=="Weibull"
	||distribution=="Pareto"||distribution=="inverse Gauss")
	&&type!="unknown"&&type!="duration"&&type!="continuous")
	stop("duration data required")
if((distribution=="logistic"||distribution=="Cauchy"||distribution=="Laplace"
	||distribution=="Levy")&&type!="unknown"&&type!="continuous"
	&&type!="duration")stop("continuous data required")
y <- response$response$y
if(dst<6&&any(y<=0))
	stop(paste("all responses must be positive for",distribution,"distribution"))
n <- length(y)
nobs <- nobs(response)
nind <- length(nobs)
nld <- max(nobs)
#
# if a data object was supplied, modify formulae or functions to read from it
#
mu3 <- sh3 <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")||inherits(envir,"tvcov")){
	if(is.null(envname))envname <- deparse(substitute(envir))
	if(inherits(mu,"formula")){
		mu3 <- if(respenv)finterp(mu,.envir=response,.name=envname)
			else finterp(mu,.envir=envir,.name=envname)}
	else if(is.function(mu)){
		if(is.null(attr(mu,"model"))){
		        tmp <- parse(text=deparse(mu)[-1])
		        mu <- if(respenv)fnenvir(mu,.envir=response,.name=envname)
		        	else fnenvir(mu,.envir=envir,.name=envname)
		        mu3 <- mu
		        attr(mu3,"model") <- tmp}
		else mu3 <- mu}
	if(inherits(shape,"formula")){
		sh3 <- if(respenv)finterp(shape,.envir=response,.name=envname)
			else finterp(shape,.envir=envir,.name=envname)}
	else if(is.function(shape)){
		if(is.null(attr(shape,"model"))){
		        tmp <- parse(text=deparse(shape)[-1])
		        shape <- if(respenv)fnenvir(shape,.envir=response,.name=envname)
		        	else fnenvir(shape,.envir=envir,.name=envname)
		        sh3 <- shape
		        attr(sh3,"model") <- tmp}
		else sh3 <- shape}}
#
# transform mu formula to function and check number of parameters
#
if(inherits(mu,"formula")){
	mu2 <- if(respenv)finterp(mu,.envir=response,.name=envname)
		else finterp(mu,.envir=envir,.name=envname)
	npt1 <- length(attr(mu2,"parameters"))
	if(is.character(attr(mu2,"model"))){
	# W&R formula
		if(length(attr(mu2,"model"))==1){
			mu1 <- function(p) p[1]*rep(1,n)
			attributes(mu1) <- attributes(mu2)
			mu2 <- NULL}}
	else {
	# formula with unknowns
		if(npr!=npt1&&!common){
			cat("\nParameters are ")
			cat(attr(mu2,"parameters"),"\n")
			stop(paste("pmu should have",npt1,"estimates"))}
		if(is.list(pmu)){
			if(!is.null(names(pmu))){
				o <- match(attr(mu2,"parameters"),names(pmu))
				pmu <- unlist(pmu)[o]
				if(sum(!is.na(o))!=length(pmu))stop("invalid estimates for mu - probably wrong names")}
			else pmu <- unlist(pmu)}}
	if(!is.null(mu2)){
		if(inherits(envir,"tccov")){
		# fix length if only time-constant covariates
			cv <- covind(response)
			mu1 <- function(p) mu2(p)[cv]
			attributes(mu1) <- attributes(mu2)}
		else {
			mu1 <- mu2
			rm(mu2)}}}
else if(is.function(mu))mu1 <- mu
else {
	if(npr!=1)stop("pmu should have one estimate")
	mu1 <- function(p) rep(p[1],n)}
#
# give appropriate attributes to mu1 for printing
#
if(!is.null(mu1)&&is.null(attr(mu1,"parameters"))){
	attributes(mu1) <- if(is.function(mu)){
		if(!inherits(mu,"formulafn")){
			if(respenv)attributes(fnenvir(mu,.envir=response))
			else attributes(fnenvir(mu,.envir=envir))}
		else attributes(mu)}
		else {
			if(respenv)attributes(fnenvir(mu1,.envir=response))
			else attributes(fnenvir(mu1,.envir=envir))}}
#
# if possible, check that correct number of estimates was supplied
#
nlp <- if(is.function(mu1))length(attr(mu1,"parameters"))
	else if(is.null(mu1))NULL
	else npt1
if(!is.null(nlp)&&!common&&nlp!=npr)
	stop(paste("pmu should have",nlp,"initial estimates"))
if(length(mu1(pmu))!=sum(nobs(response)))
	stop("The location function must provide an estimate for each observation")
n1 <- if(common&&!inherits(shape,"formula")) 1 else npr1
n2 <- npr+npshape
#
# transform shape formula to function and check number of parameters
#
if(inherits(shape,"formula")){
	old <- if(common)mu1 else NULL
	mufn <- if(shfn)"mu" else NULL
	shape4 <- if(respenv)finterp(shape,.envir=response,.start=n1,.name=envname,.old=old,.args=mufn)
		else finterp(shape,.envir=envir,.start=n1,.name=envname,.old=old,.args=mufn)
	tmp <- attributes(shape4)
	shape2 <- if(shfn)function(p) shape4(p,mu1(p)) else shape4
	attributes(shape2) <- tmp
	npt2 <- length(attr(shape2,"parameters"))
	if(is.character(attr(shape2,"model"))){
	# W&R formula
		if(length(attr(shape2,"model"))==1){
			shape1 <- function(p) p[n1]*rep(1,n)
			attributes(shape1) <- attributes(shape2)
			shape2 <- NULL}}
	else {
	# formula with unknowns
		if(npshape!=npt2&&!common){
			cat("\nParameters are ")
			cat(attr(shape2,"parameters"),"\n")
			stop(paste("pshape should have",npt2,"estimates"))}
		if(is.list(pshape)){
			if(!is.null(names(pshape))){
				o <- match(attr(shape2,"parameters"),names(pshape))
				pshape <- unlist(pshape)[o]
				if(sum(!is.na(o))!=length(pshape))stop("invalid estimates for shape - probably wrong names")}
			else pshape <- unlist(pshape)}}
	if(!is.null(shape2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			shape1 <- function(p) shape2(p)[cv]
			attributes(shape1) <- attributes(shape2)}
		else {
			shape1 <- shape2
			rm(shape2)}}}
else if(is.function(shape))shape1 <- if(shfn) function(p) shape(p[n1:n2],mu1(p))
		else function(p) shape(p[n1:n2])
else if(dst>1){
	if(npshape!=1)stop("pshape should have one estimate")
	shape1 <- function(p) rep(p[npr1],n)}
else shape1 <- NULL
#
# give appropriate attributes to shape1 for printing
#
if(dst>1){
	if(!is.null(shape1)&&is.null(attr(shape1,"parameters")))
		attributes(shape1) <- if(is.function(shape)){
			if(!inherits(shape,"formulafn")){
				if(respenv)attributes(fnenvir(shape,.envir=response))
				else attributes(fnenvir(shape,.envir=envir))}
			else attributes(shape)}
	else {
		if(respenv)attributes(fnenvir(shape1,.envir=response))
		else attributes(fnenvir(shape1,.envir=envir))}
	nlp <- if(is.function(shape))length(attr(shape1,"parameters"))-shfn
		else if(is.null(shape)||is.character(shape))NULL
		else if(!inherits(shape,"formula"))stop("shape must be a function or formula")
		else npt2}
#
# set up nesting indicator
#
if(!is.null(response$response$nest))lnest <- max(response$response$nest)
else {
	lnest <- 0
	response$response$nest <- rep(1,length(y))}
#
# check times
#
if(is.null(response$response$times)&&!is.null(par))
	stop("No times. AR cannot be fitted")
nm <- sum(nobs(response))
#
# include delta in Jacobian
#
jacob <- 0
if(!is.null(response$response$delta)){
	if(length(response$response$delta)==1)jacob <- jacob-nm*log(response$response$delta)
	else jacob <- jacob -sum(log(response$response$delta))}
#
# set up shape function if there is one and check initial estimates
#
theta <- c(pmu,pshape)
#
# check AR and random effect
#
np <- npr+npshape+npre+npar
if(!indep&&sum(c(pre,par))>=1)
	stop("Sum of correlations (pre+par) must be less than 1")
if(npre>0){
	if(any(pre<=0)||any(pre>=1))
		stop("All variance components must lie between 0 and 1")
	theta <- c(theta,log(pre/(1-pre)))}
if(npar>0){
	if(par<=0||(par>=1))
		stop("Estimate of autocorrelation must lie between 0 and 1")
	theta <- c(theta,log(par/(1-par)))}
#
# estimate model
#
if(fscale==1)fscale <- gcopula(theta)
z0 <- nlm(gcopula, p=theta, hessian=TRUE, print.level=print.level,
	typsize=typsize, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
p <- z0$estimate
z <- gcopula(p)
like <- z+jacob
pred <- if(dst<6)exp(mu1(p)) else mu1(p)
#
# get parameters
#
if(npre>0){
	tausq <- exp(p[npd1:(npr+npshape+npre)])
	tausq <- tausq/(1+tausq)}
else tausq <- 0
if(npar>0){
	rho <- exp(p[npd1+npre])
	rho <- rho/(1+rho)}
else rho <- 0
#
# calculate se's
#
if(np==1){
	nlcov <- 1/z0$hessian
	nlse <- sqrt(nlcov)}
else {
	a <- if(any(is.na(z0$hessian))||any(abs(z0$hessian)==Inf))0
		else qr(z0$hessian)$rank
	if(a==np)nlcov <- solve(z0$hessian)
	else nlcov <- matrix(NA,ncol=np,nrow=np)
	nlse <- sqrt(diag(nlcov))}
if(length(nlse)>1)nlcorr <- nlcov/(nlse%o%nlse)
else nlcorr <- as.matrix(nlcov/nlse^2)
dimnames(nlcorr) <- list(seq(1,np),seq(1,np))
#
# return list
#
if(!is.null(mu3))mu1 <- mu3
if(!is.null(sh3))shape1 <- sh3
z <- list(
	call=call,
	distribution=distribution,
	mu1=mu1,
	shape1=shape1,
	common=common,
	shfn=shfn,
	autocorr=autocorr,
	response=response$response,
	maxlike=like,
	aic=like+np,
	df=nm-np,
	np=np,
	npr=npr,
	npshape=npshape,
	npar=npar,
	npre=npre,
	coefficients=p,
	nlse=nlse,
	nlcov=nlcov,
	nlcorr=nlcorr,
	tausq=tausq,
	rho=rho,
#	residuals=z$res,
	pred=pred,
	grad=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z) <- "gausscop"
return(z)}

### standard methods
###

deviance.gausscop <- function(z) 2*z$maxlike

fitted.gausscop <- function(z) z$pred

residuals.gausscop <- function(z) z$response$y-z$pred

### print method
###
print.gausscop <- function(z,digits=max(3,.Options$digits-3),correlation=TRUE){
cat("\nMultivariate normal copula with",z$distribution,"marginals\n")
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
cat("Number of subjects    ",length(nobs(z)),"\n")
cat("Number of observations",length(z$response$y),"\n")
cat("-Log likelihood   ",z$maxlike,"\n")
cat("Degrees of freedom",z$df,"\n")
cat("AIC               ",z$aic,"\n")
cat("Iterations        ",z$iterations,"\n\n")
if(z$common)cat("Location function\n")
else cat("Location function parameters\n")
if(!is.null(attr(z$mu1,"formula")))
	cat(deparse(attr(z$mu1,"formula")),sep="\n")
else if(!is.null(attr(z$mu1,"model"))){
	t <- deparse(attr(z$mu1,"model"))
	t[1] <- sub("expression\\(","",t[1])
	t[length(t)] <- sub("\\)$","",t[length(t)])
	cat(t,sep="\n")}
cname <- if(is.character(attr(z$mu1,"model")))attr(z$mu1,"model")
	else attr(z$mu1,"parameters")
coef.table <- cbind(z$coef[1:z$npr],z$nlse[1:z$npr])
if(!z$common){
	dimnames(coef.table) <- list(cname,c("estimate","se"))
	print.default(coef.table,digits=digits,print.gap=2)}
if(inherits(z$shape1,"formulafn")){
	if(z$common)cat("\nShape function\n")
	else {
		cat("\nShape function parameters\n")
		cname <- NULL}
	if(!is.null(attr(z$shape1,"formula")))
		cat(deparse(attr(z$shape1,"formula")),sep="\n")
	else if(!is.null(attr(z$shape1,"model"))){
		t <- deparse(attr(z$shape1,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	cname <- c(cname,if(is.character(attr(z$shape1,"model")))
		attr(z$shape1,"model")
		else attr(z$shape1,"parameters")[1:z$npshape])
	if(!z$common)coef.table <- cbind(z$coef[(z$npr+1):(z$npr+z$npshape)], z$nlse[(z$npr+1):(z$npr+z$npshape)])
	else cat("\nCommon parameters\n")
	dimnames(coef.table) <- list(unique(cname),c("estimate","se"))
	print.default(coef.table,digits=digits,print.gap=2)}
else if(!is.null(z$shape1)) {
	if(z$npshape==1)cat("\nShape\n")
	else cat("\nShape parameters\n")
	if(z$npshape>1)
		vname <- c("(Intercept)",paste("t^",1:(z$npshape-1),sep=""))
	else vname <- ""
	coef.table <- cbind(z$coef[(z$npr+1):(z$npr+z$npshape)],
		z$nlse[(z$npr+1):(z$npr+z$npshape)])
	dimnames(coef.table) <- list(vname,c("estimate","se"))
	print.default(coef.table,digits=digits,print.gap=2)}
if(z$npre>0){
	cat("\nVariance components\n")
	coef.table <- cbind(z$coef[(z$npr+z$npshape+1):(z$npr+z$npshape+z$npre)],
		z$nlse[(z$npr+z$npshape+1):(z$npr+z$npshape+z$npre)],z$tausq)
	if(z$npre==1)cname <- "tausq"
	else cname <- c("Level 1","Level 2")
	dimnames(coef.table) <- list(cname,c("estimate","se",""))
	print.default(coef.table,digits=digits,print.gap=2)}
if(z$rho!=0){
	cat("\n",z$autocorr," autocorrelation\n",sep="")
	coef.table <- cbind(z$coef[z$npr+z$npshape+z$npre+1],
		z$nlse[z$npr+z$npshape+z$npre+1],z$rho)
	dimnames(coef.table) <- list("rho",c("estimate","se",""))
	print.default(coef.table,digits=digits,print.gap=2)}
if(z$np>1&&correlation){
	cat("\nCorrelation matrix\n")
	print.default(z$nlcorr,digits=digits)}}
