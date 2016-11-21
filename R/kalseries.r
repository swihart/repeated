#
#  repeated : A Library of Repeated Measurements Models
#  Copyright (C) 1998, 1999, 2000, 2001 J.K. Lindsey
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
#     kalseries(response=NULL, times=NULL, intensity="exponential",
#	depend="independence", mu=NULL, shape=NULL, density=FALSE, ccov=NULL,
#	tvcov=NULL, torder=0, interaction=NULL, preg=NULL, ptvc=NULL,
#	pintercept=NULL, pshape=NULL, pinitial=1, pdepend=NULL,
#	pfamily=NULL, delta=NULL, transform="identity", link="identity",
#	envir=parent.frame(), print.level=0, ndigit=10, gradtol=0.00001,
#	steptol=0.00001, fscale=1, iterlim=100, typsize=abs(p),
#	stepmax=10*sqrt(p%*%p))
#
#  DESCRIPTION
#
#    A function to fit various distributions inserted into a Pareto
#  distribution with serial dependence or gamma frailties using
#  Kalman-type update for continuous longitudinal data.

.First.lib <- function(lib, pkg)
	library.dynam("repeated", pkg, lib)
require(rmutil)

kalseries <- function(response=NULL, times=NULL, intensity="exponential",
	depend="independence", mu=NULL, shape=NULL, density=FALSE, ccov=NULL,
	tvcov=NULL, torder=0, interaction=NULL, preg=NULL, ptvc=NULL,
	pintercept=NULL, pshape=NULL, pinitial=1, pdepend=NULL,
	pfamily=NULL, delta=NULL, transform="identity", link="identity",
	envir=parent.frame(), print.level=0, ndigit=10, gradtol=0.00001,
	steptol=0.00001, fscale=1, iterlim=100, typsize=abs(p),
	stepmax=10*sqrt(p%*%p)){
#
# Pareto (beta) likelihood function for serial dependence
#
series <- function(p){
	if(rf)b <- mu1(p)
	if(sf)v <- sh1(p[nps1:np])
	z <- .C("kserie",
		p=as.double(p),
		y=as.double(y),
		t=as.double(times),
		x=as.double(resp$ccov$ccov),
		nind=as.integer(nind),
		nobs=as.integer(nobs),
		nbs=as.integer(n),
		nccov=as.integer(nccov),
		npv=as.integer(npv),
		model=as.integer(mdl),
		link=as.integer(lnk),
		density=as.integer(density),
		pfamily=as.integer(!is.null(pfamily)),
		dep=as.integer(dep),
		torder=as.integer(torder),
		inter=as.integer(interaction),
		tvc=as.integer(tvc),
		tvcov=resp$tvcov$tvcov,
		fit=as.integer(0),
		pred=double(n),
		rpred=double(n),
		rf=as.integer(rf),
		bbb=as.double(b),
		sf=as.integer(sf),
		vv=as.double(v),
		like=double(1),
		DUP=FALSE,
		PACKAGE="repeated")
	z$like}
#
# Pareto (beta) likelihood function for frailty
#
serief <- function(p){
	if(rf)b <- mu1(p)
	if(sf)v <- sh1(p[nps1:np])
	z <- .C("krand",
		p=as.double(p),
		y=as.double(y),
		t=as.double(times),
		x=as.double(resp$ccov$ccov),
		nind=as.integer(nind),
		nobs=as.integer(nobs),
		nbs=as.integer(n),
		nccov=as.integer(nccov),
		npv=as.integer(npv),
		model=as.integer(mdl),
		link=as.integer(lnk),
		density=as.integer(density),
		torder=as.integer(torder),
		inter=as.integer(interaction),
		tvc=as.integer(tvc),
		tvcov=resp$tvcov$tvcov,
		fit=as.integer(0),
		pred=double(n),
		rpred=double(n),
		rf=as.integer(rf),
		bbb=as.double(b),
		sf=as.integer(sf),
		vv=as.double(v),
		frser=as.integer(frser),
		like=double(1),
		DUP=FALSE,
		PACKAGE="repeated")
	z$like}
call <- sys.call()
#
# check model
#
tmp <- c("exponential","Weibull","gamma","gen logistic","normal",
	"logistic","Cauchy","Laplace","log normal","log logistic",
	"log Cauchy","log Laplace","inverse Gauss")
mdl <- match(intensity <- match.arg(intensity,tmp),tmp)
tmp <- c("independence","Markov","serial","frailty")
dep <- match(depend <- match.arg(depend,tmp),tmp)-1
transform <- match.arg(transform,c("identity","exp","square","sqrt","log"))
tmp <- c("identity","exp","square","sqrt","log")
lnk <- match(link <- match.arg(link,tmp),tmp)
rf <- !is.null(mu)
sf <- !is.null(shape)
#
# check if a data object is being supplied
#
type <- "unknown"
respenv <- exists(deparse(substitute(response)),envir=parent.frame())&&
	inherits(response,"repeated")&&!inherits(envir,"repeated")
if(respenv){
	if(dim(response$response$y)[2]>1)
		stop("kalseries only handles univariate responses")
	if(!is.null(response$NAs)&&any(response$NAs))
		stop("kalseries does not handle data with NAs")
	type <- response$response$type}
envname <- if(respenv)deparse(substitute(response))
	else if(!is.null(class(envir)))deparse(substitute(envir))
	else NULL
#
# if envir, remove extra (multivariate) responses
#
if(!respenv&&inherits(envir,"repeated")){
	if(!is.null(envir$NAs)&&any(envir$NAs))
		stop("kalseries does not handle data with NAs")
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
#
# find the covariates
#
if(respenv||inherits(envir,"repeated")){
	if(rf)resp <- rmna(response$response)
	else {
		resp <- response
		# time-constant covariates
		if(is.null(ccov))resp$ccov <- NULL
		else if(inherits(ccov,"formula")){
			if(any(is.na(match(rownames(attr(terms(ccov,data=response),"factors")),colnames(response$ccov$ccov)))))
				stop("ccov covariate(s) not found")
			tmp <- wr(ccov,data=response,expand=FALSE)$design
			resp$ccov$ccov <- tmp[,-1,drop=FALSE]
			rm(tmp)}
		else stop("ccov must be a W&R formula")
		# time-varying covariates
		if(is.null(tvcov))resp$tvcov <- NULL
		else if(inherits(tvcov,"formula")){
			if(any(is.na(match(rownames(attr(terms(tvcov,data=response),"factors")),colnames(response$tvcov$tvcov)))))
				stop("tvcov covariate(s) not found")
			tmp <- wr(tvcov,data=response)$design
			resp$tvcov$tvcov <- tmp[,-1,drop=FALSE]
			rm(tmp)}
		else stop("tvcov must be a W&R formula")}
	nccov <- if(rf||is.null(resp$ccov$ccov)) 0
		 else  dim(resp$ccov$ccov)[2]
	ttvc <- if(rf||is.null(resp$tvcov$tvcov)) 0
		 else  dim(resp$tvcov$tvcov)[2]}
else {
	# make response object
	if(inherits(response,"response")){
		resp <- response
		type <- response$type}
	else resp <- restovec(response,times,delta=delta)
	# time-constant covariates
	if(is.null(ccov))nccov <- 0
	else {
		if(!inherits(ccov,"tccov")){
			ccname <- deparse(substitute(ccov))
			if((is.matrix(ccov)&&is.null(colnames(ccov)))){
				ccname <- deparse(substitute(ccov))
				if(dim(ccov)[2]>1){
					tmp <- NULL
					for(i in 1:dim(ccov)[2])
						tmp <- c(tmp,paste(ccname,i,sep=""))
					ccname <- tmp}}
			ccov <- tcctomat(ccov,names=ccname)}
		nccov <- if(rf) 0 else dim(ccov$ccov)[2]}
	# time-varying covariates
	if(is.null(tvcov))ttvc <- 0
	else {
		if(!inherits(tvcov,"tvcov")){
			tvcname <- deparse(substitute(tvcov))
			if(is.list(tvcov)&&dim(tvcov[[1]])[2]>1){
				if(is.null(colnames(tvcov[[1]]))){
					tvcname <- deparse(substitute(tvcov))
					tmp <- NULL
					for(i in 1:dim(tvcov[[1]])[2])
						tmp <- c(tmp,paste(tvcname,i,sep=""))
					tvcname <- tmp}
				else tvcname <- colnames(tvcov[[1]])}
			tvcov <- tvctomat(tvcov,names=tvcname)}
		ttvc <- if(rf) 0 else dim(tvcov$tvcov)[2]}
	resp <- rmna(response=resp, tvcov=tvcov, ccov=ccov)
	if(!is.null(ccov))rm(ccov)
	if(!is.null(tvcov))rm(tvcov)}
n <- dim(resp$response$y)[1]
nobs <- nobs(resp)
nind <- length(nobs)
if((inherits(envir,"repeated")&&(length(nobs)!=length(nobs(envir))||
	any(nobs!=nobs(envir))))||(inherits(envir,"tvcov")&&
	(length(nobs)!=length(envir$tvcov$nobs)||any(nobs!=envir$tvcov$nobs))))
	stop("response and envir objects are incompatible")
if((intensity=="exponential"||intensity=="Weibull"||intensity=="gamma"
	||intensity=="log normal"||intensity=="log logistic"
	||intensity=="log Cauchy"||intensity=="log Laplace"
	||intensity=="inverse Gauss")&&type!="unknown"&&type!="duration"
	&&type!="continuous")stop("duration data required")
if((intensity=="gen logistic"||intensity=="normal"||intensity=="logistic"
	||intensity=="Cauchy"||intensity=="Laplace")&&type!="unknown"
	&&type!="continuous"&&type!="duration")stop("continuous data required")
if(!is.null(resp$response$censor))stop("kalseries does not handle censoring")
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
# transform location formula to function and check number of parameters
#
npreg <- length(preg)
mu1 <- sh1 <- v <- b <- NULL
if(inherits(mu,"formula")){
	pr <- if(npreg>0)preg else ptvc
	npr <- length(pr)
	mu2 <- if(respenv)
		finterp(mu,.envir=response,.name=envname,.expand=is.null(preg))
		else finterp(mu,.envir=envir,.name=envname,.expand=is.null(preg))
	npt1 <- length(attr(mu2,"parameters"))
	if(is.character(attr(mu2,"model"))){
	# W&R formula
		if(length(attr(mu2,"model"))==1){
			mu1 <- function(p) p[1]*rep(1,n)
			attributes(mu1) <- attributes(mu2)
			mu2 <- NULL}}
	else {
	# formula with unknowns
		if(npr!=npt1){
			cat("\nParameters are ")
			cat(attr(mu2,"parameters"),"\n")
			stop(paste("preg or ptvc should have",npt1,"estimates"))}
		if(is.list(pr)){
			if(!is.null(names(pr))){
				o <- match(attr(mu2,"parameters"),names(pr))
				pr <- unlist(pr)[o]
				if(sum(!is.na(o))!=length(pr))
					stop("invalid estimates for mu - probably wrong names")}
			else pr <- unlist(pr)
			if(npreg>0)preg <- pr else ptvc <- pr}}
	if(!is.null(mu2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			mu1 <- function(p) mu2(p)[cv]
			attributes(mu1) <- attributes(mu2)}
		else {
			mu1 <- mu2
			rm(mu2)}}}
else if(is.function(mu))mu1 <- mu
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
# check that correct number of estimates was supplied
#
nlp <- if(is.function(mu1))length(attr(mu1,"parameters"))
	else if(is.null(mu1))NULL
	else npt1
if(!is.null(nlp)){
	if(is.null(ptvc)&&nlp!=npreg)
		stop(paste("preg should have",nlp,"initial estimates"))
	else if(!is.null(ptvc)&&length(ptvc)!=nlp)
		stop(paste("ptvc should have",nlp,"initial estimates"))}
#
# transform shape formula to function and check number of parameters
#
nps <- length(pshape)
if(inherits(shape,"formula")){
	sh2 <- if(respenv)finterp(shape,.envir=response,.name=envname)
		else finterp(shape,.envir=envir,.name=envname)
	npt2 <- length(attr(sh2,"parameters"))
	# W&R formula
	if(is.character(attr(sh2,"model"))){
		if(length(attr(sh2,"model"))==1){
			sh1 <- function(p) p[npl1]*rep(1,n)
			attributes(sh1) <- attributes(sh2)
			sh2 <- NULL}}
	else {
	# formula with unknowns
		if(nps!=npt2){
			cat("\nParameters are ")
			cat(attr(sh2,"parameters"),"\n")
			stop(paste("pshape should have",npt2,"estimates"))}
		if(is.list(pshape)){
			if(!is.null(names(pshape))){
				o <- match(attr(sh2,"parameters"),names(pshape))
				pshape <- unlist(pshape)[o]
				if(sum(!is.na(o))!=length(pshape))
					stop("invalid estimates for shape - probably wrong names")}
			else pshape <- unlist(pshape)}}
	if(!is.null(sh2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			sh1 <- function(p) sh2(p)[cv]
			attributes(sh1) <- attributes(sh2)}
		else {
			sh1 <- sh2
			rm(sh2)}}}
else if(is.function(shape))sh1 <- shape
#
# give appropriate attributes to sh1 for printing
#
if(!is.null(sh1)&&is.null(attr(sh1,"parameters")))
	attributes(sh1) <- if(is.function(shape)){
		if(!inherits(shape,"formulafn")){
			if(respenv)attributes(fnenvir(shape,.envir=response))
			else attributes(fnenvir(shape,.envir=envir))}
		else attributes(shape)}
		else {
			if(respenv)attributes(fnenvir(sh1,.envir=response))
			attributes(fnenvir(sh1,.envir=envir))}
#
# check that correct number of estimates was supplied
#
nlp <- if(is.function(shape))length(attr(sh1,"parameters"))
	else if(is.null(shape))NULL
	else npt2
if(!is.null(nlp)&&nlp!=nps)
	stop(paste("pshape should have",nlp,"initial estimates"))
#
# check that appropriate functions supplied
#
if(rf&&!is.function(mu1))stop("mu must be a formula or function")
if(sf&&!is.function(sh1))stop("shape must be a formula or function")
tvc <- length(ptvc)
if(!rf&&(ttvc>0&&tvc!=ttvc||ttvc==0&&tvc>0))
	stop(paste(ttvc,"initial estimates of coefficients for time-varying covariates must be supplied"))
if(intensity=="exponential"){
	sf <- FALSE
	pshape <- NULL}
else {
	if(is.null(pshape))
		stop("Initial value of the shape parameter must be supplied")
	if(!sf){
		if(pshape<=0)stop("Shape must be positive")
		else pshape <- log(pshape)}}
if(intensity=="gen logistic"){
	if(is.null(pintercept))
		stop("Initial value of the intercept parameter must be supplied")}
else pintercept <- NULL
if(pinitial<=0)stop("Estimate of initial parameter must greater than 0")
else pinitial <- log(pinitial)
#
# check what dependence is required
#
frser <- FALSE
if(depend=="independence")pdepend <- NULL
else if(depend=="frailty"){
	frser <- !is.null(pdepend)
	if(frser){
		if(pdepend<=0)stop("Dependence parameter must be positive")
		pdepend <- log(pdepend)}}
else {
	if(is.null(pdepend))
		stop("An estimate of the dependence parameter must be supplied")
	else if(pdepend<=0||pdepend>=1)
		stop("Dependence parameter must be between 0 and 1")
	else pdepend <- log(pdepend/(1-pdepend))}
#
# check times
#
if(is.null(resp$response$times)){
	if(depend=="serial"||depend=="Markov")
		stop("No times. Serial and Markov dependence cannot be fitted.")
	ave <- times <- 0}
else {
	if(torder>length(unique(resp$response$times)))
		stop("torder is too large for the number of distinct times")
	ave <- mean(resp$response$times)
	times <- resp$response$times-ave}
#
# check interactions with time
#
if(!is.null(interaction)){
	if(length(interaction)!=nccov)
		stop(paste(nccov,"interactions with time must be specified"))
	else if(any(interaction>torder))
		stop(paste("Interactions can be at most of order ",torder))}
else interaction <- rep(0,nccov)
#
# check number of parameter estimates
#
if(!is.null(pfamily)&&depend=="frailty")
		stop("pfamily cannot be used with frailty dependence")
npv <- torder+sum(interaction)
if(rf&&npreg>0)nccov <- npreg-1
if(!rf&&1+nccov+npv!=npreg)
	stop(paste(1+nccov+npv,"regression estimates must be supplied"))
y <- resp$response$y
#
# check Jacobian and transformations
if(!is.null(resp$response$delta))jacob <- if(length(resp$response$delta)==1)
		-n*log(resp$response$delta)
 	else -sum(log(resp$response$delta))
else jacob <- 0
if((mdl<=3||mdl>=8)&&any(y<=0))stop("All responses must be positive")
if(transform=="exp"){
	jacob <- jacob-sum(y)
	y <- exp(y)}
else if(transform!="identity"){
	if(any(y==0))stop("Zero response values: invalid transformation")
	if(transform=="square"){
		jacob <- jacob-sum(log(abs(y)))-n*log(2)
		y  <- y^2}
	if(any(y<0))stop("Nonpositive response values: invalid transformation")
	if(transform=="sqrt"){
		jacob <- jacob+sum(log(y))/2+n*log(2)
		y <- sqrt(y)}
	else if(transform=="log"){
		jacob <- jacob+sum(log(y))
		y <- log(y)}}
#
# set up for location function and make sure it gives correct number
# of values
#
if(rf){
	if(tvc>0&&nccov>0)
		stop("With a mean function or formula, initial estimates must be supplied either in preg or in ptvc")
	if(tvc>0){
		if(length(mu1(ptvc))!=n)
			stop("The mu function or formula must provide an estimate for each observation")
		tvc <- tvc-1}
	else {
		lp <- length(mu1(preg))
		if(lp==1){
			if(nccov==0)mu1 <- function(p) rep(p[1],nind)
			else stop("Number of estimates does not correspond to mu function")}
		else if(lp!=nind)
			stop("The mu function or formula must provide an estimate for each individual")}}
if(sf&&length(sh1(pshape))!=n)
	stop("The shape function must provide an estimate for each observation")
#
# check that the likelihood function gives appropriate values and call nlm
#
np <- 1+nccov+npv+tvc+1+(depend=="serial"||depend=="Markov")+
	nps+(!is.null(pintercept))+frser+(!is.null(pfamily))
nps1 <- np-nps+1
p <- c(preg,ptvc,pinitial,pdepend,pfamily,pintercept,pshape)
if(dep==3)serie <- serief
else serie <- series
if(fscale==1)fscale <- serie(p)
if(is.na(serie(p)))
	stop("Likelihood returns NAs: probably invalid initial values")
z0 <- nlm(serie, p=p, hessian=TRUE, print.level=print.level,
	typsize=typsize, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
p <- z0$estimate
like <- z0$minimum+jacob
#
# calculate se's
#
a <- if(any(is.na(z0$hessian))||any(abs(z0$hessian)==Inf))0
	else qr(z0$hessian)$rank
if(a==np)cov <- solve(z0$hessian)
else cov <- matrix(NA,ncol=np,nrow=np)
se <- sqrt(diag(cov))
corr <- cov/(se%o%se)
dimnames(corr) <- list(1:np,1:np)
#
# calculate recursive fitted values
#
if(mdl==4)z <- list()
else {
	z <- if(depend=="frailty"){
		if(rf)b <- mu1(p)
		if(sf)v <- sh1(p[nps1:np])
		.C("krand",
			p=as.double(p),
			y=as.double(y),
			t=as.double(times),
			x=as.double(resp$ccov$ccov),
			nind=as.integer(nind),
			nobs=as.integer(nobs),
			nbs=as.integer(n),
			nccov=as.integer(nccov),
			npv=as.integer(npv),
			model=as.integer(mdl),
			link=as.integer(lnk),
			density=as.integer(density),
			torder=as.integer(torder),
			inter=as.integer(interaction),
			tvc=as.integer(tvc),
			tvcov=resp$tvcov$tvcov,
			fit=as.integer(1),
			pred=double(n),
			rpred=double(n),
			rf=as.integer(rf),
			bbb=as.double(b),
			sf=as.integer(sf),
			vv=as.double(v),
			frser=as.integer(frser),
			like=double(1),
			DUP=FALSE,
			PACKAGE="repeated")}
	else {
		if(rf)b <- mu1(p)
		if(sf)v <- sh1(p[nps1:np])
		z <- .C("kserie",
			p=as.double(p),
			y=as.double(y),
			t=as.double(times),
			x=as.double(resp$ccov$ccov),
			nind=as.integer(nind),
			nobs=as.integer(nobs),
			nbs=as.integer(n),
			nccov=as.integer(nccov),
			npv=as.integer(npv),
			model=as.integer(mdl),
			link=as.integer(lnk),
			density=as.integer(density),
			pfamily=as.integer(!is.null(pfamily)),
			dep=as.integer(dep),
			torder=as.integer(torder),
			inter=as.integer(interaction),
			tvc=as.integer(tvc),
			tvcov=resp$tvcov$tvcov,
			fit=as.integer(1),
			pred=double(n),
			rpred=double(n),
			rf=as.integer(rf),
			bbb=as.double(b),
			sf=as.integer(sf),
			vv=as.double(v),
			like=double(1),
			DUP=FALSE,
			PACKAGE="repeated")}}
if(rf&&tvc>0){
	nccov <- tvc
	tvc <- 0}
#
# return appropriate attributes on functions
#
if(!is.null(mu3))mu1 <- mu3
if(!is.null(sh3))sh1 <- sh3
z <- list(
	call=call,
	intensity=intensity,
	pfamily=!is.null(pfamily),
	mdl=mdl,
	frser=frser,
	mu=mu1,
	npr=1+nccov+tvc+torder+sum(interaction),
	shape=sh1,
	nps=nps,
	density=density,
	depend=depend,
	torder=torder,
	interaction=interaction,
	response=resp$response,
	pred=z$pred,
	rpred=z$rpred,
	transform=transform,
	ccov=resp$ccov,
	tvcov=resp$tvcov,
	link=link,
	maxlike=like,
	aic=like+np,
	df=n-np,
	npt=np,
	npv=npv,
	coefficients=p,
	se=se,
	cov=cov,
	corr=corr,
	grad=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z) <- if(mdl==4)"kalseries" else c("kalseries","recursive")
return(z)}

### standard methods
###

deviance.kalseries <- function(z) 2*z$maxlike

fitted.kalseries <- function(z, recursive=TRUE)
if(recursive) z$rpred else z$pred

residuals.kalseries <- function(z, recursive=TRUE){
if(z$transform=="exp")z$response$y <- exp(z$response$y)
else if(z$transform=="square")z$response$y  <- z$response$y^2
else if(z$transform=="sqrt")z$response$y <- sqrt(z$response$y)
else if(z$transform=="log")z$response$y <- log(z$response$y)
res <- if(recursive)z$response$y-z$rpred else z$response$y-z$pred
class(res) <- "residuals"
res}

### print method
###
print.kalseries <- function(z,digits=max(3,.Options$digits-3),correlation=TRUE){
if(!is.null(z$ccov))nccov <- dim(z$ccov$ccov)[2]
else nccov <- 0
expm <- z$intensity!="exponential"&&!is.function(z$shape)
glm <- z$intensity=="gen logistic"
npt <- if(is.function(z$shape)) z$npt-z$nps else z$npt
deppar <- (z$depend=="serial"||z$depend=="Markov")
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
cat("Number of subjects    ",length(nobs(z)),"\n")
cat("Number of observations",length(z$response$y),"\n")
if(!is.null(z$response$time))
	cat("Times centred at      ",mean(z$response$time),"\n")
cat("Transformation        ",z$trans,"\n")
cat("Link function         ",z$link,"\n\n")
if(z$density)cat(z$intensity," density",sep="")
else cat(z$intensity," intensity",sep="")
if(z$depend=="independence")cat(" with independence\n")
else cat(" with",z$depend,"dependence",if(z$frser)"and AR","\n")
cat("\n-Log likelihood   ",z$maxlike,"\n")
cat("Degrees of freedom",z$df,"\n")
cat("AIC               ",z$aic,"\n")
cat("Iterations        ",z$iterations,"\n\n")
cat("Location parameters\n")
if(!is.null(attr(z$mu,"formula")))
	cat(deparse(attr(z$mu,"formula")),sep="\n")
else if(!is.null(attr(z$mu,"model"))){
	t <- deparse(attr(z$mu,"model"))
	t[1] <- sub("expression\\(","",t[1])
	t[length(t)] <- sub("\\)$","",t[length(t)])
	cat(t,sep="\n")}
coef.table <- cbind(z$coef[1:z$npr],z$se[1:z$npr])
if(inherits(z$mu,"formulafn"))
	cname <- if(is.character(attr(z$mu,"model")))attr(z$mu,"model")
		else attr(z$mu,"parameters")
else {
	cname <- "(Intercept)"
	if(nccov)cname <- c(cname,colnames(z$ccov$ccov))
	if(z$torder){
		cname <- c(cname,paste("t^",1:z$torder,sep=""))
		if(length(z$interaction)>0){
			for(i in 1:nccov)if(z$interaction[i]>0){
				cname <- c(cname,paste(colnames(z$ccov$ccov)[i],".t^",1:z$interaction[i],sep=""))}}}
	if(!is.null(z$tvcov))cname <- c(cname,colnames(z$tvcov$tvcov))}
dimnames(coef.table) <- list(cname, c("estimate","se"))
print.default(coef.table, digits=digits, print.gap=2)
if(is.function(z$shape))cat("\nDependence parameters\n")
else cat("\nNonlinear parameters\n")
coef <- exp(z$coef[(npt-deppar-z$frser-z$pfamily-glm-expm):npt])
cname <- "initial"
if(deppar){
	cname <- c(cname,"depend")
	coef[2] <- coef[2]/(1+coef[2])}
if(z$frser)cname <- c(cname,"AR")
if(z$pfamily){
	cname <- c(cname,"family")
	coef[length(coef)-glm-expm] <- z$coef[npt-glm-expm]}
if(glm){
	cname <- c(cname,"intercept")
	coef[length(coef)-expm] <- z$coef[npt-expm]
	if(expm){
		cname <- c(cname,"asymptote")
		coef[length(coef)] <- 1/coef[length(coef)]}}
else if(expm)cname <- c(cname,"shape")
coef.table <- cbind(z$coef[(npt-deppar-z$frser-z$pfamily-glm-expm):npt],
	z$se[(npt-deppar-z$frser-z$pfamily-glm-expm):npt],coef)
dimnames(coef.table) <- list(cname, c("estimate","se","parameter"))
print.default(coef.table, digits=digits, print.gap=2)
if(z$depend=="frailty"){
	tmp <- trigamma(exp(-z$coef[npt-deppar-expm-z$frser]))
	cat("Correlation =",tmp/(tmp+trigamma(1)),"\n")}
if(inherits(z$shape,"formulafn")){
	cat("\nShape parameters\n")
	if(!is.null(attr(z$shape,"formula")))
		cat(deparse(attr(z$shape,"formula")),sep="\n")
	else if(!is.null(attr(z$shape,"model"))){
		t <- deparse(attr(z$shape,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	cname <- if(is.character(attr(z$shape,"model")))
			attr(z$shape,"model")
		else attr(z$shape,"parameters")
	coef.table <- cbind(z$coef[(z$npt-z$nps+1):z$npt],
		z$se[(z$npt-z$nps+1):z$npt])
	dimnames(coef.table) <- list(cname, c("estimate","se"))
	print.default(coef.table, digits=digits, print.gap=2)}
if(correlation){
	cat("\nCorrelation matrix\n")
	print.default(z$corr, digits=digits)}}
