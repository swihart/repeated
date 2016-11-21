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
#     gar(response=NULL, distribution="normal", times=NULL, totals=NULL,
#	censor=NULL, delta=NULL, mu=NULL, shape=NULL, depend=NULL,
#	shfn=FALSE, common=FALSE, preg=NULL, pshape=NULL, pdepend=NULL,
#	parch=NULL, arch="square", transform="identity", link="identity",
#	autocorr="exponential", order=1, envir=parent.frame(), print.level=0,
#	ndigit=10, gradtol=0.00001, steptol=0.00001, fscale=1, iterlim=100,
#	typsize=abs(p), stepmax=10*sqrt(p%*%p))
#
#  DESCRIPTION
#
#    A function to fit generalized nonlinear autoregression models with
#  various distributions

gar <- function(response=NULL, distribution="normal", times=NULL, totals=NULL,
	censor=NULL, delta=NULL, mu=NULL, shape=NULL, depend=NULL, shfn=FALSE,
	common=FALSE, preg=NULL, pshape=NULL, pdepend=NULL, parch=NULL,
	arch="square", transform="identity", link="identity",
	autocorr="exponential", order=1, envir=parent.frame(), print.level=0,
	ndigit=10, gradtol=0.00001, steptol=0.00001, fscale=1, iterlim=100,
	typsize=abs(p), stepmax=10*sqrt(p%*%p)){
#
# set up likelihood to call C code
#
likekal <- function(p){
	eta <- mu1(p)
	if(sh){
		shr <- sh1(p)
		theta <- exp(p[np])
		if(!tm)theta <- c(p[npr1:nprd],theta)}
	else theta <- if(tm)p[nprd1:np] else p[npr1:np]
	if(tm)dep <- dep1(p)
	if(archtype>0)parch <- p[npa]
	z <- .C("gar",
		y=as.double(y),
		total=as.double(response$response$n),
		my=as.integer(3*max(y)),
		nobs=as.integer(nobs(response)),
		nind=as.integer(nind),
		times=as.double(response$response$times),
		censor=as.integer(censor),
		cens=as.integer(!is.null(response$response$censor)),
		eta=as.double(eta),
		tm=as.integer(tm),
		theta=as.double(theta),
		dep=as.double(dep),
		model=as.integer(mdl),
		thp=as.integer(thp),
		shape=as.double(shr),
		sh=as.integer(sh),
		link=as.integer(lnk),
		ar=as.integer(ar),
		order=as.integer(order),
		parch=as.double(parch),
		arch=as.integer(archtype),
		pred=double(n),
		rpred=double(n),
		volatility=double(n),
		like=double(1),
		DUP=FALSE,
		PACKAGE="repeated")
	z$like+jacob}
call <- sys.call()
#
# check model required
#
tmp <- c("binomial","Poisson","exponential","negative binomial",
	"mult Poisson","double Poisson","Consul","beta binomial",
	"mult binomial","double binomial","normal","logistic","Cauchy",
	"Weibull","gamma","Laplace","inverse Gauss","Pareto","Levy",
	"beta","simplex","two-sided power","gen gamma","gen logistic","Hjorth",
	"Burr","gen Weibull","gen extreme value","gen inverse Gauss",
	"power exponential","power variance function Poisson",
	"skew Laplace","Student t")
mdl <- match(distribution <- match.arg(distribution,tmp),tmp)
bmdl <- mdl==1||mdl==8||mdl==9||mdl==10
par1 <- distribution=="binomial"||distribution=="Poisson"||
	distribution=="exponential"
par3 <- mdl>22
tmp <- c("exponential","gaussian","cauchy","spherical","IOU")
ar <- match(autocorr <- match.arg(autocorr,tmp),tmp)
#
# check transformation and link function
#
transform <- match.arg(transform,c("identity","exp","square","sqrt","log"))
if(transform!="identity"&&((mdl<11&&mdl!=3)||mdl==31))
	stop("transformations not valid for discrete data models")
tmp <- c("identity","exp","square","sqrt","log","logit","cloglog","loglog")
lnk <- match(link <- match.arg(link,tmp),tmp)
if((link=="logit"||link=="cloglog"||link=="loglog")&&!bmdl)
	stop("logit, cloglog, and loglog links can only be used with binary data")
#
# check that either AR or ARCH is being fitted
#
if(is.null(pdepend)&&is.null(parch))
	stop("Either pdepend or parch must be supplied")
if(!is.null(parch)){
	if(par1)stop(paste("arch model not possible for",distribution,"distribution"))
	if(length(parch)!=1)stop("parch should contain one value")
	tmp <- c("square","absolute value","exponential")
	archtype <- match(match.arg(arch,tmp),tmp)
	if(archtype!=3&&parch<=0)
		stop("parch must be positive for square and absolute value models")}
else archtype <- 0
#
# check if a data object is being supplied
#
type <- "unknown"
respenv <- exists(deparse(substitute(response)),envir=parent.frame())&&
	inherits(response,"repeated")
if(respenv){
	if(dim(response$response$y)[2]>1)
		stop("gar only handles univariate responses")
	if(!is.null(response$NAs)&&any(response$NAs))
		stop("gar does not handle data with NAs")
	type <- response$response$type}
envname <- if(respenv)deparse(substitute(response))
	else if(!is.null(class(envir)))deparse(substitute(envir))
	else NULL
#
#if response is not a data object, find or make one
#
if(!respenv){
	if(inherits(envir,"repeated")){
	# if envir, remove extra (multivariate) responses
		if(!is.null(envir$NAs)&&any(envir$NAs))
			stop("gar does not handle data with NAs")
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
				if(all(response$response$delta==1)||all(is.na(response$response$delta)))response$response$delta <- NULL}
			if(!is.null(response$response$n)){
				response$response$n <- response$response$n[,col,drop=FALSE]
				if(all(is.na(response$response$n)))response$response$n <- NULL}
			if(!is.null(envir$response$censor)){
				response$response$censor <- response$response$censor[,col,drop=FALSE]
				if(all(is.na(response$response$censor)))response$response$censor <- NULL}}}
	else {
		if(!inherits(response,"response")){
			if(is.null(times)){
				if(is.matrix(response))times <- matrix(1:dim(response)[2],nrow=dim(response)[1],ncol=dim(response)[2],byrow=T)
				else if(is.vector(response))times <- 1:length(response)}
			response <- if(bmdl)
				restovec(response, times=times, totals=totals)
				else restovec(response, times=times, censor=censor, delta=delta)}
		type <- response$type
		response <- rmna(response=response)}}
if(is.null(times(response)))stop("No times available")
if((inherits(envir,"repeated")&&
	(length(nobs(response))!=length(nobs(envir))||
	any(nobs(response)!=nobs(envir))))||(inherits(envir,"tvcov")&&
	(length(nobs(response))!=length(envir$tvcov$nobs)||
	any(nobs(response)!=envir$tvcov$nobs))))
	stop("response and envir objects are incompatible")
if(bmdl){
	if(type!="unknown"&&type!="nominal")stop("nominal data required")
	if(is.null(response$response$n)){
		if(any(response$response$y!=0&&response$response$y!=1))
			stop("responses must be binary if totals are not supplied")
		else response$response$n <- rep(1,dim(response$response$y)[1])}}
y <- response$response$y
n <- dim(y)[1]
#
# if a data object was supplied, modify formulae or functions to read from it
#
mu3 <- sh3 <- dep3 <- dep <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")||inherits(envir,"tvcov")){
	if(is.null(envname))envname <- deparse(substitute(envir))
	# location function or formula
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
	# shape function or formula
	if(inherits(shape,"formula")){
		sh3 <- if(respenv)finterp(shape,.envir=response,.name=envname)
			else finterp(shape,.envir=envir,.name=envname)}
	else if(is.function(shape)){
		if(is.null(attr(shape,"model"))){
		        tmp <- parse(text=deparse(shape)[-1])
		        shape <- if(respenv)fnenvir(shape,.envir=response,.name=envname)
		        	else fnenvir(shape,.envir=envir,.name=envname)
		        sh3 <- shape
			if(shfn){
				pos <- grep("mu",attr(sh3,"parameters"))
				attr(sh3,"parameters") <- attr(sh3,"parameters")[-pos]}
		        attr(sh3,"model") <- tmp}
		else sh3 <- shape}
	# AR function or formula
	if(inherits(depend,"formula")){
		dep3 <- if(respenv)finterp(depend,.envir=response,.name=envname)
			else finterp(depend,.envir=envir,.name=envname)}
	else if(is.function(depend)){
		if(is.null(attr(depend,"model"))){
		        tmp <- parse(text=deparse(depend)[-1])
		        depend <- if(respenv)fnenvir(depend,.envir=response,.name=envname)
		        	else fnenvir(depend,.envir=envir,.name=envname)
		        dep3 <- depend
		        attr(dep3,"model") <- tmp}
		else dep3 <- depend}}
#
# transform location formula to function and check number of parameters
#
npr <- length(preg)
npr1a <- npr1 <- npr+1
if(inherits(mu,"formula")){
	mu2 <- if(respenv)finterp(mu,.envir=response,.name=envname)
		else finterp(mu,.envir=envir,.name=envname)
	npt1 <- length(attr(mu2,"parameters"))
	if(is.character(attr(mu2,"model"))){
	# W&R formula
		if(length(attr(mu2,"model"))==1){
		# intercept model
			mu1 <- function(p) p[1]*rep(1,n)
			attributes(mu1) <- attributes(mu2)
			mu2 <- NULL}}
	else {
	# formula with unknowns
		if(npr!=npt1&&!common){
			cat("\nParameters are ")
			cat(attr(mu2,"parameters"),"\n")
			stop(paste("preg should have",npt1,"estimates"))}
		if(is.list(preg)){
			if(!is.null(names(preg))){
				o <- match(attr(mu2,"parameters"),names(preg))
				preg <- unlist(preg)[o]
				if(sum(!is.na(o))!=length(preg))stop("invalid estimates for mu - probably wrong names")}
			else preg <- unlist(preg)}}
	if(!is.null(mu2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			mu1 <- function(p) mu2(p)[cv]
			attributes(mu1) <- attributes(mu2)}
		else {
			mu1 <- mu2
			rm(mu2)}}}
else if(is.null(mu)){
	mu1 <- function(p) p[1]*rep(1,n)
	npt1 <- 1}
else mu1 <- mu
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
if(!is.null(nlp)&&!common&&nlp!=npr)
	stop(paste("preg should have",nlp,"initial estimates"))
#
# transform shape formula to function and check number of parameters
#
nprd <- npr+length(pdepend)
nprd1 <- if(common&&!inherits(shape,"formula")) 1 else nprd+1
nps1 <- nps <- length(pshape)
if(par3)nps1 <- nps1-1
if(inherits(shape,"formula")){
	old <- if(common)mu1 else NULL
	mufn <- if(shfn)"mu" else NULL
	sh3 <- if(respenv)finterp(shape,.envir=response,.start=nprd1,.name=envname,.old=old,.args=mufn)
		else finterp(shape,.envir=envir,.start=nprd1,.name=envname,.old=old,.args=mufn)
	tmp <- attributes(sh3)
	sh2 <- if(shfn)function(p) sh3(p,mu1(p)) else sh3
	attributes(sh2) <- tmp
	npt2 <- length(attr(sh2,"parameters"))
	if(is.character(attr(sh2,"model"))){
	# W&R formula
		if(length(attr(sh2,"model"))==1){
		# intercept model
			sh1 <- function(p) p[nprd1]*rep(1,n)
			attributes(sh1) <- attributes(sh2)
			sh2 <- NULL}}
	else {
	# formula with unknowns
		if(nps1!=npt2&&!common){
			cat("\nParameters are ")
			cat(attr(sh2,"parameters"),"\n")
			stop(paste("pshape should have",npt2+par3,"estimates"))}
		if(is.list(pshape)){
			if(!is.null(names(pshape))){
				o <- match(attr(sh2,"parameters"),names(pshape))
				pshape <- unlist(pshape)[o]
				if(sum(!is.na(o))!=length(pshape))stop("invalid estimates for shape - probably wrong names")}
			else pshape <- unlist(pshape)}}
	if(!is.null(sh2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			sh1 <- function(p) sh2(p)[cv]
			attributes(sh1) <- attributes(sh2)}
		else {
			sh1 <- sh2
			rm(sh2)}}}
else if(!is.function(shape)&&!par1){
	sh1 <- function(p) p[nprd1]*rep(1,n)
	npt2 <- 1
	if(length(pshape)!=1+par3)
		stop(paste("pshape must contain",1+par3,"estimate(s)"))}
else sh1 <- if(shfn)function(p) shape(p[nprd1:np], mu1(p))
	else function(p) shape(p[nprd1:np])
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
	else attributes(fnenvir(sh1,.envir=envir))}
#
# check that correct number of estimates was supplied
#
nlp <- if(is.function(shape))length(attr(sh1,"parameters"))-shfn
	else if(is.null(shape))NULL
	else npt2
if(!is.null(nlp)&&!common&&nlp!=nps1)
	stop(paste("pshape should have",nlp+par3,"initial estimates"))
if(common){
	nlp <- length(unique(c(attr(mu1,"parameters"),attr(sh1,"parameters"))))
	if(nlp!=npr)stop(paste("with a common parameter model, preg should contain",nlp,"estimates"))}
#
# transform AR formula to function and check number of parameters
#
npd <- length(pdepend)
if(inherits(depend,"formula")){
	dep2 <- if(respenv)finterp(depend,.envir=response,.start=npr1,.name=envname)
		else finterp(depend,.envir=envir,.start=npr1,.name=envname)
	npt3 <- length(attr(dep2,"parameters"))
	if(is.character(attr(dep2,"model"))){
	# W&R formula
		if(length(attr(dep2,"model"))==1){
		# intercept model
			dep1 <- function(p) p[npr1]*rep(1,n)
			attributes(dep1) <- attributes(dep2)
			dep2 <- NULL}}
	else {
	# formula with unknowns
		if(npd!=npt3){
			cat("\nParameters are ")
			cat(attr(dep2,"parameters"),"\n")
			stop(paste("pdepend should have",npt3,"estimates"))}
		if(is.list(pdepend)){
			if(!is.null(names(pdepend))){
				o <- match(attr(dep2,"parameters"),names(pdepend))
				pdepend <- unlist(pdepend)[o]
				if(sum(!is.na(o))!=length(pdepend))stop("invalid estimates for depend - probably wrong names")}
			else pdepend <- unlist(pdepend)}}
	if(!is.null(dep2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			dep1 <- function(p) dep2(p)[cv]
			attributes(dep1) <- attributes(dep2)}
		else {
			dep1 <- dep2
			rm(dep2)}}}
else if(is.function(dep))dep1 <- dep
else {
	dep1 <- NULL
	npt3 <- 1}
#
# give appropriate attributes to dep1 for printing
#
if(!is.null(dep1)&&is.null(attr(dep1,"parameters")))
	attributes(dep1) <- if(is.function(depend)){
		if(!inherits(depend,"formulafn")){
			if(respenv)attributes(fnenvir(depend,.envir=response))
			else attributes(fnenvir(depend,.envir=envir))}
		else attributes(depend)}
		else {
			if(respenv)attributes(fnenvir(dep1,.envir=response))
			else attributes(fnenvir(dep1,.envir=envir))}
#
# check that correct number of estimates was supplied
#
nlp <- if(is.function(depend))length(attr(dep1,"parameters"))
	else if(is.null(depend))NULL
	else npt3
if(!is.null(nlp)&&nlp!=npd)
	stop(paste("pdepend should have",nlp,"initial estimates"))
#
# check response and functions
#
if(mdl<11){
	if((distribution=="Poisson"||distribution=="negative binomial"
		||distribution=="mult Poisson"
		||distribution=="double Poisson"||distribution=="Consul")
		&&type!="unknown"&&type!="discrete")
		stop("discrete data required")
	if(distribution=="exponential"&&type!="unknown"&&type!="duration")
		stop("duration data required")
	if(any(y<0))stop("All responses must be non-negative")}
else if(mdl==20||mdl==21||mdl==22){
	if(type!="unknown"&&type!="continuous")stop("continuous data required")
	if(any(y<=0)||any(y>=1))stop("All responses must be between 0 and 1")}
else if(mdl==31){
	if(type!="unknown"&&type!="discrete")stop("discrete data required")
	if(any(y<0))stop("All responses must be non-negative")}
else if(mdl!=11&&mdl!=12&&mdl!=13&&mdl!=16&&mdl!=19&&mdl!=24&&mdl!=30&&
	mdl!=32&&mdl!=33){
	if(type!="unknown"&&type!="duration"&&type!="continuous")
		stop("duration data required")
	if(any(y<=0))stop("All responses must be positive")}
else {
	if(type!="unknown"&&type!="continuous"&&type!="duration")
		stop("continuous data required")
	if(distribution=="Levy"&&any(response$response$y<=mu1(preg)))
		stop("Location function must give values strictly less than corresponding observations")}
if(distribution=="Pareto"&&any(sh1(pshape)<=1))
	stop("shape function must be > 1")
#
# set up censoring indicator
#
censor <- response$response$censor
if(is.null(censor))censor <- rep(1,n)
else if(bmdl||mdl==2||mdl==4||mdl==5||mdl==6||mdl==7||mdl==31)
	stop(paste("Censored data not allowed for",distribution,"distribution"))
#
# set up transformation and Jacobian
#
nind <- length(nobs(response))
if(transform=="identity")jacob <- 0
else if(transform=="exp"){
	jacob <- -sum(y[censor==1])
	y <- exp(y)}
else if(any(y==0))stop("Zero response values: invalid transformation")
else if(transform=="square"){
	jacob <- -sum(log(abs(y[censor==1])))-n*log(2)
	y  <- y^2}
else if(any(y<0))stop("Nonpositive response values: invalid transformation")
else if(transform=="sqrt"){
	jacob <- sum(log(y[censor==1]))/2+n*log(2)
	y <- sqrt(y)}
else if(transform=="log"){
	jacob <- sum(log(y[censor==1]))
	y <- log(y)}
#
# set up unit of measurement
#
if(!is.null(response$response$delta)){
	if(length(response$response$delta)==1)
		jacob <- jacob-length(y[censor==1])*log(response$response$delta)
	else jacob <- jacob-sum(log(response$response$delta[censor==1]))}
#
# check AR model
#
if(is.null(pdepend))order <- 0
else if(order!=1&&order!=2)stop("Autoregression must have order 1 or 2")
tm <- !is.null(depend)
if(tm)order <- 0
if(order==2&&length(pdepend)!=2)
	stop("2 estimates of dependence parameters must be given")
else if(order>0&&!tm&&length(pdepend)!=1&&length(pdepend)!=2)
     stop("One or two estimates of dependence parameters must be given")
thp <- length(pdepend)==2&&order==1&&!tm
if(!tm&&any(pdepend<=0))stop("All dependence parameters must be positive")
#
# check location function
#
if(length(mu1(preg))!=dim(response$response$y)[1])
	stop("The mu function must provide an estimate for each observation")
else if(any(is.na(mu1(preg))))
	stop("Non-numerical mu: probably invalid initial values")
#
# set up parameter vector
#
p <- if(tm)c(preg,pdepend) else if(order>0)c(preg,-log(pdepend)) else preg
sh <- mdl>3&&(is.function(shape)||inherits(shape,"formula"))
if(!par1){
	if(!sh)p <- c(p,if(mdl==31)c(log(pshape[1]),pshape[2])else log(pshape))
	else {
		if(par3){
			if(pshape[nps]<=0&&mdl!=31)
				stop("extra shape parameter must be positive")
			p <- if(common)c(p,if(mdl==31)pshape[nps]
				else log(pshape[nps]))
				else c(p,pshape[1:(nps-1)],if(mdl==31)
				pshape[nps] else log(pshape[nps]))}
		else p <- c(p,pshape)}}
np <- length(p)
if(archtype>0)p <- c(p,if(archtype==1)log(parch)else parch)
npa <- length(p)
#
# check shape function
#
if(!sh){
	if((mdl<=3&&nps!=0)||(mdl>3&&mdl<23&&nps!=1)||(mdl>=23&&nps!=2))
		stop("Incorrect number of shape parameter estimates")
	else if(nps>0){
		if(mdl!=31&&any(pshape<=0)||mdl==31&&pshape[1]<=0)
			stop("All shape parameters must be positive")}
	shr <- rep(0,dim(response$response$y)[1])}
else {
	if(any(is.na(sh1(p))))
		stop("The shape model returns NAs: probably invalid initial values")
	if(length(sh1(p))!=dim(response$response$y)[1])
		stop("The shape function must provide an estimate for each observation")}
#
# check likelihood function and optimize with nlm
#
if(fscale==1)fscale <- likekal(p)
if(is.na(likekal(p)))
	stop("Likelihood returns NAs: probably invalid initial values")
z0 <- nlm(likekal, p, hessian=TRUE, print.level=print.level,
	typsize=typsize, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
like <- z0$minimum
#
# calculate se's
#
a <- if(any(is.na(z0$hessian))||any(abs(z0$hessian)==Inf))0
	else qr(z0$hessian)$rank
if(a==npa)cov <- solve(z0$hessian)
else cov <- matrix(NA,ncol=npa,nrow=npa)
se <- sqrt(diag(cov))
corr <- cov/(se%o%se)
dimnames(corr) <- list(1:npa,1:npa)
#
# back transform parameters
#
eta <- mu1(z0$estimate)
if(sh){
	shr <- sh1(z0$estimate)
	theta <- exp(z0$estimate[np])
	if(!tm)theta <- c(z0$estimate[npr1:nprd],theta)}
else theta <- if(tm)z0$estimate[nprd1:np] else z0$estimate[npr1:np]
if(tm)dep <- dep1(z0$estimate)
#
# calculate predicted values (and volatility)
#
if(archtype>0)parch <- z0$estmate[npa]
z <- .C("gar",
	y=as.double(y),
	total=as.double(response$response$n),
	my=as.integer(3*max(y)),
	nobs=as.integer(nobs(response)),
	nind=as.integer(nind),
	times=as.double(response$response$times),
	censor=as.integer(censor),
	cens=as.integer(!is.null(response$response$censor)),
	eta=as.double(eta),
	tm=as.integer(tm),
	theta=as.double(theta),
	dep=as.double(dep),
	model=as.integer(mdl),
	thp=as.integer(thp),
	shape=as.double(shr),
	sh=as.integer(sh),
	link=as.integer(lnk),
	ar=as.integer(ar),
	order=as.integer(order),
	parch=as.double(parch),
	arch=as.integer(archtype),
	pred=double(n),
	rpred=double(n),
	volatility=double(n),
	like=double(1),
	DUP=FALSE,
	PACKAGE="repeated")
if(archtype==0)z$volatility <- NULL
#
# return object
#
if(!is.null(mu3))mu1 <- mu3
if(!is.null(sh3))sh1 <- sh3
if(!is.null(dep3))dep1 <- dep3
z <- list(
	call=call,
	distribution=distribution,
	mu=mu1,
	formula=mu,
	shape=shape,
	depend=depend,
	sh1=sh1,
	shfn=shfn,
	common=common,
	dep1=dep1,
	response=response$response,
	link=link,
	order=order,
	autocorr=autocorr,
	transform=transform,
	maxlike=like,
	aic=like+npa,
	df=dim(response$response$y)[1]-npa,
	np=np,
	npa=npa,
	npr=npr,
	nps=nps,
	npd=npd,
	thp=thp,
	tm=tm,
	archtype=archtype,
	coefficients=z0$estimate,
	se=se,
	cov=cov,
	corr=corr,
	pred=z$pred,
	rpred=z$rpred,
	volatility=z$volatility,
	grad=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z) <- c("gar","recursive")
return(z)}

### standard methods
###

deviance.gar <- function(z) 2*z$maxlike

fitted.gar <- function(z, recursive=TRUE) if(recursive) z$rpred else z$pred

residuals.gar <- function(z, recursive=TRUE){
if(z$transform=="exp")z$response$y <- exp(z$response$y)
else if(z$transform=="square")z$response$y  <- z$response$y^2
else if(z$transform=="sqrt")z$response$y <- sqrt(z$response$y)
else if(z$transform=="log")z$response$y <- log(z$response$y)
if(recursive) z$response$y-z$rpred else z$response$y-z$pred}

### print method
###
print.gar <- function(z,digits=max(3,.Options$digits-3),correlation=TRUE){
np1 <- if(z$distribution=="binomial"||z$distribution=="exponential"
		||z$distribution=="Poisson") 0
	else if(z$distribution=="gen gamma"
		||z$distribution=="gen logistic"
		||z$distribution=="Hjorth"||z$distribution=="Burr"
		||z$distribution=="gen Weibull"
		||z$distribution=="gen extreme value"
		||z$distribution=="gen inverse Gauss"
		||z$distribution=="power exponential"
		||z$distribution=="power variance function Poisson"
		||z$distribution=="skew Laplace"
		||z$distribution=="Student t") 2
	else 1
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
cat("Number of subjects    ",length(nobs(z)),"\n")
cat("Number of observations",length(z$response$y),"\n")
cat("Transformation        ",z$trans,"\n")
cat("Link function         ",z$link,"\n\n")
cat(z$distribution,"distribution\n")
if(z$order<2)cat("First order ")
else cat("Second order ")
cat(z$autocorr,"dependence\n")
cat("\n-Log likelihood   ",z$maxlike,"\n")
cat("Degrees of freedom",z$df,"\n")
cat("AIC               ",z$aic,"\n")
cat("Iterations        ",z$iterations,"\n\n")
if(z$common)cat("Location model\n")
else cat("Location parameters\n")
if(!is.null(attr(z$mu,"formula")))
	cat(deparse(attr(z$mu,"formula")),sep="\n")
else if(!is.null(attr(z$mu,"model"))){
	t <- deparse(attr(z$mu,"model"))
	t[1] <- sub("expression\\(","",t[1])
	t[length(t)] <- sub("\\)$","",t[length(t)])
	cat(t,sep="\n")}
cname <- if(is.character(attr(z$mu,"model")))attr(z$mu,"model")
	else attr(z$mu,"parameters")
coef.table <- cbind(z$coef[1:z$npr],z$se[1:z$npr])
colname <- c("estimate","se")
if(!z$common){
	dimnames(coef.table) <- list(cname,colname)
	print.default(coef.table, digits=digits, print.gap=2)}
else {
	cat("\nShape model\n")
        if(!is.null(attr(z$sh1,"formula")))
        	cat(deparse(attr(z$sh1,"formula")),sep="\n")
        else if(!is.null(attr(z$sh1,"model"))){
        	t <- deparse(attr(z$sh1,"model"))
        	t[1] <- sub("expression\\(","",t[1])
        	t[length(t)] <- sub("\\)$","",t[length(t)])
        	cat(t,sep="\n")}
        cname <- c(cname,if(is.character(attr(z$sh1,"model")))
			attr(z$sh1,"model")
		else attr(z$sh1,"parameters")[1:(length(attr(z$sh1,"parameters"))-z$shfn)])
	cname <- unique(cname)
	cat("\nCommon parameters\n")
	dimnames(coef.table) <- list(cname,colname)
	print.default(coef.table, digits=digits, print.gap=2)}
if(z$order>0||z$tm){
	if(z$thp||z$order==2){
		cat("\nDependence parameters\n")
		if(z$thp)cname <- c("phi","rho")
		else cname <- c("rho1","rho2")
		coef.table <- cbind(z$coef[(z$npr+1):(z$npr+2)],
			z$se[(z$npr+1):(z$npr+2)],
			exp(-z$coef[(z$npr+1):(z$npr+2)]))}
	else if(z$tm){
		cat("\nDependence parameters\n")
		if(!is.null(attr(z$dep1,"formula")))
			cat(deparse(attr(z$dep1,"formula")),sep="\n")
		else if(!is.null(attr(z$dep1,"model"))){
			t <- deparse(attr(z$dep1,"model"))
			t[1] <- sub("expression\\(","",t[1])
			t[length(t)] <- sub("\\)$","",t[length(t)])
			cat(t,sep="\n")}
		cname <- if(is.character(attr(z$dep1,"model")))
			attr(z$dep1,"model")
			else attr(z$dep1,"parameters")
		coef.table <- cbind(z$coef[(z$npr+1):(z$npr+z$npd)],
			z$se[(z$npr+1):(z$npr+z$npd)])
		colname <- c("estimate","se")
		dimnames(coef.table) <- list(cname,colname)}
	else {
		cat("\nDependence parameter\n")
		cname <- "rho"
		coef.table <- cbind(z$coef[z$npr+1],
			z$se[z$npr+1],exp(-z$coef[z$npr+1]))}
	if(!z$tm)dimnames(coef.table) <- list(cname, c("estimate","se","parameter"))
	print.default(coef.table, digits=digits, print.gap=2)}
if(np1>0&&!z$common){
	cat("\nShape parameters\n")
	if(is.null(z$shape)){
		cname <- "shape"
		if(np1==2)cname <- c(cname,"psi")
		coef.table <- cbind(z$coef[(z$np-np1+1):z$np],
		z$se[(z$np-np1+1):z$np],if(z$distribution==
		"power variance function Poisson")
		c(exp(z$coef[(z$np-np1+1):(z$np-1)]),z$coef[z$np]) else 
		exp(z$coef[(z$np-np1+1):z$np]))
		dimnames(coef.table) <- list(cname, c("estimate","se","parameter"))}
	else {
                if(!is.null(attr(z$sh1,"formula")))
                	cat(deparse(attr(z$sh1,"formula")),sep="\n")
                else if(!is.null(attr(z$sh1,"model"))){
                	t <- deparse(attr(z$sh1,"model"))
                	t[1] <- sub("expression\\(","",t[1])
                	t[length(t)] <- sub("\\)$","",t[length(t)])
                	cat(t,sep="\n")}
                cname <- if(is.character(attr(z$sh1,"model")))
				attr(z$sh1,"model")
			else attr(z$sh1,"parameters")
		np2 <- length(cname)+(np1==2)
		coef.table <- cbind(z$coef[(z$np-np2+1):z$np],
			z$se[(z$np-np2+1):z$np])
                if(np1==2){
                	cname <- c(cname,"psi")
                	coef.table <- cbind(coef.table,
				c(coef.table[1:(dim(coef.table)[1]-1),1],
				if(z$distribution==
				"power variance function Poisson")
				z$coef[z$np] else exp(z$coef[z$np])))
                	colname <- c(colname,"parameter")}
		dimnames(coef.table) <- list(cname,colname)}
	print.default(coef.table,digits=digits,print.gap=2)}
if(z$archtype>0){
	cat(if(z$archtype==1)"\nSquare"else if(z$archtype==2)"\nAbsolute value"
		else "\nExponential","ARCH parameter\n")
	cname <- c("estimate","se",if(z$archtype!=3)"arch")
	coef.table <- cbind(z$coef[z$npa],z$se[z$npa],
		if(z$archtype!=3)exp(z$coef[z$npa]))
	dimnames(coef.table) <- list("", cname)
	print.default(coef.table,digits=digits,print.gap=2)}
if(correlation){
	cat("\nCorrelation matrix\n")
	print.default(z$corr, digits=digits)}}

### volatility method
###
volatility <- function(z, ...) UseMethod("volatility")

volatility.gar <- function(z, nind=NULL){
if(is.null(nind))nind <- 1:dim(z$response$y)[1]
else if(length(nind)>length(nobs(z))||any(nind>length(nobs(z))))
	stop("Individual not found")
else nind <- !is.na(match(covind(z),nind))
if(all(!nind))stop("No such individuals")
if(is.null(z$volatility))return(NULL)
z$volatility[nind]}
