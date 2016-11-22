#
#  repeated : A Library of Repeated Measurements Models
#  Copyright (C) 1999, 2000, 2001 J.K. Lindsey
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
#     cmcre(response, covariate=NULL, parameters, pcov=NULL, gradient=FALSE,
#	hessian=FALSE, print.level=0, ndigit=10, gradtol=0.00001,
#	steptol=0.00001, iterlim=100, fscale=1, typsize=abs(Params),
#	stepmax=Params)
#
#  DESCRIPTION
#
#    A function to fit a continuous-time two-state Markov process with a
#  random effect



##' Continuous-time Two-state Markov Processes with Random Effect
##' 
##' \code{cmcre} fits a two-state Markov process in continuous time, possibly
##' with one or two random effects and/or one covariate.
##' 
##' 
##' @param response A six-column matrix. Column 1: subject identification
##' (subjects can occupy several rows); column 2: time gap between events;
##' columns 3-6: transition matrix frequencies.
##' @param covariate An optional vector of length equal to the number of rows
##' of \code{response} upon which the equilibrium probability may depend.
##' @param parameters Initial parameter estimates. The number of them
##' determines the model fitted (minimum 2, yielding an ordinary Markov
##' process). 1: beta1=log(-log(equilibrium probability)); 2: beta2=log(sum of
##' transition intensities); 3: log(tau1)=log(random effect variance for
##' equilibrium probability); 4: log(tau2)=log(random effect variance for sum
##' of transition intensities).
##' @param pcov Initial parameter estimate for the covariate influencing the
##' equilibrium probability: exp(-exp(beta1+beta*covariate)).
##' @param gradient If TRUE, analytic gradient is used (with accompanying loss
##' of speed).
##' @param hessian If TRUE, analytic hessian is used (with accompanying loss of
##' speed).
##' @param others Arguments controlling \code{\link{nlm}}.
##' @return A list of class \code{cmcre} is returned.
##' @author R.J. Cook and J.K. Lindsey
##' @seealso \code{\link[repeated]{chidden}}, \code{\link[repeated]{hidden}}.
##' @references Cook, R.J. (1999) A mixed model for two-state Markov processes
##' under panel observations. Biometrics 55, 915-920.
##' @keywords models
##' @examples
##' 
##' # 12 subjects observed at intervals of 7 days
##' y <- matrix(c(1,7,1,2,3,5,
##' 	2,7,10,2,2,0,
##' 	3,7,7,0,1,1,
##' 	4,7,2,1,0,7,
##' 	5,7,1,1,1,11,
##' 	6,7,5,4,4,1,
##' 	7,7,1,1,1,8,
##' 	8,7,2,3,4,2,
##' 	9,7,9,0,0,0,
##' 	10,7,0,1,2,8,
##' 	11,7,8,2,2,1,
##' 	12,7,9,2,2,1),ncol=6, byrow=TRUE)
##' # ordinary Markov process
##' cmcre(y, par=c(-0.2,-1))
##' # random effect for the equilibrium probability
##' cmcre(y, par=c(-0.1,-2,-0.8))
##' # random effects for the equilibrium probability and sum of transition
##' #   intensities
##' cmcre(y, par=c(-0.1,-1.4,-0.5,-1))
##' 
##' @export cmcre
cmcre <- function(response, covariate=NULL, parameters, pcov=NULL,
	gradient=FALSE, hessian=FALSE, print.level=0, ndigit=10,
	gradtol=0.00001, steptol=0.00001, iterlim=100, fscale=1,
	typsize=abs(parameters), stepmax=parameters){
#
# set up likelihood functions
#
# ordinary Markov process
LogL1 <- function(Params){
	Params <- c(Params[1],0,Params[2],0,0)
	x <- .C("LogLikelihood1",
		as.double(Params),
		LogLikelihood=double(1),
		Error=integer(1),
		PACKAGE="repeated")
	like <- -x$LogLikelihood
	if(gradient){
		x <- .C("ScoreVector1",
			as.double(Params),
			SVector=double(2),
			PACKAGE="repeated")
		attr(like,"gradient") <- -x$SVector}
	if(hessian){
		x <- .C("Hessian1",
			as.double(Params),
			Hessian=double(4),
			PACKAGE="repeated")
		attr(like,"hessian") <- -x$Hessian}
	like}
# random effect for the equilibrium probability without covariate
LogL2 <- function(Params) {
	Params <- c(Params, 0, 0)
	x <- .C("LogLikelihood2",
		as.double(Params),
		LogLikelihood=double(1),
		Error=integer(1),
		PACKAGE="repeated")
	like <- -x$LogLikelihood
	if(gradient){
		x <- .C("ScoreVector2",
			as.double(Params),
			SVector=double(3),
			PACKAGE="repeated")
		attr(like,"gradient") <- -x$SVector}
	if(hessian){
		x <- .C("Hessian2",
			as.double(Params),
			Hessian=double(9),
			PACKAGE="repeated")
		attr(like,"hessian") <- -x$Hessian}
	like}
# random effect for the equilibrium probability with covariate
LogL3 <- function(Params) {
	Params <- if(cov) c(Params, 0) else c(Params[1],0,Params[2:3],0)
	x <- .C("LogLikelihood3",
		as.double(Params),
		LogLikelihood=double(1),
		Error=integer(1),
		PACKAGE="repeated")
	like <- -x$LogLikelihood
	if(gradient){
		x <- .C("ScoreVector3",
			as.double( Params ),
			SVector=double(np),
			cov=as.integer(cov),
			PACKAGE="repeated")
		attr(like,"gradient") <- -x$SVector}
	if(hessian){
		x <- .C("Hessian3",
			as.double(Params),
			Hessian=double(np*np),
			cov=as.integer(cov),
			PACKAGE="repeated")
		attr(like,"hessian") <- -x$Hessian}
	like}
# random effects for the equilibrium probability and sum of transition
#   intensities
LogL4 <- function(Params) {
	if(!cov) Params <- c(Params[1],0,Params[2:4])
	x <- .C("LogLikelihood4",
		as.double(Params),
		LogLikelihood=double(1),
		Error=integer(1),
		PACKAGE="repeated")
	like <- -x$LogLikelihood
	if(gradient){
		x <- .C("ScoreVector4",
			as.double(Params),
			SVector=double(np),
			cov=as.integer(cov),
			PACKAGE="repeated")
		attr(like,"gradient") <- -x$SVector}
	if(hessian){
		x <- .C("Hessian4",
			as.double(Params),
			Hessian=double(np*np),
			cov=as.integer(cov),
			PACKAGE="repeated")
		attr(like,"hessian") <-  -x$Hessian}
	like}
call <- sys.call()
if(!gradient)hessian <- FALSE
if(!is.matrix(response)||dim(response)[2]!=6)
	stop("response must be a six column matrix with subject id, time interval, and four transition frequencies")
cov <- !is.null(covariate)
if(cov&&(!is.vector(covariate)||length(covariate)!=dim(response)[1]))
	stop("covariate must be a vector with length equal to the number of rows of response")
if(cov){
	if(is.null(pcov)||length(pcov)!=1)stop("pcov must contain one initial estimate")
	parameters <- c(parameters[1],pcov,parameters[2:length(parameters)])}
#
# store data in memory
#
x <- .C("LoadData",
	response=as.double(t(cbind(response[,1],1,response[,2:6],covariate))),
	nrow=as.integer(dim(response)[1]),
	nSize=as.integer(dim(response)[2]+1+cov),
	NumSubjects=integer(1),
	Error=integer(1),
	PACKAGE="repeated")
if(x$Error>0)stop("error in storing data")
#
# define constants
#
nind <- x$NumSubjects
nobs <- sum(response[,3:6])
np <- length(parameters)
#
# choose model
#
if(np==2)mdl <- LogL1
else if(np==3)mdl <- if(cov)LogL2 else LogL3
else if(np==4)mdl <- if(cov)LogL3 else LogL4
else if(np==5&&cov)mdl <- LogL4
else stop("incorrect number of parameters")
#
# optimize with nlm
#
z0 <- nlm(mdl, parameters, hessian=TRUE, print.level=print.level,
	typsize=typsize, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
like <- z0$minimum
#
# calculate se's
#
a <- if(any(is.na(z0$hessian))||any(abs(z0$hessian)==Inf))0
	else qr(z0$hessian)$rank
if(a==np)covar <- solve(z0$hessian)
else covar <- matrix(NA,ncol=np,nrow=np)
se <- sqrt(diag(covar))
corr <- covar/(se%o%se)
dimnames(corr) <- list(1:np,1:np)
#
# remove data from memory and return results
#
.C("PurgeSubjectData", PACKAGE="repeated")
z <- list(
	call=call,
	covariate=cov,
	maxlike=like,
	aic=like+np,
	df=nobs-np,
	np=np,
	nobs=nobs,
	nind=nind,
	coefficients=z0$estimate,
	se=se,
	cov=covar,
	corr=corr,
	grad=z0$gradient,
	iterations=z0$iterations,
	code=z0$code,
	PACKAGE="repeated")
class(z) <- c("cmcre")
return(z)}

### print method
###
print.cmcre <- function(z, digits = max(3, .Options$digits - 3)){
cat("\nTwo-state Markov chain")
if(length(z$coef)-z$covariate>2)cat(" with random effect\n")
else cat("\n")
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
cat("Number of subjects    ",z$nind,"\n")
cat("Number of observations",z$nobs,"\n")
cat("\n-Log likelihood   ",z$maxlike,"\n")
cat("Degrees of freedom",z$df,"\n")
cat("AIC               ",z$aic,"\n")
cat("Iterations        ",z$iterations,"\n\n")
cat("Parameters estimates\n")
prob <- c(exp(-exp(z$coef[1])),exp(z$coef[2+z$covariate]))
int <- c(prob[1]*prob[2],(1-prob[1])*prob[2])
cname <- "beta1"
if(z$covariate)cname <- c(cname,"covariate")
cname <- c(cname,"beta2")
if(length(z$coef)-z$covariate>2)cname <- c(cname,"log(tau1)")
if(length(z$coef)-z$covariate>3)cname <- c(cname,"log(tau2)")
coef.table <- cbind(z$coef,z$se)
colname <- c("estimate","se")
dimnames(coef.table) <- list(cname,colname)
print.default(coef.table, digits=digits, print.gap=2)
cat("\nEquilibrium probability:",prob[1],"in state 1\n")
cat("\nTransition intensities:\n","1->2:",int[1],"  2->1:",int[2],"\n")
if(length(z$coef)-z$covariate>2)
	cat("\nRandom effect variance for equilibrium probability:",exp(z$coef[3+z$covariate]),"\n")
if(length(z$coef)-z$covariate>3)
	cat("\nRandom effect variance for sum of transition intensities:",exp(z$coef[4+z$covariate]),"\n")
cat("\nCorrelation matrix\n")
print.default(z$corr, digits=digits)}
