\name{cmcre}
\title{Continuous-time Two-state Markov Processes with Random Effect}
\alias{cmcre}
\usage{
cmcre(response, covariate=NULL, parameters, pcov=NULL, gradient=FALSE,
	hessian=FALSE, print.level=0, ndigit=10, gradtol=0.00001,
	steptol=0.00001, iterlim=100, fscale=1, typsiz=abs(parameters),
	stepmax=parameters)
}
\description{
\code{cmcre} fits a two-state Markov process in continuous time,
possibly with one or two random effects and/or one covariate.
}
\arguments{
\item{response}{A six-column matrix. Column 1: subject identification
(subjects can occupy several rows); column 2: time gap between
events; columns 3-6: transition matrix frequencies.}
\item{covariate}{An optional vector of length equal to the number of
rows of \code{response} upon which the equilibrium probability may
depend.}
\item{parameters}{Initial parameter estimates. The number of them
determines the model fitted (minimum 2, yielding an ordinary Markov
process). 1: beta1=log(-log(equilibrium probability)); 2: beta2=log(sum of
transition intensities); 3: log(tau1)=log(random effect variance for
equilibrium probability); 4: log(tau2)=log(random effect variance for sum
of transition intensities).}
\item{pcov}{Initial parameter estimate for the covariate influencing
the equilibrium probability: exp(-exp(beta1+beta*covariate)).}
\item{gradient}{If TRUE, analytic gradient is used (with
accompanying loss of speed).}
\item{hessian}{If TRUE, analytic hessian is used (with
accompanying loss of speed).}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\references{
Cook, R.J. (1999) A mixed model for two-state Markov processes under
panel observations. Biometrics 55, 915-920.
}
\value{
A list of class \code{cmcre} is returned.
}
\author{R.J. Cook and J.K. Lindsey}
\seealso{
\code{\link[repeated]{chidden}}, \code{\link[repeated]{hidden}}.
}
\examples{
# 12 subjects observed at intervals of 7 days
y <- matrix(c(1,7,1,2,3,5,
	2,7,10,2,2,0,
	3,7,7,0,1,1,
	4,7,2,1,0,7,
	5,7,1,1,1,11,
	6,7,5,4,4,1,
	7,7,1,1,1,8,
	8,7,2,3,4,2,
	9,7,9,0,0,0,
	10,7,0,1,2,8,
	11,7,8,2,2,1,
	12,7,9,2,2,1),ncol=6, byrow=TRUE)
# ordinary Markov process
cmcre(y, par=c(-0.2,-1))
# random effect for the equilibrium probability
cmcre(y, par=c(-0.1,-2,-0.8))
# random effects for the equilibrium probability and sum of transition
#   intensities
cmcre(y, par=c(-0.1,-1.4,-0.5,-1))
}
\keyword{models}
