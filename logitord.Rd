% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logitord.r
\name{logitord}
\alias{logitord}
\title{Ordinal Random Effects Models with Dropouts}
\usage{
logitord(
  y,
  id,
  out.ccov = NULL,
  drop.ccov = NULL,
  tvcov = NULL,
  out.tvcov = !is.null(tvcov),
  drop.tvcov = !is.null(tvcov),
  pout,
  pdrop,
  prand.out,
  prand.drop,
  random.out.int = TRUE,
  random.out.slope = !is.null(tvcov),
  random.drop.int = TRUE,
  random.drop.slope = !is.null(tvcov),
  binom.mix = 5,
  fcalls = 900,
  eps = 1e-04,
  print.level = 0
)
}
\arguments{
\item{y}{A vector of binary or ordinal responses with levels 1 to k and 0
indicating drop-out.}

\item{id}{Identification number for each individual.}

\item{out.ccov}{A vector, matrix, or model formula of time-constant
covariates for the outcome regression, with variables having the same
length as \code{y}.}

\item{drop.ccov}{A vector, matrix, or model formula of time-constant
covariates for the drop-out regression, with variables having the same
length as \code{y}.}

\item{tvcov}{One time-varying covariate vector.}

\item{out.tvcov}{Include the time-varying covariate in the outcome
regression.}

\item{drop.tvcov}{Include the time-varying covariate in the drop-out
regression.}

\item{pout}{Initial estimates of the outcome regression coefficients, with
length equal to the number of levels of the response plus the number of
covariates minus one.}

\item{pdrop}{Initial estimates of the drop-out regression coefficients,
with length equal to one plus the number of covariates.}

\item{prand.out}{Optional initial estimates of the outcome random
parameters.}

\item{prand.drop}{Optional initial estimates of the drop-out random
parameters.}

\item{random.out.int}{If TRUE, the outcome intercept is random.}

\item{random.out.slope}{If TRUE, the slope of the time-varying covariate is
random for the outcome regression (only possible if a time-varying
covariate is supplied and if out.tvcov and random.out.int are TRUE).}

\item{random.drop.int}{If TRUE, the drop-out intercept is random.}

\item{random.drop.slope}{If TRUE, the slope of the time-varying covariate
is random for the drop-out regression (only possible if a time-varying
covariate is supplied and if drop.tvcov and random.drop.int are TRUE).}

\item{binom.mix}{The total in the binomial distribution used to approximate
the normal mixing distribution.}

\item{fcalls}{Number of function calls allowed.}

\item{eps}{Convergence criterion.}

\item{print.level}{If 1, the iterations are printed out.}
}
\value{
A list of class \code{logitord} is returned.
}
\description{
\code{logitord} fits an longitudinal proportional odds model in discrete
time to the ordinal outcomes and a logistic model to the probability of
dropping out using a common random effect for the two.
}
\examples{


## y <- trunc(runif(20,max=4))
y <- c(1,1,1,1,2,
       1,1,1,2,3,
       1,1,2,2,2,
       1,1,3,0,0,
       1,2,3,3,0,
       0,1,1,3,0,
       3,0,0,0,0)
id <- gl(7,5)
age <- rpois(20,20)
age <- c(seq(16,16+5-1),
seq(16,16+5-1),
seq(36,36+5-1),
seq(16,16+5-1),
seq(86,86+5-1),
seq(16,16+5-1),
seq(96,96+5-1)
)
age.std <- (age - mean(age)) / sd(age)
times <- rep(1:5,7)
cbind(id,times,age,y)
logitord(y, 
         id=id, 
         out.ccov=~age.std, 
         drop.ccov=age.std, 
         pout=c(1,0,0),
         pdrop=c(1,0), 
         random.out.int=FALSE, 
         random.drop.int=FALSE,
         fcalls=40,
         eps=1)
\dontrun{         
logitord(y, 
         id, 
         tvcov=times, 
         pout=c(1,0,0), 
         pdrop=c(1,0), 
         random.out.int=FALSE, 
         random.drop.int=FALSE,
         fcalls=40,
         eps=1)
}

}
\references{
Ten Have, T.R., Kunselman, A.R., Pulkstenis, E.P. and Landis,
J.R. (1998) Biometrics 54, 367-383, for the binary case.
}
\author{
T.R. Ten Have and J.K. Lindsey
}
\keyword{models}
