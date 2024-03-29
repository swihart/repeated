% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chidden.r
\name{chidden}
\alias{chidden}
\title{Continuous-time Hidden Markov Chain Models}
\usage{
chidden(
  response = NULL,
  totals = NULL,
  times = NULL,
  distribution = "Bernoulli",
  mu = NULL,
  cmu = NULL,
  tvmu = NULL,
  pgamma,
  pmu = NULL,
  pcmu = NULL,
  ptvmu = NULL,
  pshape = NULL,
  pfamily = NULL,
  par = NULL,
  pintercept = NULL,
  delta = NULL,
  envir = parent.frame(),
  print.level = 0,
  ndigit = 10,
  gradtol = 1e-05,
  steptol = 1e-05,
  fscale = 1,
  iterlim = 100,
  typsize = abs(p),
  stepmax = 10 * sqrt(p \%*\% p)
)
}
\arguments{
\item{response}{A list of two or three column matrices with counts or
category indicators, times, and possibly totals (if the distribution is
binomial), for each individual, one matrix or dataframe of counts, or an
object of class, \code{response} (created by
\code{\link[rmutil]{restovec}}) or \code{repeated} (created by
\code{\link[rmutil]{rmna}} or \code{\link[rmutil]{lvna}}). If the
\code{repeated} data object contains more than one response variable, give
that object in \code{envir} and give the name of the response variable to
be used here. If there is only one series, a vector of responses may be
supplied instead.

Multinomial and ordinal categories must be integers numbered from 0.}

\item{totals}{If response is a matrix, a corresponding matrix of totals if
the distribution is binomial. Ignored if response has class,
\code{response} or \code{repeated}.}

\item{times}{If \code{response} is a matrix, a vector of corresponding
times, when they are the same for all individuals. Ignored if response has
class, \code{response} or \code{repeated}.}

\item{distribution}{Bernoulli, Poisson, multinomial, proportional odds,
continuation ratio, binomial, exponential, beta binomial, negative
binomial, normal, inverse Gauss, logistic, gamma, Weibull, Cauchy, Laplace,
Levy, Pareto, gen(eralized) gamma, gen(eralized) logistic, Hjorth, Burr,
gen(eralized) Weibull, gen(eralized) extreme value, gen(eralized) inverse
Gauss, power exponential, skew Laplace, Student t, or (time-)discretized
Poisson process. (For definitions of distributions, see the corresponding
[dpqr]distribution help.)}

\item{mu}{A general location function with two possibilities: (1) a list of
formulae (with parameters having different names) or functions (with one
parameter vector numbering for all of them) each returning one value per
observation; or (2) a single formula or function which will be used for all
states (and all categories if multinomial) but with different parameter
values in each so that pmu must be a vector of length the number of
unknowns in the function or formula times the number of states (times the
number of categories minus one if multinomial).}

\item{cmu}{A time-constant location function with three possibilities: (1)
a list of formulae (with parameters having different names) or functions
(with one parameter vector numbering for all of them) each returning one
value per individual; (2) a single formula or function which will be used
for all states (and all categories if multinomial) but with different
parameter values in each so that pcmu must be a vector of length the number
of unknowns in the function or formula times the number of states (times
the number of categories minus one if multinomial); or (3) a function
returning an array with one row for each individual, one column for each
state of the hidden Markov chain, and, if multinomial, one layer for each
category but the last. If used, this function or formula should contain the
intercept. Ignored if \code{mu} is supplied.}

\item{tvmu}{A time-varying location function with three possibilities: (1)
a list of formulae (with parameters having different names) or functions
(with one parameter vector numbering for all of them) each returning one
value per time point; (2) a single formula or function which will be used
for all states (and all categories if multinomial) but with different
parameter values in each so that ptvmu must be a vector of length the
number of unknowns in the function or formula times the number of states
(times the number of categories minus one if multinomial); or (3) a
function returning an array with one row for each time point, one column
for each state of the hidden Markov chain, and, if multinomial, one layer
for each category but the last. This function or formula is usually a
function of time; it is the same for all individuals. It only contains the
intercept if \code{cmu} does not. Ignored if \code{mu} is supplied.}

\item{pgamma}{A square \code{mxm} matrix of initial estimates of the
continuous-time hidden Markov transition matrix, where \code{m} is the
number of hidden states. Rows can either sum to zero or the diagonal
elements can be zero, in which case they will be replaced by minus the sum
of the other values on the rows. If the matrix contains zeroes off
diagonal, these are fixed and not estimated.}

\item{pmu}{Initial estimates of the unknown parameters in \code{mu}.}

\item{pcmu}{Initial estimates of the unknown parameters in \code{cmu}.}

\item{ptvmu}{Initial estimates of the unknown parameters in \code{tvmu}.}

\item{pshape}{Initial estimate(s) of the dispersion parameter, for those
distributions having one. This can be one value or a vector with a
different value for each state.}

\item{pfamily}{Initial estimate of the family parameter, for those
distributions having one.}

\item{par}{Initial estimate of the autoregression parameter.}

\item{pintercept}{For multinomial, proportional odds, and continuation
ratio models, \code{p-2} initial estimates for intercept contrasts from the
first intercept, where \code{p} is the number of categories.}

\item{delta}{Scalar or vector giving the unit of measurement (always one
for discrete data) for each response value, set to unity by default. For
example, if a response is measured to two decimals, delta=0.01. If the
response is transformed, this must be multiplied by the Jacobian. For
example, with a log transformation, \code{delta=1/response}. Ignored if
response has class, \code{response} or \code{repeated}.}

\item{envir}{Environment in which model formulae are to be interpreted or a
data object of class, \code{repeated}, \code{tccov}, or \code{tvcov}; the
name of the response variable should be given in \code{response}. If
\code{response} has class \code{repeated}, it is used as the environment.}

\item{print.level}{Arguments for nlm.}

\item{ndigit}{Arguments for nlm.}

\item{gradtol}{Arguments for nlm.}

\item{steptol}{Arguments for nlm.}

\item{fscale}{Arguments for nlm.}

\item{iterlim}{Arguments for nlm.}

\item{typsize}{Arguments for nlm.}

\item{stepmax}{Arguments for nlm.}
}
\value{
A list of classes \code{hidden} and \code{recursive} (unless
multinomial, proportional odds, or continuation ratio) is returned that
contains all of the relevant information calculated, including error codes.
}
\description{
\code{chidden} fits a two or more state hidden Markov chain model with a
variety of distributions in continuous time. All series on different
individuals are assumed to start at the same time point. If the time points
are equal, discrete steps, use \code{\link[repeated]{hidden}}.
}
\details{
The time-discretized Poisson process is a continuous-time hidden Markov
model for Poisson processes where time is then discretized and only
presence or absence of one or more events is recorded in each, perhaps
unequally-spaced, discrete interval.

For quantitative responses, specifying \code{par} allows an `observed'
autoregression to be fitted as well as the hidden Markov chain.

All functions and formulae for the location parameter are on the
(generalized) logit scale for the Bernoulli, binomial, and multinomial
distributions. Those for intensities of the discretized Poisson process are
on the log scale.

If \code{cmu} and \code{tvmu} are used, these two mean functions are
additive so that interactions between time-constant and time-varying
variables are not possible.

The algorithm will run more quickly if the most frequently occurring time
step is scaled to be equal to unity.

The object returned can be plotted to give the probabilities of being in
each hidden state at each time point. See \code{\link[repeated]{hidden}}
for details. For distributions other than the multinomial, proportional
odds, and continuation ratio, the (recursive) predicted values can be
plotted using \code{\link[rmutil]{mprofile}} and
\code{\link[rmutil]{iprofile}}.
}
\examples{

# model for one randomly-generated binary series
y <- c(rbinom(10,1,0.1), rbinom(10,1,0.9))
mu <- function(p) array(p, c(1,2))
print(z <- chidden(y, times=1:20, dist="Bernoulli",
	pgamma=matrix(c(-0.1,0.2,0.1,-0.2),ncol=2),
	cmu=mu, pcmu=c(-2,2)))
# or equivalently
print(z <- chidden(y, times=1:20, dist="Bernoulli",
	pgamma=matrix(c(-0.1,0.2,0.1,-0.2),ncol=2),
	cmu=~1, pcmu=c(-2,2)))
# or
print(z <- chidden(y, times=1:20, dist="Bernoulli",
	pgamma=matrix(c(-0.1,0.2,0.1,-0.2),ncol=2),
	mu=~rep(a,20), pmu=c(-2,2)))
mexp(z$gamma)
par(mfcol=c(2,2))
plot(z)
plot(iprofile(z), lty=2)
plot(mprofile(z), add=TRUE)
print(z <- chidden(y, times=(1:20)*2, dist="Bernoulli",
	pgamma=matrix(c(-0.05,0.1,0.05,-0.1),ncol=2),
	cmu=~1, pcmu=c(-2,2)))
mexp(z$gamma) \%*\% mexp(z$gamma)
plot(z)
plot(iprofile(z), lty=2)
plot(mprofile(z), add=TRUE)

}
\references{
MacDonald, I.L. and Zucchini, W. (1997) Hidden Markov and other
Models for Discrete-valued Time Series. Chapman & Hall.

For time-discretized Poisson processes, see

Davison, A.C. and Ramesh, N.I. (1996) Some models for discretized series of
events. JASA 91: 601-609.
}
\author{
J.K. Lindsey
}
\keyword{models}
