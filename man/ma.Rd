\name{ma}
\alias{ma}
\alias{ma.control}
\alias{gamlss.ma}

\title{Fitting Multivariate Adaptive Regression Splines within GAMLSS}
\description{The function \code{ma()} is an interface for fitting  Friedman's  Multivariate Adaptive Regression Splines within GAMLSS. The idea here is that MARS can be used also for other parameters of the distribution rather than  only for \code{mu}. The package for fitting MARS in \R{} is \pkg{earth}. Look in the help file for the function \code{earth()} and its package vignette `Notes on the earth package'  for more details for the MARS implementation.  The function \code{gamlss.ma()} is support for \code{ma()}. It is not intended to be called directly by users. The function \code{gamlss.ma()} is calling the \R{} function \code{earth()}.}
\usage{
ma(formula, control = ma.control(...), ...)

ma.control(wp = NULL,  pmethod = "backward", keepxy = FALSE, 
        trace = 0, degree = 1, nprune = NULL, ncross = 1, 
        nfold = 0, stratify = TRUE,  Scale.y = 1, ...)
        
gamlss.ma(x, y, w, xeval = NULL, ...)        
}
\arguments{
  \item{formula}{a  model formula}
  \item{control}{controlling arguments to \code{earth()}}
  \item{wp}{response weights see the help file for the function  \code{earth()}}
  \item{pmethod}{Pruning method, see the help file in \code{earth()}}
  \item{keepxy}{keeping y, x and weights see the help file in \code{earth()}}  
  \item{trace}{trace earth's execution see the help file in \code{earth()}}  
  \item{degree}{maximum degree of interaction (Friedman's mi). Default is 1, meaning build an additive model (i.e., no interaction terms) see  the help in \code{earth()}}  
  \item{nprune}{Maximum number of terms (including intercept) in the pruned model see \code{earth()}}  
  \item{ncross}{Only applies if nfold>1 see the help file in \code{earth()}} 
  \item{nfold}{Number of cross-validation folds see the help file in \code{earth()}} 
  \item{stratify}{Only applies if nfold>1 see the help file in \code{earth()}} 
  \item{Scale.y}{Scalinf y in the forward pass for better numeric stability see the hlp file for \code{earth()}} 
  
  \item{x}{passing the x variable information to backfitting}
  \item{y}{iterative response variable}
  \item{w}{iterative weights}
  \item{xeval}{used in prediction}
  \item{\dots}{extra arguments}
}

\value{a fitted earth object}
\references{

Friedman (1991)  Multivariate Adaptive Regression Splines (with discussion) Annuals of Statistics 19/1. 1-141 \url{https://statistics.stanford.edu/research/multivariate-adaptive-regression-splines}

Friedman (1993) Fast MARS, Stanford University Department of Statistics. Technical report 110 \url{https://statistics.stanford.edu/research/fast-mars}  

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. \doi{10.1201/9780429298547}
	An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M. Rigby R.A. Heller G. Voudouris V. De Bastiani F. (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}. March 30, 2017 Forthcoming by Chapman and Hall/CRC, see \url{http://www.gamlss.org/}) for an older version.

Stasinopoulos, M. D., Rigby, R. A., and De Bastiani F., (2018) GAMLSS: a distributional regression approach, \emph{Statistical Modelling}, Vol. \bold{18}, 
pp, 248-273, SAGE Publications Sage India: New Delhi, India. \doi{10.1177/1471082X18759144}

(see also  \url{https://www.gamlss.com/}).

}
\author{Mikis Stasinopoilos, Vlasios Voudouris and Daniil Kiose}

\seealso{\code{\link[gamlss.add]{ga}}}

\examples{
m1 <- gamlss(R~ma(~Fl+A+B+loc), data=rent, family=GA)
library(plotmo)
plotmo(getSmo(m1))
}

\keyword{regression}
