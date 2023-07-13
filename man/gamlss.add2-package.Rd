\name{gamlss.add2-package}
\alias{gamlss.add2-package}
\alias{gamlss.add2}
\docType{package}
\title{More Additve Terms for GAMLSS models}
\description{
\packageDescription{gamlss.add2}
}
\details{This package is for extra additive terms within a GAMLSS model. At the moment it contains only the interface with the \pkg{earth} package.}
\author{Mikis Stasinopoulos, Vlasios Voudouris and Daniil Kiose

Maintainer: \packageMaintainer{gamlss.add2}
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R. Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M. Rigby R.A. Heller G. Voudouris V. De Bastiani F. (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}. March 30, 2017 Forthcoming by Chapman and Hall/CRC, see \url{http://www.gamlss.org/}) for an older version.

}

\keyword{package}
\seealso{\code{\link[gamlss.add]{ga}}} 
\examples{
m1 <- gamlss(R~ma(~A+Fl), data=rent, family=GA)
plot(getSmo(m1))
}
