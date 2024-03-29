\name{SOP}
\alias{SOP}
\alias{so}
\alias{gamlss.so}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{P-Spline Regression Models with GAMLSS}
\description{
The function \code{SOP()} (or \code{so()})  is an interface for fitting an additive P-Spline regression models using the separation of overlapping precision matrices (SOP) method of Rodriguez-Alvarez \emph{at al.} (2015, 2019) (implemented in the  package \pkg{SOP}) within GAMLSS.

The idea of conection  SOP with GAMLSS. is that other parameters of the distribution rather than  only for \code{mu} can be also modelled as smooth functions of the explaratory variables.  
Look in the help file of the package \pkg{SOP} for the available functions
\pkg{SOP} contains two usefull functions: tensor product P-spline smoothing i.e. \code{~f(x1,x2)} and adaptive smoothing \code{~ad(x1)}.  

The function \code{gamlss.so()} is support for \code{so()}. It is not intended to be called directly by users. The function \code{gamlss.so()} is calling the function \code{sop()} of the package \pkg{SOP}.

}
\usage{
SOP(formula, control = sop.control(...), ...)

so(formula, control = sop.control(...), ...)

gamlss.so(x, y, w, xeval = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A SOP formula i.e. \code{~f(x1)+ad(x2)+f(x3,x4)}} 
  \item{control}{the control function \code{sop.control()} from the package \pkg{SOP} }
  \item{x}{passing the x variable information to backfitting}
  \item{y}{iterative response variable}
  \item{w}{iterative weights}
  \item{xeval}{used in prediction}
  \item{\dots}{extra arguments}
}

\value{
Returms fitted models for each parameters which can be retreived using the function \code{gerSmo()}.
}
\references{
Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines and penalties. \emph{Statistical Science}, \bold{11 (2)}, 89--121.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. \doi{10.1201/9780429298547}
	An older version can be found in \url{https://www.gamlss.com/}.
	
Rodriguez-Alvarez, M.X., Lee, D. J., Kneib, T., Durban, M., and Eilers, P. (2015). Fast smoothing parameter separation in multidimensional generalized P-splines: the SAP algorithm. \emph{Statistics and Computing}, \bold{25 (5)}, 941--957.

Rodriguez-Alvarez, M.X., Durban, M., Lee, D. J. and Eilers, P. (2019). On the estimation of variance parameters in non-standard generalised linear mixed models: application to penalised smoothing. \emph{Statistics and Computing}, \bold{29 (3)}, 483--500.	
	
	
Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M. Rigby R.A. Heller G. Voudouris V. De Bastiani F. (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}. March 30, 2017 Forthcoming by Chapman and Hall/CRC, see \url{http://www.gamlss.org/}) for an older version.

Stasinopoulos, M. D., Rigby, R. A., and De Bastiani F., (2018) GAMLSS: a distributional regression approach, \emph{Statistical Modelling}, Vol. \bold{18}, 
pp, 248-273, SAGE Publications Sage India: New Delhi, India. \doi{10.1177/1471082X18759144}


(see also  \url{https://www.gamlss.com/}).
}

\author{Mikis Stasinopoulos with the help of Xose Rodriguez-Alvarez and M Lee }


\seealso{\code{\link[gamlss.add]{ga}}}

\examples{
m1 <- gamlss(R~SOP(~f(Fl)+f(A)), data=rent)
m2 <- gamlss(R~SOP(~f(Fl)+ad(A)), data=rent)
m3 <- gamlss(R~SOP(~f(Fl, A)), data=rent)
par(mfrow=c(2,2))
plot(getSmo(m1), ask=FALSE)
plot(getSmo(m2), ask=FALSE)
par(mfrow=c(1,1))
plot(getSmo(m3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{regression}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
