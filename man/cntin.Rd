\name{cntin}
\alias{cntin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cntin}
\description{
Continue function, when variable was continue.}
\usage{
cntin(y, z, l, h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The potential variables. It is a matrix with column of values of a variables. It should be standardized in this application.
}
  \item{z}{
The gold standard variable. It should be standardized.}
  \item{l}{
Linear combination. A vector.}
  \item{h}{
The value of h falls into (n^(-1/2), n^(-1/5)).}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{theta.sh.h.p}{The estimate of the theta of Chang(2012).}
\item{var}{The variance of estimate of the theta of Chang(2012).}

}
\references{
Chang, YCI. Maximizing an ROC type measure via linear combination of markers when the gold reference is continuous. Statistics in Medicine 2012.\cr
Obuchowski NA. An ROC-type measure of diagnostic accuracy when the gold standard is continuousscale. Statistics in Medicine 2006; 25:481 493.\cr
Obuchowski N. Estimating and comparing diagnostic tests accuracy when the gold standard is not binary. Statistics in Medicine 2005; 20:3261 3278.\cr
Friedman JH, Popescu BE. Gradient directed regularization for linear regression and classification. Technical Report, Department of Statistics, Stanford University, 2004.\cr
}
\author{
Yu-chia Chang}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, z, l, h) 
{
    n = dim(y)[1]
    psii = rep(0, n)
    temp = 0
    for (i in 1:n) {
        psij = 0
        for (j in 1:n) {
            if (i != j) {
                temp = s.h(t(y[i, ] - y[j, ]) \%*\% l, h) * (2 * 
                  s.h((z[i] - z[j]), h) - 1) + 1 - s.h((z[i] - 
                  z[j]), h)
                psij = psij + temp
            }
        }
        psii[i] = psij
    }
    theta.sh.h.p = sum(psii)/(n * (n - 1))
    var = sum(((psii/(n - 1)) - theta.sh.h.p)^2)/((n/2) * ((n/2) - 
        1))
    return(list(theta.sh.h.p = theta.sh.h.p, var = var))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cntin }
