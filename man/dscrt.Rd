\name{dscrt}
\alias{dscrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
dscrt}
\description{
discrete function, when variable is discrete.}
\usage{
dscrt(y, z, l)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
The potential variables. It is a matrix with column of values of a variables. It should be standardized in this application.}
  \item{z}{
The gold standard variable. It should be standardized.}
  \item{l}{
Linear combination. A vector.}
}
\details{
Discrete function, when variable is discrete.}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{theta.h.p}{The estimate of theta when variable is discrete.}
\item{var}{The variance of estimate of theta.}
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
function (y, z, l) 
{
    n = length(z)
    psii = rep(0, n)
    temp = 0
    y = y \%*\% l
    for (i in 1:n) {
        psij = 0
        for (j in 1:n) {
            if (i != j) {
                temp = ifelse(((y[i] - y[j]) * (z[i] - z[j])) > 
                  0, 1, ifelse(((y[i] - y[j]) == 0 || (z[i] - 
                  z[j]) == 0), 0.5, 0))
                psij = psij + temp
            }
        }
        psii[i] = psij
    }
    theta.h.p = sum(psii)/(n * (n - 1))
    var = sum(((psii/(n - 1)) - theta.h.p)^2)/((n/2) * ((n/2) - 
        1))
    return(list(theta.h.p = theta.h.p, var = var))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dscrt }
