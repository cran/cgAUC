\name{cgAUC-package}
\alias{cgAUC-package}
\docType{package}
\title{
Calculate AUC when gold standard is continuous with large variables.}
\description{
In this package, the cgAUC is used to calculate the AUC-type measure raised in Obuchowski(2006) when gold standard is continuous.}
\details{
\tabular{ll}{
Package: \tab cgAUC\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-12-21\cr
License: \tab What license is it under?\cr
}
}
\author{
Yu-chia Chang
Maintainer: Yu-chia Chang <curare7177@gmail.com>
}
\references{
Chang, YCI. Maximizing an ROC type measure via linear combination of markers when the gold reference is continuous. Statistics in Medicine 2012.\cr
Obuchowski NA. An ROC-type measure of diagnostic accuracy when the gold standard is continuousscale. Statistics in Medicine 2006; 25:481 493.\cr
Obuchowski N. Estimating and comparing diagnostic tests accuracy when the gold standard is not binary. Statistics in Medicine 2005; 20:3261 3278.\cr
Friedman JH, Popescu BE. Gradient directed regularization for linear regression and classification. Technical Report, Department of Statistics, Stanford University, 2004.\cr
}
\examples{
# n = 30; p = 3;
# r.x = matrix(rnorm(n * p), , p) # raw data
# r.z = r.x[ ,1] + rnorm(n) # gold standard
# x = scale(r.x) # standardized of raw data
# z = scale(r.z) # standardized of gold standard
# h = n^(-1 / 2)
# t1 = cgAUC(r.x, r.z, h, delta = 1, auto = FALSE, tau = 1) # the delta be constant
# t1
# t2 = cgAUC(r.x, r.z, h, delta = 1, auto = TRUE, tau = 1) # the delta be variable
# t2
}
