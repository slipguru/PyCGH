\name{FLLat}

\alias{FLLat}

\title{Fused Lasso Latent Feature Model}

\description{Fits the Fused Lasso Latent Feature (FLLat) model for given
  values of \var{J} (the number of features), and \eqn{\lambda_1} and
  \eqn{\lambda_2} (the two fused lasso tuning parameters).}

\usage{
FLLat(Y, J=min(15,floor(ncol(Y)/2)), B="pc", lam1, lam2,
      thresh=10^(-4), maxiter=100, maxiter.B=1, maxiter.T=1)
}

\arguments{
  \item{Y}{A matrix of data from an aCGH experiment (usually in the form
  of log intensity ratios) or some other type of copy number data.  Rows
  correspond to the probes and columns correspond to the samples.}
  \item{J}{The number of features in the FLLat model.  The default is
  the minimum of 15 and the number of samples divided by 2.}
  \item{B}{The initial values for the features.  Can be one of
  \code{"pc"} (the first \code{J} principal components of
  \code{Y}), \code{"rand"} (a random selection of \code{J}
  columns of \code{Y}), or a user specified matrix of initial values,
  where rows correspond to the probes and columns correspond to the
  features.  The default is \code{"pc"}.}
  \item{lam1}{\eqn{\lambda_1}, the tuning parameter in the fused lasso
  penalty that controls the level of sparsity in the features.}
  \item{lam2}{\eqn{\lambda_2}, the tuning parameter in the fused lasso
  penalty that controls the level of smoothness in the features.}
  \item{thresh}{The threshold for determining when the solutions have
  converged.  The default is \eqn{10^{-4}}{10^(-4)}.}
  \item{maxiter}{The maximum number of outer iterations for the
  algorithm.  The default is 100.}
  \item{maxiter.B}{The maximum number of iterations for the inner loop
  of the algorithm for estimating the features \eqn{B}.  The default is
  1.  Increasing this may decrease the number of outer iterations but
  may still increase total run time.}
  \item{maxiter.T}{The maximum number of iterations for the inner loop
  of the algorithm for estimating the weights \eqn{\Theta}.  The default
  is 1.  Increasing this may decrease the number of outer iterations but
  may still increase total run time.}
}

\details{
  This function fits the Fused Lasso Latent Feature model to
  multi-sample aCGH data, as described in Nowak and others (2011), for
  given values of \eqn{J}, \eqn{\lambda_1} and \eqn{\lambda_2}.  Given
  aCGH data consisting of \eqn{S} samples and \eqn{L} probes, the model
  is given by: \deqn{Y=B\Theta,} where \eqn{Y} is an \eqn{L} by \eqn{S}
  matrix denoting the aCGH data (with samples in columns), \eqn{B} is a
  \eqn{L} by \eqn{J} matrix denoting the features (with features in
  columns), and \eqn{\Theta} is a \eqn{J} by \eqn{S} matrix denoting the
  weights.  Each feature describes a pattern of copy number variation
  and the weights describe the composition of each sample.
  Specifically, each sample (or column of \eqn{Y}) is modeled as a
  weighted sum of the features (or columns of \eqn{B}), with the weights
  given by the corresponding column of \eqn{\Theta}.
  
  The model is fitted by minimizing a penalized version of the residual
  sum of squares: \deqn{RSS + \sum_{j=1}^J PEN_j}{RSS + sum_{j = 1
  to J} PEN_j} where the penalty is given by: \deqn{PEN_j =
  \lambda_1\sum_{l=1}^L\left|\beta_{lj}\right| +
  \lambda_2\sum_{l=2}^L\left|\beta_{lj} -
  \beta_{l-1,j}\right|.}{PEN_j = \lambda_1*sum_{l = 1 to L}
  abs[\beta_{l,j}] + \lambda_2*sum_{l = 1 to L-1} abs[\beta_{l,j} -
  \beta_{l-1,j}].}  Here \eqn{\beta_{lj}}{\beta_{l,j}} denotes the
  \eqn{(l,j)^{th}}{(l,j)th} element of \eqn{B}.  We also constrain the \eqn{L_2}
  norm of each row of \eqn{\Theta} to be less than or equal to 1.
  
  For more details, please see Nowak and others (2011) and the package
  vignette.
}

\value{An object of class "FLLat" with components:
  \item{Beta}{The estimated features \eqn{B}.}
  \item{Theta}{The estimated weights \eqn{\Theta}.}
  \item{niter}{The number of outer iterations taken by the algorithm.}
  \item{rss}{The residual some of squares of the fitted model.}
  \item{bic}{The BIC for the fitted model.  See \code{\link{FLLat.BIC}}
    for more details.}
  There is a \code{\link{plot}} method for \code{FLLat} objects.
}

\references{G. Nowak, T. Hastie, J. R. Pollack and R. Tibshirani.  A
  Fused Lasso Latent Feature Model for Analyzing Multi-Sample aCGH
  Data. \emph{Biostatistics}, 2011, doi: 10.1093/biostatistics/kxr012}

\author{Gen Nowak \email{gen.nowak@gmail.com}, Trevor Hastie, Jonathan
  R. Pollack, Robert Tibshirani and Nicholas Johnson.}

\seealso{
  \code{\link{FLLat.BIC}}, \code{\link{FLLat.PVE}}, \code{\link{plot.FLLat}}
}

\examples{
## Load simulated aCGH data.
data(simaCGH)

## Run FLLat for J = 5, lam1 = 1 and lam2 = 9.
result <- FLLat(simaCGH,J=5,lam1=1,lam2=9)

## Plot the estimated features.
plot(result)

## Plot a heatmap of the estimated weights.
plot(result,what="weights")
}

\keyword{models}
\keyword{multivariate}
\keyword{lasso}
