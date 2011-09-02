\name{plot.FLLat}

\alias{plot.FLLat}

\title{Plots Results from the Fused Lasso Latent Feature Model}

\description{Plots either the features or a heatmap of the weights from
  an object of class \code{FLLat}.}

\usage{
\method{plot}{FLLat}(x, what=c("features","weights"), f.mar=c(5,3,4,2),
     f.xlab="Probe", w.mar=c(3,5,0,2), samp.names=1:ncol(x$Theta),
     hc.meth="complete", \dots)
}

\arguments{
  \item{x}{An object of class \code{FLLat}, as returned by
  \code{\link{FLLat}}.}
  \item{what}{Whether to plot the features or a heatmap of the weights.
  Default is features.}
  \item{f.mar}{The margins for the plot of each feature.}
  \item{f.xlab}{The label for the x axis for the plot of each feature.}
  \item{w.mar}{The margins for the heatmap of the weights.}
  \item{samp.names}{The sample names used to label the columns in the
  heatmap of the weights.}
  \item{hc.meth}{The agglomeration method to be used in the hierarchical
  clustering of the columns of \eqn{\Theta}.  See \code{\link{hclust}}.}
  \item{\dots}{Further graphical parameters, for the \code{\link{plot}}
  function when \code{what="features"} and for the \code{\link{image}}
  function when \code{what="weights"}.}
}

\details{
  This function plots the features (\eqn{B}) or a heatmap of the weights
  (\eqn{\Theta}) from a \code{FLLat} object.  The features are plotted
  in order of decreasing total magnitude, given by
  \eqn{\sum_{l=1}^L\beta_{lj}^2}{sum_{l = 1 to L} \beta_{l,j}^2}, where
  \eqn{\beta_{lj}}{\beta_{l,j}} for \eqn{l=1,\ldots,L}{l = 1 to L} denotes the
  \eqn{j^{th}}{jth} feature (or column of \eqn{B}).  Similarly, the rows of the
  heatmap of the weights are re-ordered in the same way.  Further, the
  heatmap of the weights includes a dendrogram of a hierarchical
  clustering of the samples (or columns of \eqn{\Theta}).
  
  For more details, please see Nowak and others (2011) and the package
  vignette.
}

\references{G. Nowak, T. Hastie, J. R. Pollack and R. Tibshirani.  A
  Fused Lasso Latent Feature Model for Analyzing Multi-Sample aCGH
  Data. \emph{Biostatistics}, 2011, doi: 10.1093/biostatistics/kxr012}

\author{Gen Nowak \email{gen.nowak@gmail.com}, Trevor Hastie, Jonathan
  R. Pollack, Robert Tibshirani and Nicholas Johnson.}

\seealso{
  \code{\link{FLLat}}
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
