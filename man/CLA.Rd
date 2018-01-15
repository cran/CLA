\name{CLA}
\alias{CLA}
\title{Critical Line Algorithm}
\description{
  The Critical Line Algorithm was first proposed by Markowitz to solve
  the mean-variance optimal portfolio problem.

  ... TODO: write more; give references (below) ...
}
\usage{
CLA(mu, covar, lB, uB, tol.lambda = 1e-07, give.MS = TRUE, keep.names = TRUE)
}
\arguments{
  \item{mu}{numeric vector of length \code{n} containing the expected
    return \eqn{E[R_i]} for \eqn{1=1,2,\dots,n}.}
  \item{covar}{the \eqn{n \times n}{n x n} covariance matrix of the
    returns, must be positive definite.}
  \item{lB, uB}{vectors of length \code{n} with lower and upper bounds
    for the asset weights.}
  \item{tol.lambda}{the tolerance when checking for lambda changes or
    being zero.}
  \item{give.MS}{\code{\link{logical}} indicating if \code{\link{MS}()}
    should be computed (and returned) as well.}
  \item{keep.names}{\code{\link{logical}} indicating if the
    \code{weights_set} matrix should keep the (asset) \code{names(mu)}.}
}
%% \details{
%% }
\value{
  a \code{\link{list}} with components
  \item{weights_set}{ ...}

  \item{free_indices}{ ... }
  \item{gammas}{ ... }
  \item{lambdas}{ .... }
  \item{MS_weights}{the \eqn{\mu(W)} and \eqn{\sigma(W)} corresponding
    to the asset weights \code{weights_set}, i.e., simply the same as
    \code{MS(weights_set = weights_set, mu = mu, covar = covar)}.}
}
\references{
TODO
}
\author{Alexander Norring did the very first version. Current
  implementation: Yanhao Shi and Martin Maechler}
%% \note{
%% }
\seealso{
  \code{\link{MS}}.
}
\examples{
data(muS.sp500)
## Full data taking too much time for example
set.seed(47)
iS <- sample.int(length(muS.sp500$mu), 24)

CLsp.24 <- CLA(muS.sp500$mu[iS], muS.sp500$covar[iS, iS], lB=0, uB=1/10)
CLsp.24 # using the print() method for class "CLA"

if(require(Matrix)) { ## visualize how weights change "along turning points"
  show(image(Matrix(CLsp.24$weights_set, sparse=TRUE),
             main = "CLA(muS.sp500 <random_sample(size=24)>) $ weights_set",
             xlab = "turning point", ylab = "asset number"))
}

}
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
