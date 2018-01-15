\name{findSig}
\alias{findSig} %% keep in sync with ./findMu.Rd
\title{Find sigma(W) and W, given mu(W) and CLA result}
\description{
  Find \eqn{\sigma(W)} and \eqn{W}, given \eqn{\mu(W}) and \code{\link{CLA}} result.
}
\usage{
findSig(Mu0, result, covar, equal.tol)
}
\arguments{
  \item{Mu0}{numeric vector}
  \item{result}{ .. }
  \item{covar}{  ... }
  \item{equal.tol}{numeric tolerance to be used in
    \code{\link{all.equal}(.., tolerance = equal.tol)}.}% TODO:  is this correct?
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
  a \code{\link{list}} with components
  \item{Sig}{.....}
  \item{weight}{.....}
}
\references{
  Master thesis, p.33
}

\seealso{
 \code{\link{findMu}},
 \code{\link{CLA}}, \code{\link{MS}}.
}
\examples{

 ##  ...

}
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
