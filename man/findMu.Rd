\name{findMu}
\alias{findMu} %% keep in sync with ./findSig.Rd
%%                                    ^^^^^^^^^^ %% ??? TODO ??? one help page also for  findSig() ???
\title{Find  mu(W) and W, given sigma(W) and CLA result}
\description{
  Find \eqn{\mu(W)} and \eqn{W}, given \eqn{\sigma(W}) and \code{\link{CLA}} result.
}
\usage{
findMu(Sig0, result, covar, tol.unir = 1e-06, equal.tol = 1e-06)
}
\arguments{
  \item{Sig0}{numeric vector of \eqn{\sigma(W)} values.}
  \item{result}{ .. }
  \item{covar}{  ... }
  \item{tol.unir}{numeric tolerance passed to \code{\link{uniroot}}.}
  \item{equal.tol}{numeric tolerance to be used in
    \code{\link{all.equal}(.., tolerance = equal.tol)}.}% TODO:  is this correct?
}
\value{
  a \code{\link{list}} with components
  \item{Mu}{.....}
  \item{weight}{.....}
}
\references{
  Master thesis, p.33
}

\seealso{
 \code{\link{findSig}},
 \code{\link{CLA}}, \code{\link{MS}}.
}
\examples{
## TODO
}
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
