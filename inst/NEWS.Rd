% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/CLA/inst/NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{CLA}}% MM: look into ../svn-log-from.all
\newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in version 0.95-0 (2019-03-13)}{
  \subsection{New Features}{
    \itemize{
      \item \emph{Not} back-compatible \emph{changed} result: The very
      first \dQuote{turning point} from the algorithm with (\code{lambda
	= NA} etc) is no longer stored as part of the result, compatibly
      with the litterature and other implementations.

      \item Added data set \code{muS.10ex} orginally from Markowitz and Todd.

      \item Added \file{NEWS} file.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item needed \code{, drop=FALSE} in internal \code{getMatrices()},
      e.g., for the \code{*10ex} example.

      \item Test works with old R 3.3.x.

      \item In some cases such as \code{muS.10ex}, the last turning
      point with \code{lambda=0} was not computed previously.

      \item adapt \file{man/findSig.Rd} example to R 3.6.0's new default
      \code{sample()} kind.

      \item Run the \file{tests/*.R} with R versions that use one
      OpenBLAS, ATLAS, or MKL (and skip some tests on these).
    }
  }
}

\section{Changes in version 0.90-1 (2018-02-05)}{
  \subsection{New Features}{
    \itemize{
      \item Add nice \code{plot()} method base on Shi's \code{MS_plot()}
      \item Literature references in \code{?CLA}
      \item Vectorized \code{findSig()} and \code{findMu()}.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{findSig.1()} now works.
      \item test also works in \R configured with \command{--no-long-double}.
    }
  }
}

\section{Changes in version 0.90-0 (2018-01-15)}{
  \subsection{New Features}{
    \itemize{
      \item First release to CRAN
    }
  }
}

