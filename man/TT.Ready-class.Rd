\name{TT.Ready-class}
\Rdversion{1.1}
\docType{class}
\alias{TT.Ready-class}

\title{Class \code{"TT.Ready"}}

\description{
Abstract base class inherited by TT.Ready.Gene and TT.Ready.Tx.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TT.Ready", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{gene.id}:}{Object of class \code{"character"} ~~ }
    \item{\code{no.train}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{no.test}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{no.probes}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{hts.train}:}{Object of class \code{"numericORmatrix"} ~~ }
    \item{\code{hts.test}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{probes.train}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{probes.test}:}{Object of class \code{"matrix"} ~~ }
  }
}

\section{Methods}{
No methods defined with class "TT.Ready" in the signature.
}

\author{
Paul K. Korir
}

\examples{
## Not run:
showClass("TT.Ready")
## End(Not run)
}

\keyword{classes}
