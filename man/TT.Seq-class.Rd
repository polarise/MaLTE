\name{TT.Seq-class}
\Rdversion{1.1}
\docType{class}
\alias{TT.Seq-class}

\title{Class \code{"TT.Seq"}}

\description{
Abstract base class inherited by TT.Seq.Gene and TT.Seq.Tx.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TT.Seq", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{gene.id}:}{Object of class \code{"character"} ~~ }
    \item{\code{no.samples}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{trues}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{predictions}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{predictions.lower}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{predictions.upper}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{cor.S}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{cor.S.pv}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{cor.P}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{cor.P.pv}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{means}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{vars}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{OOB}:}{Object of class \code{"logical"} ~~ }
  }
}

\section{Methods}{
No methods defined with class "TT.Seq" in the signature.
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{TT.Seq.Gene}}, \code{\link{TT.Seq.Tx}}
}

\examples{
## Not run:
showClass("TT.Seq")
## End(Not run)
}

\keyword{classes}
