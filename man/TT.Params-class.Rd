\name{TT.Params-class}
\Rdversion{1.1}
\docType{class}
\alias{TT.Params-class}
\alias{show,TT.Params-method}

\title{Class \code{"TT.Params"}}

\description{
TT.Params objects store the parameters needed for training and prediction. These are used by the run method for objects of class TT.Ready.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TT.Params", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{mtry}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{ntree}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{feature.select}:}{Object of class \code{"logical"} ~~ }
    \item{\code{min.probes}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cor.thresh}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{OOB}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "TT.Params")}: ... }
	 }
}

\author{
Paul K. Korir
}

\examples{
showClass("TT.Params")
}

\keyword{classes}
