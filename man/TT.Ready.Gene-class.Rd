\name{TT.Ready.Gene-class}
\Rdversion{1.1}
\docType{class}
\alias{TT.Ready.Gene-class}
\alias{oob.run,TT.Ready.Gene-method}
\alias{run,TT.Ready.Gene-method}
\alias{show.hts.test,TT.Ready.Gene-method}
\alias{show.hts.train,TT.Ready.Gene-method}
\alias{show,TT.Ready.Gene-method}
\alias{show.probes.test,TT.Ready.Gene-method}
\alias{show.probes.train,TT.Ready.Gene-method}

\title{Class \code{"TT.Ready.Gene"}}

\description{
Each instance of a TT.Ready.Gene object holds the training and test data for a single gene. Objects of this class are used by the run method to predict gene expression.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TT.Ready.Gene", ...)}.

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

\section{Extends}{
Class \code{"\linkS4class{TT.Ready}"}, directly.
}

\section{Methods}{
  \describe{
    \item{oob.run}{\code{signature(object = "TT.Ready.Gene")}: Performs OOB prediction. }
    \item{run}{\code{signature(object = "TT.Ready.Gene")}: Performs training and prediction. }
    \item{show.hts.test}{\code{signature(object = "TT.Ready.Gene")}: Accessor function for HTS test data. }
    \item{show.hts.train}{\code{signature(object = "TT.Ready.Gene")}: Accessor function for HTS training data. }
    \item{show}{\code{signature(object = "TT.Ready.Gene")}: Shows the gene ID, number of training and test samples, number of probes for the gene and whether or not the HTS (RNA-Seq) data for test samples is available for comparison to MaLTE predictions. }
    \item{show.probes.test}{\code{signature(object = "TT.Ready.Gene")}: Accessor function for microarray probes test data. }
    \item{show.probes.train}{\code{signature(object = "TT.Ready.Gene")}: Accessor function for microarray probes training data. }
	 }
}

\author{
Paul K. Korir
}

\examples{
## Not run:
showClass("TT.Ready.Gene")
## End(Not run)
}

\keyword{classes}
