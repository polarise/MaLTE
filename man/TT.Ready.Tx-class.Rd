\name{TT.Ready.Tx-class}
\Rdversion{1.1}
\docType{class}
\alias{TT.Ready.Tx-class}
\alias{oob.run,TT.Ready.Tx-method}
\alias{run,TT.Ready.Tx-method}
\alias{show.hts.test,TT.Ready.Tx-method}
\alias{show.hts.train,TT.Ready.Tx-method}
\alias{show,TT.Ready.Tx-method}
\alias{show.probes.test,TT.Ready.Tx-method}
\alias{show.probes.train,TT.Ready.Tx-method}

\title{Class \code{"TT.Ready.Tx"}}

\description{
Each instance of a TT.Ready.Tx object holds the training and test data for the set of transcript isoforms of a single gene. Objects of this class are used by the run method to predict transcript isoform expression.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TT.Ready.Tx", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{tx.id}:}{Object of class \code{"character"} ~~ }
    \item{\code{no.txs}:}{Object of class \code{"numeric"} ~~ }
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
    \item{oob.run}{\code{signature(object = "TT.Ready.Tx")}: Performs OOB prediction. }
    \item{run}{\code{signature(object = "TT.Ready.Tx")}: Performs training and prediction. }
    \item{show.hts.test}{\code{signature(object = "TT.Ready.Tx")}: Accessor function for HTS test data. }
    \item{show.hts.train}{\code{signature(object = "TT.Ready.Tx")}: Accessor function for HTS training data. }
    \item{show}{\code{signature(object = "TT.Ready.Tx")}: Shows the gene and transcript(s) ID, number of transcripts, number of training and test samples, number of probes for the gene and whether or not the HTS (RNA-Seq) data for test samples is available for comparison to MaLTE predictions. }
    \item{show.probes.test}{\code{signature(object = "TT.Ready.Tx")}: Accessor function for microarray probes test data. }
    \item{show.probes.train}{\code{signature(object = "TT.Ready.Tx")}: Accessor function for microarray probes training data. }
	 }
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{TT.Ready-class}}, \code{\link{TT.Ready.Tx}}, \code{link{TT.Ready.Gene}}
}

\examples{
## Not run:
showClass("TT.Ready.Tx")
## End(Not run)
}

\keyword{classes}
