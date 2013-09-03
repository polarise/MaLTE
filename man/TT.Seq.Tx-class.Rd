\name{TT.Seq.Tx-class}
\Rdversion{1.1}
\docType{class}
\alias{TT.Seq.Tx-class}
\alias{gene.id,TT.Seq.Tx-method}
\alias{tx.id,TT.Seq.Tx-method}
\alias{predictions,TT.Seq.Tx-method}
\alias{trues,TT.Seq.Tx-method}
\alias{show,TT.Seq.Tx-method}

\title{Class \code{"TT.Seq.Tx"}}

\description{
Class for objects that hold results of transcript expression predictions. The \code{train.and.predict.txs} or \code{run} functions return items of this class. The \code{filter} function also takes two such objects. The \code{get.predictions} function and \code{gene.id} and \code{predictions} methods are accessors for slots in this class.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TT.Seq.Tx", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{gene.id}:}{Object of class \code{"character"} ~~ }
    \item{\code{tx.id}:}{Objects of class \code{"character"} ~~ }
    \item{\code{no.txs}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{no.samples}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{trues}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{predictions}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{cor.S}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{cor.S.pv}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{cor.P}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{cor.P.pv}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{means}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{vars}:}{Object of class \code{"numericORlogicalORmatrix"} ~~ }
    \item{\code{OOB}:}{Object of class \code{"logical"} ~~ }
  }
}

\section{Extends}{
Class \code{"\linkS4class{TT.Seq}"}, directly.
}

\section{Methods}{
  \describe{
    \item{gene.id}{\code{signature(object = "TT.Seq.Tx")}: Accessor function for gene identifier. }
    \item{tx.id}{\code{signature(object = "TT.Seq.Tx")}: Accessor function for transcript identifier(s). }
    \item{predictions}{\code{signature(object = "TT.Seq.Tx")}: Accessor function for predicted expression. }
    \item{trues}{\code{signature(object = "TT.Seq.Tx")}: Accessor function for true expression. }
    \item{show}{\code{signature(object = "TT.Seq.Tx")}: Shows the gene and transcript(s) ID, number of transcripts present, mean and variance of predictions, number of samples and whether or not this was an OOB prediction. }
	 }
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{TT.Seq-class}}, \code{\link{TT.Seq.Gene}}, \code{\link{TT.Seq.Tx}}
}

\examples{
showClass("TT.Seq.Tx")
}

\keyword{classes}
