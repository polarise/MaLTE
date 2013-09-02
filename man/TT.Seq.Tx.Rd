\name{TT.Seq.Tx}
\alias{TT.Seq.Tx}

\title{
Constructor for objects of class \code{"TT.Seq.Tx"}.
}

\description{
A function that takes a list item and creates an object of class \code{"TT.Seq.Tx"}. It is used internally by the \code{run} method of \code{"TT.Ready.Tx"} objects.
}

\usage{
TT.Seq.Tx(m)
}

\arguments{
  \item{m}{A list element. See Details for a complete description.}
}

\details{
Here follows a detailed description of \code{m}:
  \enumerate{
    \item{\code{gene.id}:}{a unique gene identifier}
    \item{\code{tx.id}:}{a vector of transcript identifiers}
    \item{\code{no.txs}:}{the number of transcripts with expression estimates; the length of tx.id}
    \item{\code{no.samples}:}{the number of test samples only}
    \item{\code{trues}:}{a vector of true HTS estimates of length no.samples}
    \item{\code{predictions}:}{a vector of predicted HTS estimates of length no.samples}
    \item{\code{cor.S}:}{the Spearman correlation coefficient between trues and predictions; only valid when trues exist or when predicting OOB estimates}
    \item{\code{cor.S.pv}:}{p-value for cor.S}
    \item{\code{cor.P}:}{the Pearson correlation coefficient between trues and predictions}
    \item{\code{cor.P.pv}:}{p-value fo cor.P}
    \item{\code{means}:}{mean of predicted expression estimates}
    \item{\code{vars}:}{variance of predicted expression estimates}
    \item{\code{OOB}:}{logical; whether or not predictions are out-of-bag (OOB)}
  }
}

\value{
Returns an object of class \code{TT.Seq.Tx}.
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{TT.Seq-class}}, \code{\link{TT.Seq.Gene}}
}

\keyword{transcripts}
\keyword{predictions}
