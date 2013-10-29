\name{TT.Seq.Gene}
\alias{TT.Seq.Gene}

\title{
Constructor for objects of class \code{"TT.Seq.Gene"}.
}

\description{
A function that takes a list item and creates an object of class \code{"TT.Seq.Gene"}. It is used internally by the \code{run} method of \code{"TT.Ready.Gene"} objects.
}

\usage{
TT.Seq.Gene(m)
}

\arguments{
  \item{m}{
\code{m} is a list with the following slots: \code{gene.id}, \code{no.samples}, \code{trues}, \code{predictions}, \code{predictions.lower}, \code{predictions.upper}, \code{cor.S}, \code{cor.S.pv}, \code{cor.P}, \code{cor.P.pv}, \code{means}, \code{vars}, and \code{OOB}.
}
}

\details{
Here follows a detailed description of \code{m}:
  \enumerate{
    \item{\code{gene.id}:}{a unique gene identifier}
    \item{\code{no.samples}:}{the number of test samples only}
    \item{\code{trues}:}{a vector of true HTS estimates of length no.samples}
    \item{\code{predictions}:}{a vector of predicted HTS estimates of length no.samples}
    \item{\code{predictions.lower}:}{a vector of 5\% quantiles of predicted HTS estimates obtained only when using \code{quantregForests}}
    \item{\code{predictions.upper}:}{a vector of 95\% quantiles of predicted HTS estimates obtained only when using \code{quantregForests}}
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
Returns an object of class \code{TT.Seq.Gene}.
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{TT.Seq-class}}, \code{\link{TT.Seq.Tx}}
}

\keyword{gene}
\keyword{prediction}
