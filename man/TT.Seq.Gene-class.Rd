\name{TT.Seq.Gene-class}
\Rdversion{1.1}
\docType{class}
\alias{TT.Seq.Gene-class}
\alias{gene.id,TT.Seq.Gene-method}
\alias{predictions,TT.Seq.Gene-method}
\alias{trues,TT.Seq.Gene-method}
\alias{show,TT.Seq.Gene-method}
\alias{cor.P,TT.Seq.Gene-method}
\alias{cor.S,TT.Seq.Gene-method}

\title{Class \code{"TT.Seq.Gene"}}

\description{
Class for objects that hold results of gene expression prediction. The train.and.predict or run functions return items of this class. The filter function also takes two such objects. The get.predictions function and gene.id and predictions methods take such objects.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TT.Seq.Gene", ...)}.
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

\section{Extends}{
Class \code{"\linkS4class{TT.Seq}"}, directly.
}

\section{Methods}{
  \describe{
    \item{gene.id}{\code{signature(object = "TT.Seq.Gene")}: Accessor function 
    for gene identifier. }
    \item{predictions}{\code{signature(object = "TT.Seq.Gene")}: Accessor 
    function for predicted expression. }
    \item{trues}{\code{signature(object = "TT.Seq.Gene")}: Accessor function for
     true expression. }
    \item{show}{\code{signature(object = "TT.Seq.Gene")}: Shows the gene ID, 
    mean and variance of predictions, number of samples and whether or not this 
    was an OOB prediction.}
    \item{cor.P}{\code{signature(object = "TT.Seq.Gene")}: Accessor function for
    Pearson correlation coefficient.}
    \item{cor.S}{\code{signature(object = "TT.Seq.Gene")}: Accessor function for
    Spearman correlation coefficient.}
	 }
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{TT.Seq-class}}, \code{\link{TT.Seq.Gene}}, \code{\link{TT.Seq.Tx}}
}

\examples{
showClass("TT.Seq.Gene")
}

\keyword{classes}
