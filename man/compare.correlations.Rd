\name{compare.correlations}
\alias{compare.correlations}

\title{
Compare MaLTE to summarisation correlations
}

\description{
This function produces a data frame of correlation coefficients (Pearson and 
Spearman) of expression estimates using MaLTE and a summarisation algorithm. 
Each row of this data frame is a gene with the respective correlations. The 
columns prefixed by 'affy' denote summarisation; the last two are due to MaLTE.
}

\usage{
compare.correlations(tt.seq, affy.fn, raised=FALSE, filter.fpkm=0, 
	filter.count=0)
}

\arguments{
  \item{tt.seq}{a list with objects of class \code{TT.Seq.Gene}}
  \item{affy.fn}{the filename of the summarisation data; see Notes below}
  \item{raised}{logical; whether to expression estimates to the power of two 
  (for RMA); default is \code{FALSE}}
  \item{filter.fpkm}{the minimum FPKM to consider as expressed in RNA-Seq; used
  	together with \code{filter.count}; if there are at least \code{filter.count}
  	samples with at least \code{filter.fpkm} FPKM then consider this gene; 
  	otherwise give it \code{NA}s}
  \item{filter.count}{the minimum number of samples needed to cross the 
  \code{filter.fpkm} FPKM threshold for comparisons to be carried out}
}

\value{
Returns a data frame with the following header: "all.genes, affy.cor.S, 
affy.cor.P, cor.S, cor.P", where columns prefixed by "affy" denoting 
correlations for summarisation (median-polish, PLIER) expression estimates.
}

\author{
Paul K. Korir
}

\note{
The file represented by \code{affy.fn} must have the same sample labels in the 
header as those defined in \code{samples.txt}. The first column must also use 
the same gene identifiers as those in RNA-Seq training/test data.

This function only applies to gene expression; comparisons cannot be done for 
transcripts because summarisation to transcript isoform expression is not 
available.
}

\examples{
# location of example input files
# to be added...
}

\keyword{compare}
\keyword{correlation}
