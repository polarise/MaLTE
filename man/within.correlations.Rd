\name{within_correlations}
\alias{within_correlations}

\title{
Compute within-sample correlations for MaLTE and another summarisation
}

\description{
This function computes within-sample correlations for each individual sample for
both MaLTE and another summarisation method. It relays the results as a list 
having slots for Pearson and Spearman correlation estimates.
}

\usage{
within_correlations(tt.seq, affy.fn, affy.name="median_polish", 
	check.names=TRUE)
}

\arguments{
  \item{tt.seq}{a list with objects of class \code{TT.Seq.Gene}}
  \item{affy.fn}{the filename of the summarisation data; see Notes below}
  \item{affy.name}{the name to be used to identify this Affymetrix 
  	summarisation}
  \item{check.names}{ensure that \code{affy.name} is a valid \code{R} variable}
}

\value{
Returns a list with two slots: Pearson and Spearman. Each slot is a data frame
with the first column having MaLTE within-sample correlations and the second 
column having those for the Affymetrix summarisation.
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
## Not run:
# location of example input files
# to be added...
## End(Not run)
}

\keyword{within}
\keyword{correlation}
