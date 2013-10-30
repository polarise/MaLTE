\name{tune}
\alias{tune}

\title{
Obtain training parameters tuned to optimise out-of-bag correlations.
}

\description{
Takes a TT.Ready.Gene object and optional parameters and returns the best values of \code{mtry}, \code{ntree}, and \code{min.probes} to use for training and prediction.
}

\usage{
tune(object, OOB=FALSE, quantreg=FALSE, verbose=FALSE )
}

\arguments{
  \item{object}{an object of class \code{TT.Ready.Gene}}
  \item{OOB}{logical; if \code{TRUE} then set the \code{TT.Params} \code{OOB} slot to \code{TRUE}}
  \item{quantreg}{logical; if \code{TRUE} then will set the \code{TT.Params} \code{quantreg} slot to \code{TRUE}}
  \item{verbose}{logical; if \code{TRUE} then will print out tuning information; should not be set when running in parallel}
}

\value{
Returns an \code{TT.Params} object of class that will be used for the current gene only.
}

\author{
Paul K. Korir
}

\seealso{
\code{link{array2seq.oob}}
}

\examples{
\dontrun{
# location of example input files
samples.fn <- paste( system.file( package="MaLTE" ), "data", "samples.txt.gz", sep="/" )
ma.fn <- paste( system.file( package="MaLTE" ), "data", "ma_data.txt.gz", sep="/" )
hts.fn <- paste( system.file( package="MaLTE" ), "data", "hts_data.txt.gz", sep="/" )
g2p.fn <- paste( system.file( package="MaLTE" ), "data", "gene_probesets.txt.gz", sep="/" )

# prepare training and test data
prepare.data( samples.fn=samples.fn, ma.fn=ma.fn, hts.fn=hts.fn, g2p.fn=g2p.fn )

# read in the training and test data
tt.ready <- read.data( train.fn="train_data.txt.gz", test.fn="test_data.txt.gz" )

# get tuned parameters for each gene (ad hoc)
tt.gene.params <- lapply( tt.ready, tune )
}
}

\keyword{predict }
