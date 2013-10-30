\name{oob.run}
\alias{oob.run}

\title{
OOB prediction.
}

\description{
Performs OOB prediction given an object of class \code{TT.Ready}.
}

\usage{
oob.run(object, params.object, ...)
}

\arguments{
  \item{object}{an object of class \code{TT.Ready}}
  \item{params.object}{an object of class TT.Params that contains the parameters to be passed to the conditional random forest}
  \item{...}{other parameters passed onto \code{train.and.predict}}
}

\value{
Returns an object of class \code{TT.Seq} containing OOB predictions.
}

\author{
Paul K. Korir
}

\seealso{
\code{link{run}}
}

\examples{
## Not run:
# location of example input files
samples.fn <- paste( system.file( package="MaLTE" ), "data", "samples.txt.gz", sep="/" )
ma.fn <- paste( system.file( package="MaLTE" ), "data", "ma_data.txt.gz", sep="/" )
hts.fn <- paste( system.file( package="MaLTE" ), "data", "hts_data.txt.gz", sep="/" )
g2p.fn <- paste( system.file( package="MaLTE" ), "data", "gene_probesets.txt.gz", sep="/" )

# prepare training and test data
prepare.data( samples.fn=samples.fn, ma.fn=ma.fn, hts.fn=hts.fn, g2p.fn=g2p.fn )

# read in training and test data
tt.ready <- read.data( train.fn="train_data.txt.gz", test.fn="test_data.txt.gz" )

# set params
tt.params <- TT.Params( )

# OOB run
prediction <- oob.run( tt.ready[[1]], tt.params )
## End(Not run)
}

\keyword{predict}
\keyword{out-of-bag}
