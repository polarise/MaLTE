\name{show.hts.train}
\alias{show.hts.train}

\title{
Show HTS training data.
}

\description{
Accessor function for HTS training data \code{TT.Ready} objects.
}

\usage{
show.hts.train(object)
}

\arguments{
  \item{object}{an object of class \code{TT.Ready}}
}

\value{
Returns a matrix of HTS training data.
}

\author{
Paul K. Korir
}

\seealso{
\code{link{show.hts.test}}, \code{\link{show.probes.train}}, \code{\link{show.probes.test}}
}

\examples{
# location of example input files
samples.fn <- paste( system.file( package="MaLTE" ), "data", "samples.txt.gz", sep="/" )
ma.fn <- paste( system.file( package="MaLTE" ), "data", "ma_data.txt.gz", sep="/" )
hts.fn <- paste( system.file( package="MaLTE" ), "data", "hts_data.txt.gz", sep="/" )
g2p.fn <- paste( system.file( package="MaLTE" ), "data", "gene_probesets.txt.gz", sep="/" )

# prepare training and test data
prepare.data( samples.fn=samples.fn, ma.fn=ma.fn, hts.fn=hts.fn, g2p.fn=g2p.fn )

# read in the training and test data
tt.ready <- read.data( train.fn="train_data.txt.gz", test.fn="test_data.txt.gz" )

hts.test <- show.hts.train( tt.ready[[1]] )
}

\keyword{hts}
\keyword{train}
