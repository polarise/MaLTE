\name{show.probes.test}
\alias{show.probes.test}

\title{
Show microarray probes test data.
}

\description{
Accessor function for microarray probes test data \code{TT.Ready} objects.
}

\usage{
show.probes.test(object)
}

\arguments{
  \item{object}{an object of class \code{TT.Ready}}
}

\value{
Returns a matrix of microarray probes test data.
}

\author{
Paul K. Korir
}

\seealso{
\code{link{show.hts.test}}, \code{\link{show.hts.train}}, \code{\link{show.probes.train}}
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

hts.test <- show.probes.test( tt.ready[[1]] )
}

\keyword{probes}
\keyword{test}
