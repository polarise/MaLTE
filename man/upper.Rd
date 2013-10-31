\name{upper}
\alias{upper}

\title{
Extract upper bound of predictions
}

\description{
Accessor function for upper bound for predicted expression given objects of 
class \code{TT.Seq}.
}

\usage{
upper(object)
}

\arguments{
  \item{object}{an object of class \code{TT.Seq}}
}

\value{
Returns either a vector or matrix of upper bound of predicted values.
}

\note{
May be run on both filtered and non-filtered \code{tt.seq} lists. Only returns 
values if \code{quantregForest} is used (i.e. 
\code{array2seq(tt.ready, params.object=TT.Params(quantreg=TRUE))}).
}

\author{
Paul K. Korir
}

\seealso{
\code{link{get.upper}}
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

# params
tt.params <- TT.Params( quantreg=TRUE )

tt.seq <- array2seq( tt.ready, tt.params )

# get predictions
predictions.upper <- upper( tt.seq[[1]] )
}
}

\keyword{predict}
\keyword{out-of-bag}
