\name{array2seq.oob}
\alias{array2seq.oob}

\title{
Predict the out-of-bag (OOB) expression estimates using MaLTE on a list of TT.Ready.Gene or TT.Ready.Tx objects.
}

\description{
Takes a list of TT.Ready.Gene or TT.Ready.Tx objects and for each trains a conditonal random forest model then predicts an out-of-bag (OOB) estimate.
}

\usage{
array2seq.oob(tt.ready, params.object)
}

\arguments{
  \item{tt.ready}{a list of TT.Ready.Gene or TT.Ready.Tx objects}
  \item{params.object}{an object of class TT.Params that contains the parameters to be passed to the conditional random forest}
}

\value{
Returns OOP predictions as a list of objects of classes \code{TT.Seq.Gene} or \code{TT.Seq.Tx}.
}

\author{
Paul K. Korir
}

\seealso{
\code{link{array2seq}}
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

# set default parameters
params.object <- TT.Params()

# predict
tt.seq <- array2seq( tt.ready, params.object )

# perform OOB predictions
tt.seq.oob <- array2seq.oob( tt.ready, params.object )
}

\keyword{predict}
