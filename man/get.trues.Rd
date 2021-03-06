\name{get.trues}
\alias{get.trues}

\title{
Get a data.frame of the true expression estimates.
}

\description{
The sample names will be required in the output data when carrying out downstream analyes. This function takes the samples file ('samples.txt') and extracts the required names.
}

\usage{
get.trues(tt.seq, sample.names)
}

\arguments{
  \item{tt.seq}{an list of objects of class TT.Seq.Gene or TT.Seq.Tx}
  \item{sample.names}{a vector of the sample names; the output of get.names(samples.fn)}
}

\value{
Returns a data.frame of true expression estimates (RNA-Seq). The first column is the gene identifier field.
}

\author{
Paul K. Korir
}

\note{
May be run on both filtered or non-filtered \code{tt.seq} lists.
}

\seealso{
\code{\link{get.names}}
}

\examples{
\dontrun{
# location of example input files
samples.fn <- paste( system.file( package="MaLTE" ), "data", "samples.txt.gz", sep="/" )
ma.fn <- paste( system.file( package="MaLTE" ), "data", "ma_data.txt.gz", sep="/" )
raw_ma.fn <- paste( system.file( package="MaLTE" ), "data", "raw_ma_data.txt.gz", sep="/" )
hts.fn <- paste( system.file( package="MaLTE" ), "data", "hts_data.txt.gz", sep="/" )
g2p.fn <- paste( system.file( package="MaLTE" ), "data", "gene_probesets.txt.gz", sep="/" )

# prepare training and test data
prepare.data( samples.fn=samples.fn, ma.fn=ma.fn, hts.fn=hts.fn, g2p.fn=g2p.fn )

# use raw APT output for microarray data
prepare.data( samples.fn=samples.fn, ma.fn=raw_ma.fn, hts.fn=hts.fn, g2p.fn=g2p.fn, raw=TRUE )

# read in the training and test data
tt.ready <- read.data( train.fn="train_data.txt.gz", test.fn="test_data.txt.gz" )

# set default parameters
params.object <- TT.Params()

# predict
tt.seq <- array2seq( tt.ready, params.object )

# perform OOB predictions
tt.seq.oob <- array2seq.oob( tt.ready, params.object )

# filter by OOB at default threshold (thresh=0)
tt.filtered <- oob.filter( tt.seq, tt.seq.oob )

# get the test sample names
sample.names <- get.names( samples.fn )
sample.names <- get.test( samples.fn )

# get the predicted expression estimates after OOB filtering
df.preds <- get.predictions( tt.filtered, sample.names )

# get the true expression estimates before OOB filtering
df.trues <- get.trues( tt.seq, sample.names )
}
}

\keyword{predictions}
