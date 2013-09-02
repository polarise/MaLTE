\name{oob.filter}
\alias{oob.filter}

\title{
Filter a list of genes using a threshold on the correlation coefficient of OOB estimates
}

\description{
Takes two lists of TT.Seq.Gene objects: using test data and using OOB estimates. It derives a list of genes in the former that exceed some threshold of correlation in the OOB expression estimates.
}

\usage{
oob.filter(list.objects, list.objects.oob, thresh = 0)
}

\arguments{
  \item{list.objects}{a list of TT.Seq.Gene or TT.Seq.Tx objects}
  \item{list.objects.oob}{a list of TT.Seq.Gene or TT.Seq.Tx objects computed using the OOB estimates only}
  \item{thresh}{the threshold for correlation estimates; the default value thresh=0}
}

\value{
Returns a list of objects of classes \code{TT.Seq.Gene} or \code{TT.Seq.Tx} filtered by OOB correlation.
}

\author{
Paul K. Korir
}

\seealso{
\code{\link{array2seq}}, \code{\link{array2seq.oob}}
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

# filter by OOB at default threshold (thresh=0)
tt.filtered <- oob.filter( tt.seq, tt.seq.oob )
}

\keyword{filter}
\keyword{out-of-bag}
