\name{gene.id}
\alias{gene.id}

\title{
Get the gene identifier of objects of class \code{TT.Seq}.
}

\description{
Accessor function for gene identifier for \code{TT.Seq} objects.
}

\usage{
gene.id(object)
}

\arguments{
  \item{object}{an object of class \code{TT.Seq}}
}

\value{
Returns the gene identifier.
}

\author{
Paul K. Korir
}

\seealso{
\code{link{tx.id}}
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

# read in the training and test data
tt.ready <- read.data( train.fn="train_data.txt.gz", test.fn="test_data.txt.gz" )

# params 
tt.params <- TT.Params()

tt.seq <- array2seq( tt.ready, tt.params )

id <- gene.id( tt.seq[[1]] )
ids <- as.vector( sapply( tt.seq, gene.id ))
## End(Not run)
}

\keyword{gene}
\keyword{identifier}
\keyword{ID}
