\name{train.and.predict}
\alias{train.and.predict}

\title{
Function that performs training and prediction for each gene.
}

\description{
Training is currently implemented using the conditional random forest algorithm but is open to improvements in machine learning technology. For each gene, this function learns a model based on HTS expression estimates and microarray probe intensities then either predicts the HTS expression or an out-of-bag estimate for a completely new set of probes of the same gene. The function \code{run} is a member function of TT.Ready.Gene class and calls \code{train.and.predict}.
}

\usage{
train.and.predict(object, params.object=NULL, gene.tuned=TRUE, tune.quantreg=FALSE, tune.verbose=FALSE, OOB=FALSE)
}

\arguments{
  \item{object}{an object of class \code{TT.Ready.Gene}}
  \item{params.object}{\code{NULL} or an object of class \code{TT.Params}}
  \item{gene.tuned}{logical; if \code{TRUE} then tuned parameters for this gene will be determined}
  \item{tune.quantreg}{logical; if \code{TRUE} then will use \code{quantregForest} instead of \code{cforest} random forest implementation}
  \item{tune.verbose}{logical; if \code{TRUE} then will print tuning information; should not be used when running in parallel}
  \item{OOB}{logical; if \code{TRUE} then will estimate out-of-bag values instead of using test data}
}

\value{
Returns an object of class \code{TT.Seq.Gene}.
}

\author{
Paul K. Korir
}

\note{
This function is used by \code{array2seq} expression prediction function.
}

\seealso{
\code{\link{train.and.predict.txs}}, \code{\link{array2seq}}
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

# default parameters
params.object <- TT.Params()

# train and predict
a.gene <- train.and.predict( tt.ready[[1]], params.object )
a.gene
predictions( a.gene )
sample.names <- get.names( samples.fn )
}
}

\keyword{train}
\keyword{predict}
\keyword{gene}
