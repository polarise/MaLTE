\name{MaLTE-package}
\alias{MaLTE-package}
\alias{MaLTE}
\docType{package}
\title{Package to perform Machine Learning of Transcript Expression (MaLTE)}
\description{In routine analysis of Affymetrix microarrays, probe intensities are first subjected to background correction and quantile-normalisation before being summarised into single gene expression estimates. Two of the most widely used summarisation algorithms are \bold{median-polish} used in robust multiarray averaging (RMA) and \bold{probe logarithmic intensity error} (PLIER). These methods suffer from two major limitations: one, they provide relative, as opposed to absolute, expression estimates--the estimates for two genes within the same sample should not be compared; and two, they do not extend to quantifying the abundance of individual transcript isoforms. Both of these are overcome by high-throughput sequencing (HTS) quantification techniques such as RNA-Seq.

We have supplemented the summarisation algorithm by using a statistical learning approach. Our method, named Machine Learning of Transcript Expression (MaLTE), learns a regression model between RNA-Seq expression estimates and the microarray probe intensities that may later be used to predict the expression from probe intensities alone. Our current implementation uses conditional random forest (CRF) as the learning framework but is open to developments in machine learning algorithms.}
\details{
  \tabular{ll}{
    Package: \tab MaLTE\cr
    Type: \tab Package\cr
    Version: \tab 0.1\cr
    Date: \tab 2013-05-13\cr
    License: \tab GPL (>= 3)\cr
  }
}
\examples{
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
}
\author{
Paul Korir <paul.korir@gmail.com>

Maintainer: Paul Korir <paul.korir@gmail.com>
}
