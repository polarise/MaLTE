#
# create a list of TT.Ready.Gene objects
#
.numerise <- function( m )	# function to make all variables into numbers
{
	if ( is.na( m$hts.test ) )
	{
		hts.test <- NA
	} else 
	{
		hts.test <- as.vector( unlist( lapply( strsplit( m$hts.test, "," ), as.numeric )))
	}
	hts.train <- as.vector( unlist( lapply( strsplit( m$hts.train, "," ), as.numeric )))
	
	if ( is.na( m$probes.train ) | is.na( m$probes.test ) )
	{	
	   probes.train <- NA
	   probes.test <- NA
	} else
	{
		probes.train <- matrix( as.vector( unlist( lapply( strsplit( m$probes.train, "," ), as.numeric ))), nrow=m$no.train )
		probes.test <- matrix( as.vector( unlist( lapply( strsplit( m$probes.test, "," ), as.numeric ))), nrow=m$no.test )
	}
	
	m.list <- list( gene.id=m$gene.id, no.train=m$no.train, no.test=m$no.test, no.probes=m$no.probes, hts.train=hts.train, hts.test=hts.test, probes.train=probes.train, probes.test=probes.test )
	
	m.obj <- TT.Ready.Gene( m.list )
	
	return( m.obj )
}

#
# create a list of TT.Ready.Tx objects
#
.numerise.txs <- function( m )
{
	if ( is.na( m$hts.test ) )
	{
		hts.test <- NA
	} else
	{
		hts.test <- matrix( as.vector( unlist( lapply( strsplit( m$hts.test, "," ), as.numeric ))), nrow=m$no.test )
	}
	hts.train <- matrix( as.vector( unlist( lapply( strsplit( m$hts.train, "," ), as.numeric ))), nrow=m$no.train )
	
	if ( is.na( m$probes.train ))
	{	
	   probes.train <- NA
	   probes.test <- NA
	} 
	else
	{
		probes.train <- matrix( as.vector( unlist( lapply( strsplit( m$probes.train, "," ), as.numeric ))), nrow=m$no.train )
		probes.test <- matrix( as.vector( unlist( lapply( strsplit( m$probes.test, "," ), as.numeric ))), nrow=m$no.test )
	}
	tx.id <- as.vector( unlist( strsplit( m$tx.id, "," )))
	
	m.list <- list( gene.id=m$gene.id, tx.id=tx.id, no.txs=m$no.txs, no.train=m$no.train, no.test=m$no.test, no.probes=m$no.probes, hts.train=hts.train, hts.test=hts.test, probes.train=probes.train, probes.test=probes.test )
	
	m.obj <- TT.Ready.Tx( m.list )
	
	return( m.obj )
}

#
# genes
#
read.data <- function( train.fn, test.fn )
{
	cat( "Reading training data from '", train.fn, "'\n" )
	train.data <- read.table( train.fn, stringsAsFactors=F )
	
	cat( "Reading testing data from '", test.fn, "'\n" )
	test.data <- read.table( test.fn, stringsAsFactors=F )
	
	tt <- cbind( train.data[,c(1,2,3,5,6)], test.data[,c(2,5,6)] )
	colnames( tt ) <- c( "gene.id", "no.train", "no.probes", "hts.train", "probes.train", "no.test", "hts.test", "probes.test" )

	# perform the following in parallel
	cat( "Creating list of objects ready for prediction...\n" )
	tt.data <- mclapply( split(tt, 1:nrow(tt)), as.list )	# convert the data.frame to a list of lists
	tt.ready <- mclapply( tt.data, .numerise )								# convert strings to vectors of numbers

	return( tt.ready )
}

#
# transcripts
#
read.txs.data <- function( train.fn, test.fn )
{
	cat( "Reading training data from '", train.fn, "'\n" )
	train.data <- read.table( train.fn, stringsAsFactors=F )
	
	cat( "Reading testing data from '", test.fn, "'\n" )
	test.data <- read.table( test.fn, stringsAsFactors=F )
	
	tt <- cbind( train.data[,c(1,2,3,4,5,7,8)], test.data[,c(4,7,8)] )
	colnames( tt ) <- c( "gene.id", "tx.id", "no.txs", "no.train", "no.probes", "hts.train", "probes.train", "no.test", "hts.test", "probes.test" )

	# perform the following in parallel
	cat( "Creating list of objects ready for prediction...\n" )
	tt.data <- mclapply( split(tt, 1:nrow(tt)), as.list )	# convert the data.frame to a list of lists
	tt.ready.txs <- mclapply( tt.data, .numerise.txs )				# convert strings to vectors of numbers

	return( tt.ready.txs )
}

