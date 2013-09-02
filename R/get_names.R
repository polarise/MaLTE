get.names <- function( samples.fn="samples.txt", test=TRUE )
{
	# validation?
	# does the file exist?
	
	samples <- read.table( samples.fn, header=T, stringsAsFactors=F )
	
	# either there are NA's or not
	if ( sum( is.na( samples$hts )) != 0 )
	{
		test.indexes <- is.na( samples$hts )
	} else if ( sum( is.na( samples$hts )) == 0 )
	{
		test.indexes <- as.vector( sapply( samples$hts, function( m ){ grepl( "^\\*", m )} ))
	} else
	{
		stop( "Possibly malformed sample names file. Check that all ma samples are present and that hts samples are either 'NA' or '*<sample_name>'." )
	}
	
	if ( test )
	{
		names <- samples$ma[ test.indexes ]
	} else
	{
		names <- samples$ma[ !test.indexes ]
	}

	return( names )
}

get.train <- function( samples.fn="samples.txt", test=FALSE ) get.names( samples.fn, test )

get.test <- function( samples.fn="samples.txt", test=TRUE ) get.names( samples.fn, test )
