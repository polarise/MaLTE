get.predictions <- function( tt.seq, sample.names, qnorm=FALSE )
{
	# check that tt.seq is not empty

	# validation

	if ( is( tt.seq[[1]], "TT.Seq.Gene" ))
	{
		# first deal with predictions
		tt.p1 <- lapply( tt.seq,
			function( m ){
				if ( length( m@predictions ) != length( sample.names ) )
				{
					if ( is.na( m@predictions ))
					{
						rep( NA, m@no.samples )
					} else
					{
						stop( paste( "The number of samples predicted for gene ", m@gene.id, " is different from the number of sample names. Ensure that when you run get.names( ) you set test to either 'TRUE' for test samples or 'FALSE' for training samples (e.g. OOB). Alternatively, if you are running get.predictions() on tt.seq instead of on tt.filtered it is likely that this gene has no predictions (NA).\n", sep="" ) )
					}
				} else
				{
					predictions( m )
				}
			}
		)
		tt.p2 <- unlist( tt.p1, recursive=FALSE, use.names=FALSE )
		tt.p3 <- as.vector( tt.p2 )
		tt.p4 <- matrix( tt.p3, nrow=length( tt.seq ), byrow=TRUE )
		tt.predictions <- data.frame( tt.p4 )
		
		if ( qnorm )
		{
			cat( "Note: predictions will be quantile normalised.\n" )
			tt.predictions = normalizeQuantiles( tt.predictions, ties=FALSE ) # do not allow ties
		}
		colnames( tt.predictions ) <- sample.names

		# next deal with the gene names
		tt.g1 <- lapply( tt.seq, gene.id )
		tt.g2 <- unlist( tt.g1 )
		tt.gene.id <- as.vector( tt.g2 )

		# put them all together
		tt.predictions.df <- cbind( gene_id=tt.gene.id, tt.predictions )
	} else
	if ( is( tt.seq[[1]], "TT.Seq.Tx" ))
	{
		# first deal with the transcript names
		tt.t1 <- lapply( tt.seq, function( m ){ m@tx.id } )
		tt.t2 <- unlist( tt.t1 )
		txs.names <- as.vector( tt.t2 )
				
		# next deal with the predictions
		tt.p1 <- lapply( tt.seq,
			function( m )
			{
				if ( m@no.txs == 1 )
				{
					no <- length( m@predictions )
				} else
				{
					no <- nrow( m@predictions )
				}
				
				if ( no != length( sample.names ) )
				{
					if ( is.na( m@predictions ) )
					{
						rep( NA, m@no.samples*m@no.txs )
					} else
					{
						stop( paste( "The number of samples predicted for gene ", m@gene.id, " is different from the number of sample names. Ensure that when you run get.names( ) you set test to either 'TRUE' for test samples or 'FALSE' for training samples (e.g. OOB). Alternatively, if you are running get.predictions() on tt.seq instead of on tt.filtered it is likely that this gene has no predictions (NA).\n", sep="" ) )
					}
				} else
				{
					predictions( m )
				}
			}
		)
		tt.p2 <- unlist( tt.p1, recursive=FALSE, use.names=FALSE )
		tt.p3 <- as.vector( tt.p2 )
		tt.p4 <- matrix( tt.p3, nrow=length( txs.names ), byrow=TRUE )
		tt.predictions <- data.frame( tt.p4 )
		
		if ( qnorm )
		{
			cat( "Note: predictions will be quantile normalised.\n" )
			tt.predictions = normalizeQuantiles( tt.predictions, ties=FALSE )
		}
		
		colnames( tt.predictions ) <- sample.names
		
		# put them all together
		tt.predictions.df <- cbind( tx_id=txs.names, tt.predictions )
	}

	return( tt.predictions.df )
}

get.lower <- function( tt.seq, sample.names )
{
	# validation

	if ( is( tt.seq[[1]], "TT.Seq.Gene" ))
	{
		# first deal with predictions
		tt.p1 <- lapply( tt.seq,
			function( m ){
				if ( length( m@predictions.lower ) != length( sample.names ) )
				{
					if ( is.na( m@predictions.lower ))
					{
						rep( NA, m@no.samples )
					} else
					{
						stop( paste( "The number of samples predicted for gene ", m@gene.id, " is different from the number of sample names. Ensure that when you run get.names( ) you set test to either 'TRUE' for test samples or 'FALSE' for training samples (e.g. OOB). Alternatively, if you are running get.predictions() on tt.seq instead of on tt.filtered it is likely that this gene has no predictions (NA).\n", sep="" ) )
					}
				} else
				{
					lower( m )
				}
			}
		)
		tt.p2 <- unlist( tt.p1, recursive=FALSE, use.names=FALSE )
		tt.p3 <- as.vector( tt.p2 )
		tt.p4 <- matrix( tt.p3, nrow=length( tt.seq ), byrow=TRUE )
		tt.predictions <- data.frame( tt.p4 )
		colnames( tt.predictions ) <- sample.names

		# next deal with the gene names
		tt.g1 <- lapply( tt.seq, gene.id )
		tt.g2 <- unlist( tt.g1 )
		tt.gene.id <- as.vector( tt.g2 )

		# put them all together
		tt.predictions.df <- cbind( gene_id=tt.gene.id, tt.predictions )
	}
	return( tt.predictions.df )
}

get.upper <- function( tt.seq, sample.names )
{
	# validation

	if ( is( tt.seq[[1]], "TT.Seq.Gene" ))
	{
		# first deal with predictions
		tt.p1 <- lapply( tt.seq,
			function( m ){
				if ( length( m@predictions.upper ) != length( sample.names ) )
				{
					if ( is.na( m@predictions.upper ))
					{
						rep( NA, m@no.samples )
					} else
					{
						stop( paste( "The number of samples predicted for gene ", m@gene.id, " is different from the number of sample names. Ensure that when you run get.names( ) you set test to either 'TRUE' for test samples or 'FALSE' for training samples (e.g. OOB). Alternatively, if you are running get.predictions() on tt.seq instead of on tt.filtered it is likely that this gene has no predictions (NA).\n", sep="" ) )
					}
				} else
				{
					upper( m )
				}
			}
		)
		tt.p2 <- unlist( tt.p1, recursive=FALSE, use.names=FALSE )
		tt.p3 <- as.vector( tt.p2 )
		tt.p4 <- matrix( tt.p3, nrow=length( tt.seq ), byrow=TRUE )
		tt.predictions <- data.frame( tt.p4 )
		colnames( tt.predictions ) <- sample.names

		# next deal with the gene names
		tt.g1 <- lapply( tt.seq, gene.id )
		tt.g2 <- unlist( tt.g1 )
		tt.gene.id <- as.vector( tt.g2 )

		# put them all together
		tt.predictions.df <- cbind( gene_id=tt.gene.id, tt.predictions )
	}
	return( tt.predictions.df )
}
