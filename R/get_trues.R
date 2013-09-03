get.trues <- function( tt.seq, sample.names )
{
	if ( is( tt.seq[[1]], "TT.Seq.Gene" ))
	{
		# first deal with the trues
		tt.p1 <- lapply( tt.seq,
			function( m )
			{
				if ( length( m@trues ) != length( sample.names ) )
				{
					if ( is.na( m@trues ))
					{
						rep( NA, m@no.samples )
					} else
					{
						stop( paste( "The number of true samples for gene ", m@gene.id, " is different from the number of sample names. Ensure that when you run get.names( ) you set test to either 'TRUE' for test samples or 'FALSE' for training samples (e.g. OOB).\n", sep="" ) )
					}
				} else
				{
					trues( m )
				}
			}
		)
		tt.p2 <- unlist( tt.p1, recursive=FALSE, use.names=FALSE )
		tt.p3 <- as.vector( tt.p2 )
		tt.p4 <- matrix( tt.p3, nrow=length( tt.seq ), byrow=TRUE )
		tt.trues <- data.frame( tt.p4 )
		colnames( tt.trues ) <- sample.names

		# next deal with the gene names
		tt.g1 <- lapply( tt.seq, gene.id )
		tt.g2 <- unlist( tt.g1 )
		tt.gene.id <- as.vector( tt.g2 )

		# put them all together
		tt.trues.df <- cbind( gene_id=tt.gene.id, tt.trues )
	} else
	if ( is( tt.seq[[1]], "TT.Seq.Tx" ))
	{
		# first deal with the transcript names
		tt.t1 <- lapply( tt.seq, function( m ){ m@tx.id } )
		tt.t2 <- unlist( tt.t1 )
		txs.names <- as.vector( tt.t2 )
				
		# next deal with the trues
		tt.p1 <- lapply( tt.seq,
			function( m )
			{
				if ( m@no.txs == 1 )
				{
					no <- length( m@trues )
				} else
				{
					no <- nrow( m@trues )
				}
				
				if ( no != length( sample.names ) )
				{
					if ( is.na( m@trues ) )
					{
						rep( NA, m@no.samples*m@no.txs )
					} else
					{
						stop( paste( "The number of true samples for gene ", m@gene.id, " is different from the number of sample names. Ensure that when you run get.names( ) you set test to either 'TRUE' for test samples or 'FALSE' for training samples (e.g. OOB).\n", sep="" ) )
					}
				} else
				{
					trues( m )
				}
			}
		)
	
		tt.p2 <- unlist( tt.p1, recursive=FALSE, use.names=FALSE )
		tt.p3 <- as.vector( tt.p2 )
		tt.p4 <- matrix( tt.p3, nrow=length( txs.names ), byrow=TRUE )
		tt.trues <- data.frame( tt.p4 )
		colnames( tt.trues ) <- sample.names
		
		# put them all together
		tt.trues.df <- cbind( tx_id=txs.names, tt.trues )
	}

	return( tt.trues.df )
}
