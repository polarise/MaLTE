#
#
#
.subsetter <- function( object, n, trues.present=FALSE )
{
	# convert all NAs to FALSE
	n[ is.na( n ) ] = FALSE
	
	m <- list(
		gene.id = object@gene.id,
		tx.id = object@tx.id[n],
		no.txs = length( object@tx.id[n] ),
		no.samples = object@no.samples,
		trues = object@trues,
		predictions = if ( sum( n ) == 0 ){ object@predictions[n] } else if ( sum( n ) > 0 ){ object@predictions[,n] },
		cor.S = if ( trues.present ){ object@cor.S[n] } else { NA },
		cor.S.pv = if ( trues.present ){ object@cor.S.pv[n] } else { NA },
		cor.P = if ( trues.present ){ object@cor.P[n] } else { NA },
		cor.P.pv = if ( trues.present ){ object@cor.P.pv[n] } else { NA },
		means = object@means[n],
		vars = object@vars[n],
		OOB = object@OOB
	)
	
	m.obj <- TT.Seq.Tx( m )
	
	return( m.obj )
}

#
#
#
oob.filter <- function( list.objects, list.objects.oob, thresh=0, method="pearson" )
{
	# make sure that both lists are non-empty
	
	# make sure that both lists are pairs
	
	# make sure that list.objects is not OOB
	
	
	# make sure that list.objects.oob is OOB

	# check the class of object
	if ( is( list.objects[[1]], "TT.Seq.Gene" ) & is( list.objects.oob[[1]], "TT.Seq.Gene" ) )
	{
		if ( method == "pearson" )
		{
			cor.P.oob <- as.vector( sapply( list.objects.oob, function( m ){ m@cor.P })) # get the OOB Pearsons
			cor.P.oob.pv <- as.vector( sapply( list.objects.oob, function( m ){ m@cor.P.pv })) # get the OOB Pearson p-values
			tt.filtered <- list.objects[ which( cor.P.oob > thresh & cor.P.oob.pv <= 0.05 ) ]
		} else if ( method == "spearman" )
		{
			cor.S.oob <- as.vector( sapply( list.objects.oob, function( m ){ m@cor.S })) # get the OOB Spearman
			cor.S.oob.pv <- as.vector( sapply( list.objects.oob, function( m ){ m@cor.S.pv })) # get the OOB Spearman p-values	
			tt.filtered <- list.objects[ which( cor.S.oob > thresh & cor.S.oob.pv <= 0.05 ) ]
		} else
		{
			stop( "Error: incorrect value in method argument; use either 'pearson' (default) or 'spearman'" )
		}
	} else
	if ( is( list.objects[[1]], "TT.Seq.Tx" ) &	is( list.objects.oob[[1]], "TT.Seq.Tx" ) )
	{
		cor.P.oob <- lapply( list.objects.oob, function( m ){ m@cor.P > thresh })
		
		trues.present <- all( !is.na( list.objects[[1]]@trues ))
		
		# subset
		tt.filtered1 <- mapply( .subsetter, list.objects, cor.P.oob, trues.present=trues.present, SIMPLIFY=FALSE )
		
		# TRUE - retain; FALSE - exclude
		tt.filtered2 <- as.vector( sapply( tt.filtered1, function( m )
					{
						if ( m@no.txs > 0 )
						{
							return( TRUE )
						} else
						{
							return( FALSE )
						}
					}
				)
			)
		
		# filter
		tt.filtered <- tt.filtered1[ tt.filtered2 ]
	}
	
	return( tt.filtered )
}

