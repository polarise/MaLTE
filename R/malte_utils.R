.paired_correlations <- function( gene, other.genes, affy.genes, other.list, affy.list )
{
	#
	# input:
	# - gene of interest
	# - vector of MaLTE genes
	# - vector of median-polish genes
	# - list of MaLTE data with predictions per gene
	# - list of median-polish estimates per gene
	#
	# output:
	# - Pearson and Spearman correlation coefficients for median-polish and MaLTE against RNA-Seq
	#
	
	# get indexes
	i <- which( other.genes == gene )
	j <- which( affy.genes == gene )

	if ( length( i ) == 0 )	# absent in MaLTE; present in Affy
	{
		cor.S = NA		# because the gene is absent
		cor.P = NA
		affy.cor.S = NA	# because there are no trues
		affy.cor.P = NA
	}	else if ( length( j ) == 0 ) # absent in Affy
	{
		cor.P <- other.list[[i]]@cor.P
		cor.S <- other.list[[i]]@cor.S
		affy.cor.S = NA 	# because it's absent in Affy
		affy.cor.P = NA
	} else
	{		# both present
		trues <- other.list[[i]]@trues
		
		affy <- affy.list[[j]]$affy
		cor.P <- other.list[[i]]@cor.P
		cor.S <- other.list[[i]]@cor.S
		
		# affy correlations
		affy.cor.P <- cor.test( affy, trues )$estimate
		affy.cor.S <- cor.test( affy, trues, me="sp" )$estimate
	}
	
	list( genes=gene, affy.cor.S=affy.cor.S, affy.cor.P=affy.cor.P, cor.S=cor.S, cor.P=cor.P )
}

.df_to_list <- function( m )
{
	list( genes=as.vector( m[1] ), affy=as.vector( as.numeric( m[-1] )))
}

.add.trues.affy.list <- function( affy.mixed.row, mixed, genes )
{
	gene <- affy.mixed.row$genes
	i <- which( genes == gene )
	
	affy <- affy.mixed.row$affy
	
	if ( length( i ) == 0 )
	{
		list( genes=gene, trues=NA, predictions=NA )
	} else
	{
		trues <- mixed[[i]]@trues
		list( genes=gene, trues=trues, predictions=affy )
	}
}

compare.correlations <- function( tt.seq, affy.fn )
{
	# read in the Affymetrix RMA/PLIER (or other summarisation) results
	affy.mixed.test <- read.table( affy.fn, header=T, stringsAsFactors=F )
	
	# convert it to a list
	affy.mixed.test.list <- apply( affy.mixed.test, 1, .df_to_list )
	
	# get the gene names in the Affymetrix data
	affy.genes <- as.vector( sapply( affy.mixed.test.list, function( m ){ m$genes }))

	# get the MaLTE gene names
	genes <- as.vector( sapply( tt.seq, function( m ){ m@gene.id } ))

	# add the true (RNA-Seq) expression values to the Affymetrix data
	aug.affy.mixed.test.list <- mclapply( affy.mixed.test.list, .add.trues.affy.list, tt.seq, genes  )

	# perform gene-wise comparisons pair-wise
	tt.seq.compared <- mclapply( intersect( genes, affy.genes ), .paired_correlations, genes, affy.genes, tt.seq, affy.mixed.test.list )

	# collate the data
	all.genes <- as.vector( sapply( tt.seq.compared, function( m ){ m$genes }))

	affy.cor.S <- as.vector( sapply( tt.seq.compared, function( m ){ m$affy.cor.S }))
	affy.cor.P <- as.vector( sapply( tt.seq.compared, function( m ){ m$affy.cor.P }))

	cor.S <- as.vector( sapply( tt.seq.compared, function( m ){ m$cor.S }))
	cor.P <- as.vector( sapply( tt.seq.compared, function( m ){ m$cor.P }))
	
	# make the data.frame
	tt.seq.compared.df <- data.frame( all.genes, affy.cor.S, affy.cor.P, cor.S, cor.P )

	# return
	return( tt.seq.compared.df )
}

cors.P <- function( tt.seq )
{
	cors <- as.vector( unlist( lapply( tt.seq, cor.P )))
	return( cors )
}

cors.S <- function( tt.seq )
{
	cors <- as.vector( unlist( lapply( tt.seq, cor.S )))
	return( cors )
}
