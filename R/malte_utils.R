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
	affy.mixed.test <- read.table( affy.fn, header=T, stringsAsFactors=F, check.names=F, comment.char="#" )
	
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

within_correlations <- function( tt.seq, affy.fn, affy.name="median_polish", check.names=TRUE  )
{
	# make sure the affy.name is valid
	if ( check.names )
		affy.name <- make.names( affy.name ) # 'median-polish' becomes 'median.polish'
	
	# read in the Affymetrix RMA/PLIER (or other summarisation) results
	affy.mixed.test <- read.table( affy.fn, header=T, stringsAsFactors=F, check.names=F, comment.char="#" )
	
	# convert it to a list
	affy.mixed.test.list <- apply( affy.mixed.test, 1, .df_to_list )
	
	# get the gene names in the Affymetrix data
	affy.genes <- as.vector( sapply( affy.mixed.test.list, function( m ){ m$genes }))

	# get the MaLTE gene names
	genes <- as.vector( sapply( tt.seq, function( m ){ m@gene.id } ))

	# add the true (RNA-Seq) expression values to the Affymetrix data
	aug.affy.mixed.test.list <- mclapply( affy.mixed.test.list, .add.trues.affy.list, tt.seq, genes  )

	# data frames for results
	within.P <- data.frame()
	within.S <- data.frame()

	# vectors for intermediate results
	within.cors.P <- c()
	within.cors.S <- c()
	affy.within.cors.P <- c()
	affy.within.cors.S <- c()

	# determine the number of samples involved
	no.samples <- tt.seq[[1]]@no.samples

	# calculate the within-sample correlation
	# for each sample...
	for ( i in 1:no.samples )
	{
		# get the true estimates (RNA-Seq) and MaLTE predictions from tt.seq list
		samp.trues <- as.vector( unlist( lapply( tt.seq, function( m ){ m@trues[i] })))
		samp.preds <- as.vector( unlist( lapply( tt.seq, function( m ){ m@predictions[i] })))
	
		# get the true estimates and summarisations from aug.affy.mixed.test.list
		affy.trues <- as.vector( unlist( lapply( aug.affy.mixed.test.list, function( m ){ m$trues[i] })))
		affy.preds <- as.vector( unlist( lapply( aug.affy.mixed.test.list, function( m ){ m$predictions[i] })))
	
		# compute respective correlations for MaLTE
		cor.P <- as.vector( cor.test( samp.trues, samp.preds )$estimate )
		cor.S <- as.vector( cor.test( samp.trues, samp.preds, me="sp" )$estimate )
	
		# compute respective correlations for summarisation
		affy.cor.P <- as.vector( cor.test( affy.trues, affy.preds )$estimate )
		affy.cor.S <- as.vector( cor.test( affy.trues, affy.preds, me="sp" )$estimate )
	
		# append correlations to the data frame for MaLTE
		within.cors.P <- c( within.cors.P, cor.P )
		within.cors.S <- c( within.cors.S, cor.S )
	
		# append correlations to the data frame for summarisation
		affy.within.cors.P <- c( affy.within.cors.P, affy.cor.P )
		affy.within.cors.S <- c( affy.within.cors.S, affy.cor.S )
	}

	# some ad hoc stuff
	# put the Pearsons together
	within.P <- rbind( within.P, within.cors.P )
	within.P <- rbind( within.P, affy.within.cors.P )

	# put the Spearmans together
	within.S <- rbind( within.S, within.cors.S )
	within.S <- rbind( within.S, affy.within.cors.S )

	# add sample labels
	colnames( within.P ) <- paste( "S", 1:no.samples, sep="" )
	colnames( within.S ) <- paste( "S", 1:no.samples, sep="" )

	# transpose
	within.P <- data.frame( t( within.P ))
	within.S <- data.frame( t( within.S ))

	# add method labels
	names <- c( "MaLTE", affy.name )

	colnames( within.P ) <- names
	colnames( within.S ) <- names

	# return as a list
	list( Pearson=within.P, Spearman=within.S )
}

