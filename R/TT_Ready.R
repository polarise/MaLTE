#-------------------------------------------------------------------------------
#
# Union of classes numeric, logical and matrix
#
#-------------------------------------------------------------------------------
setClassUnion( "numericORlogicalORmatrix", c( "numeric", "logical", "matrix" ))
# hts.test is either numeric or logical hence the class union

setClassUnion( "numericORmatrix", c( "numeric", "matrix" ))
# hts.train is either numeric or logical hence the class union

#-------------------------------------------------------------------------------
#
# TT.Ready Abstract Base Class
#
#-------------------------------------------------------------------------------
setClass( "TT.Ready", representation( gene.id="character", no.train="numeric", no.test="numeric", no.probes="numeric", hts.train="numericORmatrix", hts.test="numericORlogicalORmatrix", probes.train="matrix", probes.test="matrix" ), "VIRTUAL" )

#-------------------------------------------------------------------------------
#
# TT.Ready.Gene Class
#
#-------------------------------------------------------------------------------
setClass( "TT.Ready.Gene", contains="TT.Ready" )

#
# TT.Ready.Gene constructor
#
TT.Ready.Gene <- function( m )
{
	# VALIDATIONS	
	# length of hts.train == no.train
	if ( length( m$hts.train) != m$no.train )
		stop( paste( "Invalid HTS data for '", m$gene.id, "': number of HTS data and samples differ.", sep="" ))
	
	# dim probes.train == no.train x no.probes
	if ( all( dim( m$probes.train ) != c( m$no.train, m$no.probes )))
		stop( paste( "Invalid probe data for '", m$gene.id ,"': matrix of training probes is of dimension ", dim( m$probes.train ) ," but should be ", c( m$no.train, m$no.probes ),". ", sep="" ))
		
	# dim probes.test == no.test x no.probes
	if ( all( dim( m$probes.test ) != c( m$no.test, m$no.probes )))
		stop( paste( "Invalid probe data for '", m$gene.id ,"': matrix of test probes is of dimension ", dim( m$probes.test ) ," but should be ", c( m$no.test, m$no.probes ),". ", sep="" ))
	
	# make sure that m$hts.test is only one NA (there could be multiple)
	if ( is.na( m$hts.test )) 
		m$hts.test <- NA
	
	# now it's safe to build the object
	object <- new( "TT.Ready.Gene", gene.id=m$gene.id, no.train=m$no.train, no.test=m$no.test, no.probes=m$no.probes, hts.train=m$hts.train, hts.test=m$hts.test, probes.train=m$probes.train, probes.test=m$probes.test )
	return( object )
}

#
# TT.Ready.Gene methods
#

# show()
setMethod( "show", signature( object="TT.Ready.Gene" ), function( object )
	{
		cat( "Gene ID :", object@gene.id, "\n" )
		cat( "#Train  :", object@no.train, "\n" )
		cat( "#Test   :", object@no.test, "\n" )
		cat( "#Probes :", object@no.probes, "\n" )
		cat( "HTS-Test?", all( !is.na( object@hts.test )), "\n" ) # FALSE if at least one sample is missing
	}
)

# show.hts.train()
setGeneric( "show.hts.train", function( object ) standardGeneric( "show.hts.train" ))
setMethod( "show.hts.train", signature( object="TT.Ready.Gene" ), function( object )
	{
		cat( "Training HTS data:\n" )
		cat( object@hts.train, "\n" )
	}
)

# show.hts.test()
setGeneric( "show.hts.test", function( object ) standardGeneric( "show.hts.test" ))
setMethod( "show.hts.test", signature( object="TT.Ready.Gene" ), function( object )
	{
		cat( "Test HTS data:\n" )
		cat( object@hts.test, "\n" )
	}
)

# show.probes.train()
setGeneric( "show.probes.train", function( object ) standardGeneric( "show.probes.train" ))
setMethod( "show.probes.train", signature( object="TT.Ready.Gene"), function( object )
	{
		cat( "Training probe data:\n" )
		cat( object@probes.train, "\n" )
	}
)

# show.probes.test()
setGeneric( "show.probes.test", function( object ) standardGeneric( "show.probes.test" ))
setMethod( "show.probes.test", signature( object="TT.Ready.Gene" ), function( object )
	{
		cat( "Test probe data:\n" )
		cat( object@probes.test, "\n" )
	}
)

#setMethod( "summary", ... )

# run()
setGeneric( "run", function( object, params.object, ... ) standardGeneric( "run" ) )
setMethod( "run", signature( object="TT.Ready.Gene" ), function( object, params.object, gene.tuned=TRUE, tune.quantreg=FALSE, tune.verbose=FALSE, OOB=FALSE )
	{
		if ( !is.null( params.object ) & OOB )
			params.object@OOB <- OOB
		tt.seq <- train.and.predict( object, params.object=params.object, gene.tuned=gene.tuned, tune.quantreg=tune.quantreg, tune.verbose=tune.verbose, OOB=OOB )
		return( tt.seq )	
	}
)

# oob.run()
# alias to run( TT.Ready.Gene, OOB=TRUE )
setGeneric( "oob.run", function( object, params.object, ... ) standardGeneric( "oob.run" ) )
setMethod( "oob.run", signature( object="TT.Ready.Gene" ), function( object, params.object, gene.tuned=TRUE, tune.quantreg=FALSE, tune.verbose=FALSE )
	{
		if ( !is.null( params.object ))
			params.object@OOB <- TRUE
		tt.seq.oob <- train.and.predict( object, params.object=params.object, gene.tuned=gene.tuned, tune.quantreg=tune.quantreg, tune.verbose=tune.verbose, OOB=TRUE )
		return( tt.seq.oob )
	}
)

#-------------------------------------------------------------------------------
#
# TT.Ready.Tx Class
#
#-------------------------------------------------------------------------------
setClass( "TT.Ready.Tx", contains="TT.Ready", representation( tx.id="character", no.txs="numeric" ))

TT.Ready.Tx <- function( m )
{
	# VALIDATIONS
	# no of txs
	if ( length( m$tx.id ) != m$no.txs )
		stop( paste( "Invalid transcript data for '", m$gene.id ,"': number of transcript IDs and transcripts differ.", sep="" ))
	
	# dim of hts.train == no.train x no.txs
	if ( all( dim( m$hts.train) != c( m$no.train, m$no.txs )))
		stop( paste( "Invalid HTS data for '", m$gene.id, "': matrix of transcript HTS data is of dimension ", dim( m$hts.train ) ," but should be ", c( m$no.train, m$no.txs ) ,".", sep="" ))
		
	# dim probes.train == no.train x no.probes
	if ( all( dim( m$probes.train ) != c( m$no.train, m$no.probes )))
		stop( paste( "Invalid probe data for '", m$gene.id ,"': matrix of training probes is of dimension ", dim( m$probes.train ) ," but should be ", c( m$no.train, m$no.probes ),".", sep="" ))
		
	# dim probes.test == no.test x no.probes
	if ( all( dim( m$probes.test ) != c( m$no.test, m$no.probes )))
		stop( paste( "Invalid probe data for '", m$gene.id ,"': matrix of test probes is of dimension ", dim( m$probes.test ) ," but should be ", c( m$no.test, m$no.probes ),".", sep="" ))
	
	# now it's safe to build the object
	object <- new( "TT.Ready.Tx", gene.id=m$gene.id, tx.id=m$tx.id, no.txs=m$no.txs, no.train=m$no.train, no.test=m$no.test, no.probes=m$no.probes, hts.train=m$hts.train, hts.test=m$hts.test, probes.train=m$probes.train, probes.test=m$probes.test )
	return( object )
}

#
# TT.Ready.Tx methods
#

# show()
setMethod( "show", signature( object="TT.Ready.Tx" ), function( object )
	{
		cat( "Gene ID :", object@gene.id, "\n" )
		cat( "Tx ID   :", object@tx.id, "\n" )
		cat( "#Tx     :", object@no.txs, "\n" )
		cat( "#Train  :", object@no.train, "\n" )
		cat( "#Test   :", object@no.test, "\n" )
		cat( "#Probes :", object@no.probes, "\n" )
		cat( "HTS-Test?", all( !is.na( object@hts.test )), "\n" )
	}
)

# show.hts.train()
setMethod( "show.hts.train", signature( object="TT.Ready.Tx" ), function( object )
	{
		cat( "Training HTS data:\n" )
		cat( object@hts.train, "\n" )
	}
)

# show.hts.test()
setMethod( "show.hts.test", signature( object="TT.Ready.Tx" ), function( object )
	{
		cat( "Test HTS data:\n" )
		cat( object@hts.test, "\n" )
	}
)

# show.probes.train()
setMethod( "show.probes.train", signature( object="TT.Ready.Tx"), function( object )
	{
		cat( "Training probe data:\n" )
		cat( object@probes.train, "\n" )
	}
)

# show.probes.test()
setMethod( "show.probes.test", signature( object="TT.Ready.Tx" ), function( object )
	{
		cat( "Test probe data:\n" )
		cat( object@probes.test, "\n" )
	}
)

# run()
setMethod( "run", signature( object="TT.Ready.Tx" ), function( object, params.object, OOB=FALSE )
	{
		if ( is.null( params.object ))
			stop( "params.object must be set when predicting transcript isoform expression; re-run with params.object=TT.Params()" )
		if ( OOB )
			params.object@OOB <- OOB
		tt.seq.txs <- train.and.predict.txs( object, params.object )
		return( tt.seq.txs )	
	}
)

# oob.run()
# alias to run( TT.Ready.Tx, OOB=TRUE )
setMethod( "oob.run", signature( object="TT.Ready.Tx" ), function( object, params.object )
	{
		if ( is.null( params.object ))
			stop( "params.object must be set when predicting transcript isoform expression; re-run with params.object=TT.Params()" )
		params.object@OOB <- TRUE
		tt.seq.oob.txs <- train.and.predict.txs( object, params.object )
		return( tt.seq.oob.txs )
	}
)
