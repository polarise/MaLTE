#-------------------------------------------------------------------------------
#
# Union of classes numeric, logical and matrix
#
#-------------------------------------------------------------------------------
setClassUnion( 
	name="numericORlogicalORmatrix", 	
	members=c( "numeric", "logical", "matrix" )
) 
# hts.test is either numeric or logical hence the class union

#-------------------------------------------------------------------------------
#
# TT.Seq Abstract Base Class
#
#-------------------------------------------------------------------------------
setClass( 
	Class="TT.Seq", 
	representation( 
		gene.id="character", 
		no.samples="numericORlogicalORmatrix", 
		trues="numericORlogicalORmatrix", 
		predictions="numericORlogicalORmatrix", 
		cor.S="numericORlogicalORmatrix", 
		cor.S.pv="numericORlogicalORmatrix", 
		cor.P="numericORlogicalORmatrix", 
		cor.P.pv="numericORlogicalORmatrix", 
		means="numericORlogicalORmatrix", 
		vars="numericORlogicalORmatrix", 
		OOB="logical" ), 
	contains="VIRTUAL"
)

#-------------------------------------------------------------------------------
#
# TT.Seq.Gene Class
#
#-------------------------------------------------------------------------------
setClass( 
	Class="TT.Seq.Gene", 
	contains="TT.Seq"
)

#
# TT.Seq.Gene constructor
#
TT.Seq.Gene <- function( m )
{
	# validation?
	
	object <- new( 
		Class="TT.Seq.Gene", 
		gene.id=m$gene.id, 
		no.samples=m$no.samples, 
		trues=m$trues, 
		predictions=m$predictions, 
		cor.S=m$cor.S, 
		cor.S.pv=m$cor.S.pv, 
		cor.P=m$cor.P, 
		cor.P.pv=m$cor.P.pv, 
		means=m$means, 
		vars=m$vars, 
		OOB=m$OOB
	)
	return( object )
}

#
# TT.Seq.Gene methods
#

# show()
setMethod( 
	f="show", 
	signature( object="TT.Seq.Gene" ), 
	function( object )
	{
		cat( "Gene ID  :", object@gene.id, "\n" )
		cat( "Mean Expr:", object@means, "\n" )
		cat( "Var  Expr:", object@vars, "\n" )
		cat( "Samples  :", object@no.samples, "\n" )
		cat( "OOB      :", object@OOB, "\n" )
	}
)

# predictions()
setGeneric(
	name="predictions", 
	function( object ) 
		standardGeneric( "predictions" )
)
setMethod( 
	f="predictions", 
	signature( object="TT.Seq.Gene" ), 
	function( object )
	{
		# what do you do when it is NA?
		
		preds <- object@predictions
		return( preds )
	}
)

# gene.id()
setGeneric( 
	name="gene.id", 
	function( object ) 
		standardGeneric( "gene.id" )
)
setMethod( 
	f="gene.id", 
	signature( object="TT.Seq.Gene" ), 
	function( object )
	{
		id <- object@gene.id
		return( id )
	}
)

#-------------------------------------------------------------------------------
#
# TT.Seq.Tx Class
#
#-------------------------------------------------------------------------------
setClass( 
	Class="TT.Seq.Tx", 
	contains="TT.Seq", 
	representation( 
		tx.id="character", 
		no.txs="numericORlogicalORmatrix" )
)

#
# TT.Seq.Tx constructor
#
TT.Seq.Tx <- function( m )
{
	# VALIDATION
	# TODO
	
	object <- new( 
		Class="TT.Seq.Tx",
		gene.id=m$gene.id, 
		tx.id=m$tx.id, 
		no.txs=m$no.txs, 
		no.samples=m$no.samples, 
		trues=m$trues, 
		predictions=m$predictions, 
		cor.S=m$cor.S, 
		cor.S.pv=m$cor.S.pv, 
		cor.P=m$cor.P, 
		cor.P.pv=m$cor.P.pv, 
		means=m$means, 
		vars=m$vars, 
		OOB=m$OOB
	)
	return( object )
}

#
# TT.Seq.Tx methods
#

# show()
setMethod( 
	f="show", 
	signature( object="TT.Seq.Tx" ), 
	function( object )
	{
		cat( "Gene ID  :", object@gene.id, "\n" )
		cat( "Tx ID    :", object@tx.id, "\n" )
		cat( "#Tx      :", object@no.txs, "\n" )
		cat( "Mean Expr:", object@means, "\n" )
		cat( "Var  Expr:", object@vars, "\n" )
		cat( "Samples  :", object@no.samples, "\n" )
		cat( "OOB      :", object@OOB, "\n" )
	}
)

# predictions()
setMethod( 
	f="predictions", 
	signature( object="TT.Seq.Tx" ), 
	function( object )
	{
		# what do you do when it is NA?
	
		preds <- t( object@predictions ) # assume that nothing breaks
		return( preds )
	}
)

# gene.id()
setMethod( 
	f="gene.id", 
	signature( object="TT.Seq.Tx" ), 
	function( object )
	{
		id <- object@gene.id
		return( id )
	}
)

# tx.id()
setGeneric( 
	name="tx.id", 
	function( object ) 
		standardGeneric( "tx.id" )
)
setMethod( 
	f="tx.id", 
	signature( object="TT.Seq.Tx" ), 
	function( object )
	{
		id <- object@tx.id
		return( id )
	}
)


