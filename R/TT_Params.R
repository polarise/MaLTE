#
# TT.Params Class
#
setClass( 
	Class="TT.Params", 
	representation( 
		mtry="numeric", 
		ntree="numeric", 
		feature.select="logical", 
		min.probes="numeric", 
		cor.thresh="numeric", 
		OOB="logical" 
	)
)

#
# TT.Params constructor
#
TT.Params <- function( mtry=2, ntree=1000, feature.select=TRUE, min.probes=15, cor.thresh=0, OOB=FALSE )
{
	# VALIDATION
	# 1 <= mtry <= 10
	if ( mtry < 2 | mtry > 10 )
		stop( "Invalid value for 'mtry'. Should be between 2 and 10 inclusive." )
	
	# 1 <= ntree <= 100000
	if ( ntree < 5 | ntree > 10000 )
		stop( "Invalid value for 'ntree'. Should be between 5 and 10000 inclusive." )
	
	# 1 <= min.probes <= 250
	if ( min.probes < 1 | min.probes > 250 )
		stop( "Invalid value for 'min.probes'. Should be between 1 and 250 inclusive." )
	
	# -1 <= cor.thresh <= 1
	if ( cor.thresh < -1 | cor.thresh > 1 )
		stop( "Invalid value for 'cor.thresh'. Should be between -1 and 1 inclusive." )
	
	# now it's safe to build the object
	object <- new( "TT.Params", mtry=mtry, ntree=ntree, feature.select=feature.select, min.probes=min.probes, cor.thresh=cor.thresh, OOB=OOB )
	
	return( object )
}

#
# TT.Params show method
#
setMethod( 
	f="show", 
	signature( object="TT.Params" ), 
	function( object )
	{
	cat( "mtry           =", object@mtry, "\n" )
	cat( "ntree          =", object@ntree, "\n" )
	cat( "feature.select =", object@feature.select, "\n" )
	cat( "min.probes     =", object@min.probes, "\n" )
	cat( "cor.thresh     =", object@cor.thresh, "\n" )
	cat( "OOB            =", object@OOB, "\n" )
	}
)
