#
#
#
array2seq <- function( tt.ready, params.object=NULL, ... )
{
	tt.seq <- mclapply( tt.ready, run, params.object=params.object, ... )
	return( tt.seq )
}

#
#
#
array2seq.oob <- function( tt.ready, params.object=NULL, ... )
{
	tt.seq.oob <- mclapply( tt.ready, oob.run, params.object=params.object, ... )
	return( tt.seq.oob )
}
