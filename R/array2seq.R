#
#
#
array2seq <- function( tt.ready, params.object, ... )
{
	tt.seq <- mclapply( tt.ready, run, params.object, ... )
	return( tt.seq )
}

#
#
#
array2seq.oob <- function( tt.ready, params.object, ... )
{
	tt.seq.oob <- mclapply( tt.ready, oob.run, params.object, ... )
	return( tt.seq.oob )
}
