#
#
#
array2seq <- function( tt.ready, params.object, OOB=FALSE )
{
	tt.seq <- mclapply( tt.ready, run, params.object, OOB, quantreg )
	return( tt.seq )
}

#
#
#
array2seq.oob <- function( tt.ready, params.object )
{
	tt.seq.oob <- mclapply( tt.ready, oob.run, params.object, quantreg )
	return( tt.seq.oob )
}
