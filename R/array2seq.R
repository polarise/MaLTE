#
#
#
array2seq <- function( tt.ready, params.object=NULL, gene.tuned=TRUE, tune.quantreg=FALSE, tune.verbose=FALSE, OOB=FALSE )
{
	tt.seq <- mclapply( tt.ready, run, params.object=NULL, gene.tuned=TRUE, tune.quantreg=FALSE, tune.verbose=FALSE, OOB=OOB )
	return( tt.seq )
}

#
#
#
array2seq.oob <- function( tt.ready, params.object=NULL, gene.tuned=TRUE, tune.quantreg=FALSE, tune.verbose=FALSE )
{
	tt.seq.oob <- mclapply( tt.ready, oob.run, params.object=NULL, gene.tuned=TRUE, tune.quantreg=FALSE, tune.verbose=FALSE )
	return( tt.seq.oob )
}
