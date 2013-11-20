prepare.data <- function( samples.fn="samples.txt", ma.fn="ma_data.txt", hts.fn="hts_data.txt", g2p.fn="gene_probesets.txt", raw=FALSE, PCs=FALSE )
{
	py.options <- paste( "-s", samples.fn, "-m", ma.fn, "-t", hts.fn, "-g", g2p.fn, sep=" " )

	if ( raw )
	{
		py.options <- paste( py.options, "-r", sep=" " )
	}
	
	if ( PCs )
	{
		py.options <- paste( py.options, "-p", sep=" " )
	}
	
	path <- paste( system.file( package="MaLTE" ), "prepare_data.py", sep="/" )
	cmd <- paste( "python", path, py.options, sep=" " )
	status <- system( cmd, intern=T )
	
	if ( length( status ) == 0 )
		cat( "Successfully prepared gene-level training and test data.\n" )
	else
		cat( "There was a problem running prepare.data(). Please check the error above and retry.\n" )
}

prepare.txs.data <- function( samples.fn="samples.txt", train.fn="train_data.txt", test.fn="test_data.txt", hts.txs.fn="hts_txs_data.txt", g2tx.fn="gene_transcripts.txt" )
{
	py.options <- paste( "-s", samples.fn, "-r", train.fn, "-a", test.fn, "-e", hts.txs.fn, "-t", g2tx.fn, sep=" " )

	path <- paste( system.file( package="MaLTE" ), "prepare_txs_data.py", sep="/" )
	cmd <- paste( "python", path, py.options, sep=" " )
	status <- system( cmd, intern=TRUE, ignore.stdout=TRUE )
	
  if ( length( status ) == 0 )
		cat( "Successfully prepared transcript-level training and test data.\n" )
	else
		cat( "There was a problem running prepare.txs.data(). Please check the error above and retry.\n" )
}

