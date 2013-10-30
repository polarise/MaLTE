# start with conditional random forest
# then proceed to quantregForest

.feature.select.tune <- function( min.probes, D, cor.mat, no.probes, r )
{
	if ( no.probes > min.probes ) # else use all features
	{
		D <- D[,order( cor.mat[1,], decreasing=TRUE )[1:(min.probes+1)]]		# take the min.no most correlated
		
		# cforest
		model <- cforest( y ~ ., data=D )
		r.pred <- predict( model, OOB=TRUE )
	
		# Pearson correlation between training and OOB predictions
		cor.P <- cor.test( r, r.pred, method="pearson" )$estimate
	} else
	{
		cor.P <- -1
	}
	cat( min.probes, "=>", cor.P, "\n" )
	cor.P
}

.mtry.tune <- function( mtry, D, r )
{
	# cforest
	model <- cforest( y ~ ., data=D, controls=cforest_control( mtry=mtry ) )
	r.pred <- predict( model, OOB=TRUE )

	# Pearson correlation between training and OOB predictions
	cor.P <- cor.test( r, r.pred, method="pearson" )$estimate

	cat( mtry, "=>", cor.P, "\n" )
	cor.P	
}

.ntree.tune <- function( ntree, D, r, mtry )
{
	# cforest
	model <- cforest( y ~ ., data=D, controls=cforest_control( mtry=mtry, ntree=ntree ) )
	r.pred <- predict( model, OOB=TRUE )

	# Pearson correlation between training and OOB predictions
	cor.P <- cor.test( r, r.pred, method="pearson" )$estimate

	cat( ntree, "=>", cor.P, "\n" )
	cor.P	
}

# should return a TT.Params object
tuned.params <-
	function( object )    # a single gene's data
{
	# structure the data
	m <- object@probes.train
	r <- object@hts.train
	D <- data.frame( r, m )
	colnames( D ) <- c( "y", paste( "x", 1:object@no.probes, sep="" ))
	
	cat( "Your data:\n" )
	print( D )
	
	# min.probe domain
	M <- c( 10, 15, 17, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, object@no.probes )
	M <- sort( M )
	
	# mtry domain
	T <- 2:8
	
	# ntree domain
	N <- c( 50, 100, 250, 500, 750, 1000 )
	
	# tune feature selection
	cat( "Tuning feature selection...\n" )
	cor.mat <- cor( D ) 					# correlation matrix
	cor.P.min.probes <- as.vector( sapply( M, .feature.select.tune, D, cor.mat, object@no.probes, r ))
	if ( all( cor.P.min.probes == -1 ))
		best.min.probes <- object@no.probes
	else
		best.min.probes <- M[ which( cor.P.min.probes == max( cor.P.min.probes )) ]
	cat( "Best min.probes is", best.min.probes, ".\n\n" )
	
	# tune mtry
	cat( "Tuning mtry...\n" )
	D <- D[,order( cor.mat[1,], decreasing=TRUE )[1:(best.min.probes+1)]]		# use the value of best.FS
	cor.P.mtry <- as.vector( sapply( T, .mtry.tune, D, r ))
	best.mtry <- T[ which( cor.P.mtry == max( cor.P.mtry )) ]
	cat( "Best mtry is", best.mtry, ".\n\n" )
	
	# tune ntree
	cat( "Tuning ntree...\n" )
	cor.P.ntree <- as.vector( sapply( N, .ntree.tune, D, r, best.mtry ))
	best.ntree <- N[ which( cor.P.ntree == max( cor.P.ntree )) ]
	cat( "Best ntree is", best.ntree, ".\n\n" )
	
	# create a TT.Params object to return
	best.params <- TT.Params(
									mtry=best.mtry,
									ntree=best.ntree,
									feature.select=if ( best.min.probes < object@no.probes ){ TRUE } else { FALSE },
									min.probes=if( best.min.probes < object@no.probes ){ best.min.probes } else { object@no.probes },
									tune.cor.P=max( cor.P.ntree )
								)
	
	return( best.params )
}

