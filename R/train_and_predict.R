#
#
#

#
#
#

train.and.predict <-
	function( object,
					params.object=NULL, gene.tuned=TRUE, tune.quantreg=FALSE, tune.verbose=FALSE, OOB=FALSE )
{
	if ( is.null( params.object ) & gene.tuned )
	{
		# gene-specific tuning
		params.object <- tune( object, quantreg=tune.quantreg, verbose=tune.verbose, OOB=OOB )
	} else if ( is.null( params.object ) & !gene.tuned )
	{
		stop( "Missing training parameters: either pass defaults (params.object=TT.Params) OR set gene.tuned=TRUE." )
	}
	
	# get params
	mtry <- params.object@mtry
	ntree <- params.object@ntree
	feature.select <- params.object@feature.select
	min.probes <- params.object@min.probes
	cor.thresh <- params.object@cor.thresh
	OOB <- params.object@OOB
	quantreg <- params.object@quantreg
	
	# test whether there is data for training/testing
	if ( is.na( object@probes.train ) || is.na( object@probes.test ) || mean( object@hts.train ) == 0.00 || length( unique( object@hts.train )) < 10 || object@no.probes < 5 ) 
	# conditions for which to avoid training and predicting
	# 1 - no training probes (how?)
	# 2 - no test probes (how?)
	# 3 - constant zero response
	# 4 - fewer than ten unique training response values (to protect quantregForest)
	# 5 - fewer than five probes (to protect quantregForest)
	{
		m <- list(
			gene.id=object@gene.id,
			no.samples=object@no.test,
			trues=object@hts.test,
			predictions=NA,
			predictions.lower=NA,
			predictions.upper=NA,
			cor.S=NA,
			cor.S.pv=NA,
			cor.P=NA,
			cor.P.pv=NA,
			means=NA,
			vars=NA,
			OOB=OOB
		)
		
		m.obj <- TT.Seq.Gene( m )
		
		return( m.obj )
	} else
	{
		# structure the data
		m <- object@probes.train
		r <- object@hts.train
		D <- data.frame( r, m )
		colnames( D ) <- c( "y", paste( "x", 1:object@no.probes, sep="" ))
		
		# feature selection conditional on the number of probes
		cor.mat <- cor( D ) 					# correlation matrix
		if ( feature.select & object@no.probes > min.probes ) # else use all features
		{
				D <- D[,order( cor.mat[1,], decreasing=TRUE )[1:(min.probes+1)]]		# take the min.no most correlated
		}
		
		# the model
		if ( quantreg )
		{
			model <- quantregForest( y=D[,1], x=D[,2:ncol( D )], mtry=mtry, ntree=ntree )
		} else
		{
			model <- cforest( y ~ ., data=D, controls=cforest_control( ntree=ntree, mincriterion=0.95, mtry=mtry ) )
		}
				
		if ( OOB )
		{
			# estimate the prediction using OOB samples
			if ( quantreg )
			{
				r.pred.all <- predict( model )
				r.pred <- r.pred.all[,2]
			} else
			{
				r.pred <- predict( model, OOB=TRUE )
			}
			r.f <- object@hts.train
		} else
		{
			# use the testing data
			m.f <- object@probes.test
			r.f <- object@hts.test
			
			if ( feature.select & object@no.probes > min.probes )
			{
				D.f <- data.frame( r.f, m.f )			
				colnames( D.f ) <- c( "y", paste( "x", 1:object@no.probes, sep="" ))
				D.f <- D.f[,order( cor.mat[1,], decreasing=TRUE )[1:(min.probes+1)]] # retain only needed probes
				D.f <- D.f[,2:ncol( D.f )]		# get rid of the response
			} else
			{
				D.f <- data.frame( m.f )			
				colnames( D.f ) <- paste( "x", 1:object@no.probes, sep="" )
			}
			
		
			# predict the new response
			if ( quantreg )
			{
				r.pred.all <- predict( model, newdata=D.f )
				r.pred <- r.pred.all[,2]
			} else
			{
				r.pred <- predict( model, newdata=D.f )
			}
		}
		
		if ( OOB | ( length( object@hts.test ) == object@no.test & length( object@hts.test ) > 2 ))
		{
			# Spearman correlation
			cor_S <- cor.test( r.f, r.pred, method="spearman", continuity=TRUE )
			cor.S <- cor_S$estimate   # estimate
			cor.S.pv <- cor_S$p.value	# p-value
		
			# Pearson correlation
			cor_P <- cor.test( r.f, r.pred, method="pearson" )
			cor.P <- cor_P$estimate   # estimate
			cor.P.pv <- cor_P$p.value	# p-value
		}
		
		# mean & var
		means <- mean( r.pred, na.rm=TRUE )
		vars <- var( r.pred, na.rm=TRUE )
		
		if ( OOB | ( length( object@hts.test ) == object@no.test & length( object@hts.test ) > 2 ) )
		{
			m <- list( 
				gene.id=object@gene.id,
				no.samples=object@no.test,
				trues=r.f,
				predictions=as.vector( r.pred ),
				predictions.lower=if ( quantreg ){ as.vector( r.pred.all[,1] ) } else { NA },
				predictions.upper=if ( quantreg ){ as.vector( r.pred.all[,3] ) } else { NA },
				cor.S=cor.S,
				cor.S.pv=cor.S.pv,
				cor.P=cor.P,
				cor.P.pv=cor.P.pv,
				means=means,
				vars=as.vector( vars ),
				OOB=OOB
			)
		} else
		{
			m <- list(
				gene.id=object@gene.id,
				no.samples=object@no.test,
				trues=r.f,
				predictions=as.vector( r.pred ),
				predictions.lower=if ( quantreg ){ as.vector( r.pred.all[,1] ) } else { NA },
				predictions.upper=if ( quantreg ){ as.vector( r.pred.all[,3] ) } else { NA },
				cor.S=NA,
				cor.S.pv=NA,
				cor.P=NA,
				cor.P.pv=NA,
				means=means,
				vars=as.vector( vars ),
				OOB=OOB
			)
		}
				
		m.obj <- TT.Seq.Gene( m )
		
		return( m.obj )
	}
}


train.and.predict.txs <- 
	function( object,
					params.object )
{
	# get params
	mtry <- params.object@mtry
	ntree <- params.object@ntree
	feature.select <- params.object@feature.select
	min.probes <- params.object@min.probes
	cor.thresh <- params.object@cor.thresh
	OOB <- params.object@OOB

	if ( is.na( object@probes.train ) || is.na( object@probes.test ) || mean( object@hts.train ) == 0.00 ) # ?????
	{
		m <- list( gene.id=object@gene.id, tx.id=object@tx.id, no.txs=object@no.txs, no.samples=object@no.test, trues=object@hts.test, predictions=NA, cor.S=NA, cor.S.pv=NA, cor.P=NA, cor.P.pv=NA, means=NA, vars=NA, OOB=OOB )
		
		m.obj <- TT.Seq.Tx( m )
		
		return( m.obj )
	} else
	{
		# structure the data
		m <- object@probes.train
		r <- object@hts.train
		D <- data.frame( r, m )
		ys <- paste( "y", 1:object@no.txs, sep="" )
		xs <- paste("x", 1:object@no.probes, sep="")
		colnames(D) <- c( ys, xs )
		
		# feature selection conditional on the number of probes
		# if there is only one transcript then perform FS otherwise ignore
		cor.mat <- cor( D ) # correlation matrix
		if ( feature.select & object@no.probes > min.probes & object@no.txs == 1 ) # else use all features
		{
				D <- D[,order( cor.mat[1,], decreasing=TRUE )[1:( min.probes+1 )]] # take the min.no most correlated
		}
		
		# the model
		model.text <- parse( text=paste( "cforest(", paste( paste( ys, collapse="+" ), sep="" ), "~ ., data=D, controls=cforest_control( ntree=ntree, mincriterion=0.95, mtry=mtry ))", sep="" ))
		model <- eval( model.text )
				
		if ( OOB )
		{
			# estimate the prediction using OOB samples
			r.pred.1 <- predict( model, OOB=TRUE )
			r.pred <- matrix( as.vector( unlist( r.pred.1 )), byrow=T, nrow=length( r.pred.1 ))
			r.f <- object@hts.train
		} else
		{
			# use the testing data
			m.f <- object@probes.test
			r.f <- object@hts.test
			D.f <- data.frame( m.f )
			colnames(D.f) <- paste( "x", 1:object@no.probes, sep="" )
			
			r.pred.1 <- predict( model, newdata=D.f )
			r.pred <- matrix( as.vector( unlist( r.pred.1 )), byrow=T, nrow=length( r.pred.1 ))
		}
		
		# mean & var
		means <- colMeans( r.pred, na.rm=TRUE )
		vars <- diag( var( r.pred, na.rm=TRUE ))

#		if ( !OOB | length( object@hts.test ) == 1 )
#		{
#			M <- list( gene.id=object@gene.id, tx.id=object@tx.id, no.txs=object@no.txs, no.samples=object@no.test, trues=r.f, predictions=r.pred, cor.S=NA, cor.S.pv=NA, cor.P=NA, cor.P.pv=NA, means=means, vars=as.vector( vars ), OOB=OOB )
#		} else 
		if ( OOB | ( nrow( object@hts.test ) == object@no.test & nrow( object@hts.test ) > 2 ))
		{
			# Spearman correlation
			cor.S <- diag( cor( r.f, r.pred, method="spearman" ))
			cor.S.pv <- NA
		
			# Pearson correlation
			cor.P <- diag( cor( r.f, r.pred, method="pearson" ))
			cor.P.pv <- NA

			M <- list( gene.id=object@gene.id, tx.id=object@tx.id, no.txs=object@no.txs, no.samples=object@no.test, trues=r.f, predictions=r.pred, cor.S=cor.S, cor.S.pv=cor.S.pv, cor.P=cor.P, cor.P.pv=cor.P.pv, means=means, vars=as.vector( vars ), OOB=OOB )
		} else
		{
			M <- list( gene.id=object@gene.id, tx.id=object@tx.id, no.txs=object@no.txs, no.samples=object@no.test, trues=r.f, predictions=r.pred, cor.S=NA, cor.S.pv=NA, cor.P=NA, cor.P.pv=NA, means=means, vars=as.vector( vars ), OOB=OOB )
		}
		
		M.obj <- TT.Seq.Tx( M )
		
		return( M.obj )
	}
}


