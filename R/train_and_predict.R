#
#
#

#
#
#

train.and.predict <-
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
	
	# test whether there is data for training/testing
	if ( is.na( object@probes.train ) || is.na( object@probes.test ) || mean( object@hts.train ) == 0.00 ) # ????
	{
		m <- list( gene.id=object@gene.id, no.samples=object@no.test, trues=object@hts.test, predictions=NA, cor.S=NA, cor.S.pv=NA, cor.P=NA, cor.P.pv=NA, means=NA, vars=NA, OOB=OOB )
		
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
				D <- D[,order( cor.mat[1,], decreasing=T )[1:(min.probes+1)]]		# take the min.no most correlated
		}
		
		# the model
		model <- cforest( y ~ ., data=D, controls=cforest_control( ntree=ntree, mincriterion=0.95, mtry=mtry ) )
				
		if ( OOB )
		{
			# estimate the prediction using OOB samples
			r.pred <- predict( model, OOB=T )
			r.f <- object@hts.train
		} else
		{
			# use the testing data
			m.f <- object@probes.test
			r.f <- object@hts.test
			D.f <- data.frame( m.f )
			colnames( D.f ) <- paste( "x", 1:object@no.probes, sep="" )
		
			# predict the new response
			r.pred <- predict( model, newdata=D.f )
		}
		
		if ( OOB | length( object@hts.test ) == object@no.test )
		{
			# Spearman correlation
			cor_S <- cor.test( r.f, r.pred,method="spearman", continuity=T)
			cor.S <- cor_S$estimate		# estimate
			cor.S.pv <- cor_S$p.value	# p-value
		
			# Pearson correlation
			cor_P <- cor.test( r.f, r.pred, method="pearson" )
			cor.P <- cor_P$estimate		# estimate
			cor.P.pv <- cor_P$p.value	# p-value
		}
		
		# mean & var
		means <- mean( r.pred, na.rm=T )
		vars <- var( r.pred, na.rm=T )
		
		if ( OOB | length( object@hts.test ) == object@no.test )
		{
			m <- list( gene.id=object@gene.id, no.samples=object@no.test, trues=r.f, predictions=as.vector( r.pred ), cor.S=cor.S, cor.S.pv=cor.S.pv, cor.P=cor.P, cor.P.pv=cor.P.pv, means=means, vars=as.vector( vars ), OOB=OOB )
		} else
		{
			m <- list( gene.id=object@gene.id, no.samples=object@no.test, trues=r.f, predictions=as.vector(r.pred), cor.S=NA, cor.S.pv=NA, cor.P=NA, cor.P.pv=NA, means=means, vars=as.vector( vars ), OOB=OOB )
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
				D <- D[,order( cor.mat[1,], decreasing=T )[1:( min.probes+1 )]]		# take the min.no most correlated
		}
		
		# the model
		model.text <- parse( text=paste( "cforest(", paste( paste( ys, collapse="+" ), sep="" ), "~ ., data=D, controls=cforest_control( ntree=ntree, mincriterion=0.95, mtry=mtry ))", sep="" ))
		model <- eval( model.text )	
				
		if ( OOB )
		{
			# estimate the prediction using OOB samples
			r.pred.1 <- predict( model, OOB=T )
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
		

		if ( OOB | dim( object@hts.test )[1] == object@no.test ){
			# Spearman correlation
			cor.S <- diag( cor( r.f, r.pred, method="spearman" ))
			cor.S.pv <- NA
		
			# Pearson correlation
			cor.P <- diag( cor( r.f, r.pred, method="pearson" ))
			cor.P.pv <- NA
		}

		# mean & var
		means <- colMeans( r.pred, na.rm=T )
		vars <- diag( var( r.pred, na.rm=T ))
		
		if ( OOB | dim( object@hts.test )[1] == object@no.test )
		{
			m <- list( gene.id=object@gene.id, tx.id=object@tx.id, no.txs=object@no.txs, no.samples=object@no.test, trues=r.f, predictions=r.pred, cor.S=cor.S, cor.S.pv=cor.S.pv, cor.P=cor.P, cor.P.pv=cor.P.pv, means=means, vars=as.vector( vars ), OOB=OOB )
		} else
		{
			m <- list( gene.id=object@gene.id, tx.id=object@tx.id, no.txs=object@no.txs, no.samples=object@no.test, trues=r.f, predictions=r.pred, cor.S=NA, cor.S.pv=NA, cor.P=NA, cor.P.pv=NA, means=means, vars=as.vector( vars ), OOB=OOB )
		}
	
		m.obj <- TT.Seq.Tx( m )
		
		return( m.obj )			
	}
}


