
Zscore.ftn <- function ( dt.res, cat.fn ) {

  ctest.fn <- chisq.test( as.matrix( table( dt.res , cat.fn ) ) )

	X.fn <- ctest.fn $ statistic
	df.fn <- ctest.fn $ parameter
	p.fn <- ctest.fn $ p.value

  #if ( p.fn > 1e-10 ) {                              ###########20200101 수정함
  #      Z.fn <- qnorm( 1 - p.fn )                   ###########20200101 수정함
  #} else {                                           ###########20200101 수정함
	   if ( df.fn == 1 ) { 
            Z.fn <- sqrt( X.fn )
	   } else {
            W.fn <- X.fn - df.fn + 1
		    cal.fn <- ( df.fn - 1 ) * log ( ( df.fn - 1 ) / X.fn ) + W.fn
            Z.fn <- 1 / abs( W.fn ) * ( W.fn - 1/3 - 0.08 / df.fn )  * sqrt ( cal.fn )
            Z.fn <- abs( Z.fn ) ## 20210114  수정함 $$$$$$$$$$$$$$$$$$$$$$$
	   }
  #}
	
# 이거두줄 지금 테스트중중	
#	Z.fn <- 1-p.fn     ## 20210107  수정함 $$$$$$$$$$$$$$$$$$$$$$$ 확실히 더 나으면 위에꺼도 다 지우자!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#	if( Z.fn == 0 ) Z.fn <- runif(1) ## 20210108  수정함 $$$$$$$$$$$$$$$$$$$$$$$ 분할표 모든 행 빈도 다 똑같이 나오는 경우 처리하기 위한 것.

	
	
	  #if ( is.infinite( Z.fn ) ) { Z.fn <- NA }
	      #if ( is.infinite( Z.fn ) ) {print(table( dt.res , cat.fn )); print(Z.fn)}
	      #if ( is.nan( Z.fn ) ) {print(table( dt.res , cat.fn )); print(Z.fn)}
	      #if ( is.na( Z.fn ) ){print(table( dt.res , cat.fn )); print(Z.fn)} 

  return( Z.fn )
}
cutbyqtl <- function ( dt, n ){
    
    psq <- seq( 0, 1, length = n + 1 )
    cpt <- quantile( dt, prob = psq, type=5 )    
    if ( length( unique( cpt ) ) < (n+1) ) { cpt <- cpt + 1e-4 * sort( abs(rnorm( n+1 )) )*duplicated(cpt) } 
    
    #f.vec <- cut( dt, breaks = cpt, labels = 1 : n, include.lowest = TRUE )
    f.vec <- cut( dt, breaks = cpt, labels = 1 : n , include.lowest = TRUE ) ; 
    #if ( any(is.na(f.vec)) ) print(f.vec)
    #f.vec[is.na(f.vec)] <- n
    return( as.factor( as.numeric( f.vec ) ) )
    }
 
    
exp.table <- function ( res, fac ) {
    tab <- as.matrix( table( res , fac ) )
    exp.tab <- ( rowSums( tab ) %o% colSums( tab ) ) / sum( tab )
    return ( exp.tab )
    }

    
five.sets.of.Z.COZMOS <- function ( data, var.type, child=FALSE ) {

	colnames ( data ) <- rownames ( data ) <- NULL
  data.res <- data [ , var.type=="y" ]
	nr <- nrow ( data )
	vty <- c( "n", "c", "y" ) %in% var.type
  colNum <- 2

	var.value.zn <- var.value.znn <- var.value.zc <- var.value.zcc <- var.value.znc <- 
	cb.zn <- cb.znn <-cb.zc <- cb.zcc <- cb.znc <- NULL
  Z.chisq.zn<- Z.chisq.zc <- NA
    
  idx.n <- which( var.type == "n" ) 
  idx.c <- which( var.type == "c" )
  
  #==========================================================================
  #==========================================================================
  
	#==========  marginal test for each numerical variable ==========#
	if ( vty [1] == TRUE ) {
        data.num <- as.matrix( data [ , idx.n ] )
        n.zn <- ncol( data.num )
        Z.chisq.zn <- rep ( NA, n.zn )
        for ( t in 1 : n.zn ) {
            dt.nm <- data.num [ , t ]
            Minlevel <- length( unique( dt.nm ) )
            if ( Minlevel <= 1 ) next
              
  			      #=======================================
              start.grp <- 7
              tmp.zn <- cutbyqtl ( dt.nm , start.grp )   
              pt.tmp <-length ( levels( tmp.zn ) )  
              
              #if( pt.tmp <= 5 ) {                                      ###########20200101 수정함
              #  tmp.zn <- cutbyqtl ( dt.nm , 5 )                       ###########20200101 수정함
              #  pt.tmp <-length ( levels( tmp.zn ) )                   ###########20200101 수정함                 
              #}                                                        ###########20200101 수정함

              #  while ( any( exp.table( data.res, tmp.zn ) < 5 ) & pt.tmp > colNum  ) {                   ###########20200101 수정함
  		      	while ( any( exp.table( data.res, tmp.zn ) < 5 ) & pt.tmp > colNum  ) {
                        pt.tmp <- max( pt.tmp - 1, colNum )
                        tmp.zn <- cutbyqtl ( dt.nm, pt.tmp )
                        pt.tmp <- length ( levels( tmp.zn ) )
                } 
		          if ( pt.tmp  == 1 ) { tmp.zn <- as.numeric( dt.nm < median( dt.nm ) ) }     ##### VVVVVVVV 이런경우 있나 출력해보기!!!!! 20210108  
			        #=======================================
            
            rm( list = c( "start.grp","pt.tmp" ) )
            Z.chisq.zn [ t ] <- Zscore.ftn ( data.res, tmp.zn )
            
            rm( list = c( "dt.nm", "Minlevel", "tmp.zn" ) )        
        }
        rm( t )

        if ( all( is.na( Z.chisq.zn ) )  ) {
            var.value.zn <- cb.zn <- NA
        } else {
            var.value.zn <- max( Z.chisq.zn, na.rm=TRUE )
            cb.zn <- idx.n [ na.omit( which( Z.chisq.zn == var.value.zn ) ) ]
            if ( length( cb.zn ) > 1 ){  
               cb.zn <- cb.zn [ sample( 1 : length( cb.zn ) , 1 )  ]    ##### VVVVVVVV 이런경우 있나 출력해보기!!!!! 20210108
            } 
            
        }
	}
	
	#==========  marginal test for each categorical variable ==========#
	if ( vty[2] == TRUE ) {
		    data.cat <- as.matrix( data [ , idx.c ] )
		    n.zc <- ncol( data.cat )
		    Z.chisq.zc <- rep ( NA, n.zc )
       
        for ( b in 1 : n.zc ) {
            dt.ct <- data.cat[ , b ]
            Minlevel.ct <- length( unique( dt.ct ) )
            if ( Minlevel.ct <= 1 ) next           
	          tmp.zc <- dt.ct
            Z.chisq.zc [ b ] <- Zscore.ftn ( data.res, tmp.zc )
            rm( list = c( "dt.ct", "Minlevel.ct", "tmp.zc" ) )   
        }	
        rm( b )

        
        if ( all( is.na( Z.chisq.zc ) )  ) {
            var.value.zc <- cb.zc <- NA
        } else {
            var.value.zc <- max( Z.chisq.zc, na.rm=TRUE  )
            cb.zc <- idx.c [ na.omit( which( Z.chisq.zc == var.value.zc ) ) ]
            if ( length( cb.zc ) > 1 ) { 
                    cb.zc <- cb.zc [ sample( 1 : length( cb.zc ) , 1 ) ] 
            }
        }
  }
        
  #==========================  2 dimension ====================================
  #============================================================================
        
  #==========  interaction test for each pair of numerical variable ==========#
  if( vty [1] == TRUE ){
    if( child==FALSE & n.zn >= 2 ) {
   
		  cb.temp <- combn( n.zn, 2 )
		  Z.chisq.znn <- rep ( NA , ncol( cb.temp ) )

		  for ( a in 1 : ncol( cb.temp ) ) { 
            cb <- cb.temp [ , a ] 	
            dt.nm1 <- data.num [ , cb[ 1 ] ]
            dt.nm2 <- data.num [ , cb[ 2 ] ]
            Minlevel.1 <- length ( unique( dt.nm1 ) )
            Minlevel.2 <- length ( unique( dt.nm2 ) )
            if (  Minlevel.1 < 3  |  Minlevel.2 < 3 ) next ##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ bian simul에서는 별 차이 없었음 20210108
            

            Minlevel <- c( Minlevel.1, Minlevel.2 )

            f.nm1 <- cutbyqtl ( dt.nm1, 3 ) 
            f.nm2 <- cutbyqtl ( dt.nm2, 3) 
            if (  length ( unique( f.nm1 ) ) < 2 ) { f.nm1  <- as.numeric( dt.nm1 < median( dt.nm1 ) ) }   #### _3 이거 수정했음 #next   
            if (  length ( unique( f.nm2 ) ) < 2 ) { f.nm2  <- as.numeric( dt.nm2 < median( dt.nm2 ) ) }   #### _3 이거 수정했음 #next   
	    
            
            tmp.znn <- rep ( 0, nr )
            for ( nn1 in 1 : max( as.numeric( f.nm1 ) ) ){
                for ( nn2 in 1 : max( as.numeric( f.nm2 ) ) ) {
			            tmp.znn[ which ( f.nm1 == nn1 & f.nm2 == nn2 ) ] <- ( nn1 - 1 ) * max( as.numeric( f.nm2 ) ) + nn2
                }
            }
            rm( list = c( "f.nm1", "f.nm2" ) ) 
            
            Z.chisq.znn [ a ] <-  Zscore.ftn ( data.res, tmp.znn )  
            rm ( list = c( "cb", "dt.nm1", "dt.nm2", "Minlevel.1","Minlevel.2","Minlevel",
                           "tmp.znn" ) ) 
		  }
      rm ( a )
		  
		  if ( all( is.na( Z.chisq.znn ) ) ) {
          var.value.znn <- cb.znn <- NA
      } else {
		      var.value.znn <- max( Z.chisq.znn, na.rm=TRUE )
		      cb.znn.col <- na.omit( which( Z.chisq.znn == var.value.znn ))
		      cb.znn <- rbind( idx.n [ cb.temp [ 1, cb.znn.col ] ], idx.n [ cb.temp [ 2, cb.znn.col ] ] )
		      if ( length( cb.znn.col ) > 1 ) {
              cb.znn <- cb.znn [ , sample( 1 : length( cb.znn.col ) , 1 )  ] ##$$$$$$$$$$$$$$$ vector ??? $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ bian simul에서는 별 차이 없었음 20210108 
          }
          rm( list = c( "cb.temp" ) )
      }	
    }
  }
  #==========================================================================

	
	#==========  interaction test for each pair of categorical variable ==========#
  if( vty[2] == TRUE ){
    if( child==FALSE & n.zc >= 2 ) { 
    
      cb.temp.zcc <- combn( n.zc, 2 )
      nc.cb.zcc <- ncol( cb.temp.zcc )
      Z.chisq.zcc <- rep ( NA ,nc.cb.zcc )

		  for ( y in 1 : nc.cb.zcc ) {

			  cb.zcc <- cb.temp.zcc [ , y ] 	
			  dt.ct1 <- data.cat [ , cb.zcc[1] ]
			  dt.ct2 <- data.cat [ , cb.zcc[2] ]
        Minlevel.ct1 <- length( unique( dt.ct1 ) )
        Minlevel.ct2 <- length( unique( dt.ct2 ) )
            	
        if ( any( Minlevel.ct1 <= 1, Minlevel.ct2 <= 1 ) ) next  
			  #dt.ct1 <- table.modify.catcat ( dt.ct1, data.res ) 
			  #dt.ct2 <- table.modify.catcat ( dt.ct2, data.res ) 

			  tmp.zcc <- rep ( 0, nr )

        for ( q in 1 : length( unique( dt.ct1 ) ) ){
          for ( w in 1 : length( unique( dt.ct2 ) ) ) {
					  tmp.zcc [ which ( dt.ct1 == unique( dt.ct1 )[q] & dt.ct2 == unique( dt.ct2 )[w] ) ] <- paste( q, w )
          }
        }          
        Z.chisq.zcc [ y ] <- Zscore.ftn ( data.res, tmp.zcc )
			  rm( list = c( "cb.zcc", "dt.ct1", "dt.ct2", "Minlevel.ct1", "Minlevel.ct2", "tmp.zcc" ) )   
		  }
      rm( y )
      
      if ( all ( is.na( Z.chisq.zcc ) ) ){
          var.value.zcc <- cb.zcc <- NA
      } else {
          var.value.zcc <- max( Z.chisq.zcc, na.rm=TRUE ) 
		      cb.zcc.col <- na.omit( which( Z.chisq.zcc == var.value.zcc ) )
		      cb.zcc <- rbind( idx.c[ cb.temp.zcc [ 1, cb.zcc.col ] ], idx.c[ cb.temp.zcc [ 2, cb.zcc.col ] ])   
		      if ( length( cb.zcc.col ) > 1 ) {
		        cb.zcc <- cb.zcc [ , sample ( seq ( 1 : length( cb.zcc.col ) ), 1 ) ] ##$$$$$$$$$$$$$$$ vector ??? $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		      }
      }	
      rm( list = c( "cb.temp.zcc" ) )
    }
  }
        
	#==========================================================================
	#==========================================================================
		
	#==========  interaction test for each pair of numerical and categorical variable ==========#    
	if( child==FALSE & all( vty[ 1 : 2 ] ) ) {
	  
		cb.temp0 <- c()
		for ( p in 1 : n.zc ) { cb.temp0 <- c( cb.temp0, rep( p, n.zn ) ) }; rm( p )
		cb.temp.znc <- rbind( rep ( 1 : n.zn, n.zc ), cb.temp0 ); rm( cb.temp0 )
    nc.cb.znc <- ncol( cb.temp.znc )
    Z.chisq.znc <- rep ( NA , nc.cb.znc  )

		for ( r in 1 :  nc.cb.znc  ) {
			cb.znc <- cb.temp.znc [ , r ]
			dt.nm.nc <- data.num [ , cb.znc[1] ]
			dt.ct.nc <- data.cat [ , cb.znc[2] ] 
      Minlevel.nm.nc <- length( unique( dt.nm.nc ) )
      Minlevel.ct.nc <- length( unique( dt.ct.nc ) )
            
      if ( any( Minlevel.nm.nc < 3, Minlevel.ct.nc <= 1 ) ) next     ##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
			tmp.znc <- rep ( 0, nr )
      
      tmp.znc.n <- cutbyqtl ( dt.nm.nc, 3 )
      #tmp.znc.c <- table.modify.catcat ( dt.ct.nc, data.res )
      tmp.znc.c <- dt.ct.nc
			if ( length( unique( tmp.znc.n ) ) < 2 ) { tmp.znc.n <- as.numeric( dt.nm.nc < median( dt.nm.nc ) ) }   #next  #### _3 이거 수정했음
      
      for ( e in 1 : length( unique( tmp.znc.n ) ) ){
          for ( u in 1 : length( unique( tmp.znc.c ) ) ) {
			       tmp.znc [ which ( tmp.znc.n == unique( tmp.znc.n )[e] & tmp.znc.c == unique( tmp.znc.c )[u] ) ] <- paste( e, u )
          }
      }
      Z.chisq.znc [ r ] <- Zscore.ftn ( data.res, tmp.znc )        
      rm( list = c( "cb.znc", "dt.nm.nc", "dt.ct.nc", "Minlevel.nm.nc", "Minlevel.ct.nc",
                    "tmp.znc", "tmp.znc.n", "tmp.znc.c") ) 
		}
    rm( r ); rm( nc.cb.znc )
        
    if ( all ( is.na( Z.chisq.znc ) ) ){
      var.value.znc <- cb.znc <- NA
    } else {
		  var.value.znc <- max( Z.chisq.znc, na.rm=TRUE )
		  cb.znc.col <- na.omit( which( Z.chisq.znc == var.value.znc ) )
		  cb.znc <- rbind( idx.n [ cb.temp.znc [ 1, cb.znc.col ] ] , idx.c [ cb.temp.znc [ 2, cb.znc.col ] ] )
		  if ( length( cb.znc.col ) > 1 ) {
		    cb.znc <- cb.znc [ , sample ( seq ( 1 : length( cb.znc.col ) ), 1 ) ]   ##$$$$$$$$$$$$$$$ vector ??? $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		  }
      rm( cb.znc.col )
		}
    rm ( list = c( "cb.temp.znc") )	
	}


if( child== FALSE ) {
        Z.chisq <- c()
        Z.chisq [ idx.n ] <- Z.chisq.zn
        Z.chisq [ idx.c ] <- Z.chisq.zc
        if ( is.null( var.value.zn ) ) var.value.zn <- NA
        if ( is.null( var.value.znn ) ) var.value.znn <- NA
        if ( is.null( var.value.zc ) ) var.value.zc <- NA
        if ( is.null( var.value.zcc ) ) var.value.zcc <- NA
        if ( is.null( var.value.znc ) ) var.value.znc <- NA
         
        result.list <- list ( z = c( var.value.zn, var.value.znn, var.value.zc, var.value.zcc, var.value.znc ),
                    cb = list( cb.zn, cb.znn, cb.zc, cb.zcc, cb.znc ),
                    Z.chisq = Z.chisq )
        return ( result.list )
    } else {
        Z.chisq <- c()
        Z.chisq [ idx.n ] <- Z.chisq.zn
        Z.chisq [ idx.c ] <- Z.chisq.zc
        if ( is.null( var.value.zn ) ) var.value.zn <- NA
        if ( is.null( var.value.zc ) ) var.value.zc <- NA
        return ( list ( z = c( var.value.zn, var.value.zc ),
                    cb = list( cb.zn, cb.zc ),
                    Z.chisq = Z.chisq )  )
    }

}

