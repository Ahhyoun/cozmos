split.L.ftn <- function ( data, v, d ) {

    spsubP <- split.L.ftn.sub ( data, v, d ) 
    split.point <- spsubP$split.point
    split.criterion <- spsubP$split.criterion 
    return ( list( split.point=split.point, split.criterion=split.criterion ) )
}



split.L.ftn.sub <- function ( data, v, d ) {
    if ( !is.matrix( data ) ) return ()
    tempx <- data[ , v ] 
    tempy <- data[ , d ]
    
    n1 <-  length( tempx[ tempy == 1 ] )
    n2 <-  length( tempx[ tempy != 1 ] )

    vx <- sort( unique( as.numeric( tempx ) ) )
    if ( any( n2 < 1, n1 < 1, length(vx) < 2 )) return()
    if ( length( vx ) == 2 ) return(
            list( split.point=mean(vx), 
                  split.criterion= 100 ) )
    bdr <- ( vx [ -1 ] + vx [ - length( vx ) ] ) / 2
    n.bdr <- length ( bdr )

    srt.y <- tempy [ order ( tempx )] ; srt.x <- sort( tempx )
    ac.idx <- which(!duplicated( srt.x ))[-1] - 1

    L10 <- cumsum ( srt.y == 1 )[ ac.idx ] ; R10 <- ( n1 - L10 ) 
    L20 <- cumsum ( srt.y == 2 )[ ac.idx ] ; R20 <- ( n2 - L20 )  

    L1 <- L10 +1
    R1 <- R10 +1
    L2 <- L20 +1
    R2 <- R20 +1
    
    sim <-  ( ( L2 + R1 ) > ( L1 + R2 ) ) * ( 
			        log( pmax( ( R1 * n2 ) / (n1 * L2 ), ( n1 * L2 )/( R1 * n2 )  ) ) / 
			          sqrt( 1/R1 + 1/n2 + 1/n1 + 1/L2 ) ) + 
            ( ( L2 + R1 )<= ( L1 + R2 ) ) * ( 
			        log( pmax( ( L1 * n2 ) / (n1 * R2 ), ( n1 * R2 )/( L1 * n2 )  ) ) / 
			          sqrt( 1/L1 + 1/n2 + 1/n1 + 1/R2 ) )    

	  oddstmp <- L1 * R2 / ( L2 * R1 )

	  odds <- log(  pmax( oddstmp, 1/oddstmp ) ) / sqrt( 1/L1 + 1/L2 + 1/R1 + 1/R2 )
	  
	  listsum2 <- odds-sim
	  
	  sim <- rank( sim )
	  odds <- rank( odds )
  	listsum <- odds - sim 


    mm <- max ( listsum, na.rm=TRUE ) 
    mmid <-na.omit( which( listsum == mm ) ) 

    split.point <- if ( length( mmid ) >= 2 ){ 
                      mm2 <- max ( listsum2, na.rm=TRUE )    
                      mmid2 <-na.omit( which( listsum2 == mm2 ) )
                      if ( length(mmid2) == 1 ){
                          bdr [ mmid2 ]
                      } else{
                          bdr [ mmid2 [ round ( length( mmid2 ) / 2 ) ] ]    
                      }
                   } else {
                            bdr [ mmid ]
                   } 
    split.criterion <- mm     
    return ( list( split.point=split.point, split.criterion=split.criterion ) )
}



split.L.ftn.categorical <- function ( data, v, d ) {
    tempx <- data[ , v ]
    tempy <- data[ , d ]

    n1 <-  length( tempx[ tempy == 1 ] )
    n2 <-  length( tempx[ tempy != 1 ] )
    
    base.tab <- table ( tempy, tempx )
    tabsum <- apply( base.tab, 1, sum )
    if( tabsum[1]> tabsum[2] ){
      base.prob <- base.tab[1,]/colSums( base.tab ) 
    } else {
      base.prob <- base.tab[2,]/colSums( base.tab ) 
    }

    new.x <- match( tempx, colnames( base.tab [ , order ( base.prob, decreasing=TRUE ) ] ) )

    n.bdr <-length( unique( new.x ) )
    
   if ( n.bdr > 2 ) {
    
      srt.y <- tempy [ order ( new.x )] 
      L10 <- cumsum ( srt.y  == 1 ) [ which( !duplicated(sort(new.x)))[-1] -1 ]
      R10 <- n1 - L10
      L20 <- cumsum ( srt.y  == 2 ) [ which( !duplicated(sort(new.x)))[-1] -1 ]
      R20 <- n2 - L20

	    L1 <- L10 + 1 
	    L2 <- L20 + 1 
	    R1 <- R10 + 1
	    R2 <- R20 + 1
	
	    sim <-  ( ( L2 + R1 ) > ( L1 + R2 ) ) * ( 
	      log( pmax( ( R1 * n2 ) / (n1 * L2 ), ( n1 * L2 )/( R1 * n2 )  ) ) / 
	        sqrt( 1/R1 + 1/n2 + 1/n1 + 1/L2 ) ) + 
	      ( ( L2 + R1 )<= ( L1 + R2 ) ) * ( 
	        log( pmax( ( L1 * n2 ) / (n1 * R2 ), ( n1 * R2 )/( L1 * n2 )  ) ) / 
	          sqrt( 1/L1 + 1/n2 + 1/n1 + 1/R2 ) )    
	    
	    oddstmp <- L1 * R2 / ( L2 * R1 )
	    
	    odds <- log(  pmax( oddstmp, 1/oddstmp ) ) / sqrt( 1/L1 + 1/L2 + 1/R1 + 1/R2 )
	    
	    listsum2 <- odds - sim
	    
	    sim <- rank( sim )
	    odds <- rank( odds )
	    listsum <- odds - sim 


      mm <- max ( listsum, na.rm=TRUE )    
      
      mmid <-na.omit( which( listsum == mm ) )
    
      if ( length(mmid) > 1 ) {
        mm2 <- max ( listsum2, na.rm=TRUE )    
        mmid2 <-na.omit( which( listsum2 == mm2 ) )
        if ( length(mmid2) == 1 ){
          split.point <- unique( tempx [ which( new.x <= mmid2 ) ] ) 
        } else{
          split.point <- unique( tempx [ which( new.x <= mmid2 [ round( length( mmid2 ) / 2 ) ] ) ] )
        }
      } else {
        split.point <- unique( tempx [ which( new.x <= mmid ) ] ) 
      } 
        split.criterion <- mm 
        
    } else { 
      split.point <- unique( tempx [ which( new.x == 1 ) ] )
      split.criterion <- -dhyper( sum( new.x == 1 & tempy == 2 ), n2, n1, sum( new.x == 1 ), log=TRUE) 
    }
        
    return ( list( split.point=split.point, split.criterion=split.criterion ) )
}




#########################################################
### following is only for plotting "splitdist" : 
# called at Runsimulsplit_fordistplot_20151126(3).R
#=======================================================

split.L.ftn.sub.for.plot <- function ( data, v, d ) {
	if ( !is.matrix( data ) ) return ()
    tempx <- data[ , v ] 
    tempy <- data[ , d ]
    
    n1 <-  length( tempx[ tempy == 1 ] )
    n2 <-  length( tempx[ tempy != 1 ] )

    
    vx <- sort( unique( as.numeric( tempx ) ) )
    if ( any( n2 < 1, n1 < 1, length(vx) < 2 )) return()
    if ( length( vx ) == 2 ) return(
      list( split.point=mean(vx), 
            split.criterion= 100 ) )
    bdr <- ( vx [ -1 ] + vx [ - length( vx ) ] ) / 2
    n.bdr <- length ( bdr )
    
    srt.y <- tempy [ order ( tempx )] ; srt.x <- sort( tempx )
    ac.idx <- which(!duplicated( srt.x ))[-1] - 1
    
    L10 <- cumsum ( srt.y == 1 )[ ac.idx ] ; R10 <- ( n1 - L10 ) 
    L20 <- cumsum ( srt.y == 2 )[ ac.idx ] ; R20 <- ( n2 - L20 )  
    
    L1 <- L10 +1
    R1 <- R10 +1
    L2 <- L20 +1
    R2 <- R20 +1
    
    sim <-  ( ( L2 + R1 ) > ( L1 + R2 ) ) * ( 
      log( pmax( ( R1 * n2 ) / (n1 * L2 ), ( n1 * L2 )/( R1 * n2 )  ) ) / 
        sqrt( 1/R1 + 1/n2 + 1/n1 + 1/L2 ) ) + 
      ( ( L2 + R1 )<= ( L1 + R2 ) ) * ( 
        log( pmax( ( L1 * n2 ) / (n1 * R2 ), ( n1 * R2 )/( L1 * n2 )  ) ) / 
          sqrt( 1/L1 + 1/n2 + 1/n1 + 1/R2 ) )    
    
    #oddstmp <- L1 * R2 / ( L2 * R1 )
    
    #odds <- log(  pmax( oddstmp, 1/oddstmp ) ) / sqrt( 1/L1 + 1/L2 + 1/R1 + 1/R2 )
    
    
    sim <- rank( sim )
    #odds <- rank( odds )
    listsum <- -sim #odds - sim 
	

    mm <- max ( listsum, na.rm=TRUE ) 
    mmid <-na.omit( which( listsum == mm ) ) 

    split.point <- if ( length( mmid ) >= 2 ){ 
                            bdr [ mmid [ round ( length( mmid ) / 2) ] ]   
                  } else {
                            bdr [ mmid ]
                  } 
    split.criterion <- mm     

    return ( list( split.point=split.point, split.criterion=split.criterion ) )
}






