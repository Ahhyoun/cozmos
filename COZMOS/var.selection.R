
COZMOS.var.selection <- function ( data, var.type, minsplit ) {   
    colnames ( data ) <- rownames ( data ) <- NULL
    vty <- c( "n", "c", "y" ) %in% var.type
    
    five.z.result <- five.sets.of.Z.COZMOS( data, var.type, child=FALSE ) # child=FALSE : (1dim+2dim), child=TRUE : 1dim only

    Z <- five.z.result$z 

    idex.1dim <- c(1, 3)
    idex.2dim <- c(2, 4, 5)
    
    if ( all( is.na( Z[ idex.1dim ] ) ) ) 
            return ( list( V=NULL, idx12 = NULL, spl.pnt = NULL, skip=TRUE ) )   # 1d든 2d 든 더이상 split 불가
    
    Z.1dim <- Z[ idex.1dim ]      
    Z.2dim <- Z[ idex.2dim ]
    
    id.1dim <- which( Z.1dim == max( Z.1dim , na.rm=TRUE ))
    if ( length( id.1dim ) > 1 ) {id.1dim <- id.1dim[ sample( 1:length( id.1dim ), 1 ) ]; print("sample1")}
    
    if ( all( is.na( Z[ idex.2dim ] ) ) ) 
            return ( list( V=five.z.result$cb[[ idex.1dim[ id.1dim ] ]], idx12 = 1, spl.pnt = 0, skip=FALSE ) )   # 1d로 처리
    
    id.2dim <- which( Z.2dim == max( Z.2dim , na.rm=TRUE ))
    if ( length( id.2dim ) > 1 ) {id.2dim <- id.2dim[ sample( 1:length( id.2dim ), 1 ) ]; print("sample2")}

    Z <- Z[ c( idex.1dim[ id.1dim ], idex.2dim[ id.2dim ] )]
    
    CB <- five.z.result$cb
    CB <- list ( CB [[ idex.1dim[ id.1dim ] ]], CB [[ idex.2dim[ id.2dim ] ]] )
    rm( list=c("Z.1dim", "Z.2dim", "idex.1dim", "idex.2dim", "id.1dim", "id.2dim" ))
    
    Z.chisq <- five.z.result$Z.chisq 
    dt.res <- data [ , var.type=="y" ] 
    nx <- ncol ( data [ , var.type!="y" ] ) 
    ndata <- nrow ( data ) 
    CB.dup <- duplicated( c( CB [[ 1 ]], CB [[ 2 ]] ) )
    control <- 0
    
    #qqq <- diff( sort(Z.chisq, decreasing=TRUE) [1:2] ) ; print(qqq)
    #if( ( pnorm( Z[1], lower.tail=FALSE) / ndata  ) <= pnorm( Z[2], lower.tail=FALSE )   ) {
    #print("g")
    #print( c(Z[2], Z[1]))
    if( ( Z[2]-Z[1] ) < (control) ) {
        V <-  CB [[ 1 ]]
        idx12 <- 1 
        spl.pnt <- 0                       
    } else {

       V.candd <- c( CB [[ 1 ]], CB [[ 2 ]] )
       n.candd <- length ( V.candd )
       LL <- RR <- rep ( NA, n.candd )  
       spl.left <- as.list( rep ( 0, n.candd ) )
       
       for ( t in 1 : n.candd ) {
            
           if ( ( is.na( Z.chisq [ V.candd [ t ] ] ) ) | ( CB.dup [[ t ]] == TRUE ) )  next
           
           dt.candd <- data [ , V.candd [ t ] ] ;   
           if ( var.type[ V.candd [ t ] ]== "n" ) {
                spl.left[[t]] <- median(dt.candd) #####$$$$$$$$$$$ 20210114 바꿔봤음 split.L.ftn ( cbind ( dt.candd, dt.res ) , 1, 2 )[[1]]
                idxLeft <- which ( dt.candd < spl.left[[t]] )
                idxRight <- which ( dt.candd >=spl.left[[t]] )
           } else {
               
               #####$$$$$$$$$$$ 이게 최초버전 split.L.ftn.categorical ( data=cbind ( dt.candd, dt.res ), v=1, d=2 )[[1]] : 이게 최초버전임.
               #####$$$$$$$$$$$ 20210114 바꿔봤음 : real simul까지 잘되는 거 확인 but 논리가 약함...  spl.left[[t]] <-  sample( dt.candd, round(length(unique(dt.candd))/2), replace=FALSE) 
                rftbl <- table( dt.res, dt.candd )
                rnkfreq <- rank( apply ( rftbl, 2, sum ), ties.method='random' )
                spl.left[[t]] <-  unique(dt.candd)[ rnkfreq %% 2 != 0 ]
                
                
                Given.dt.cat <- levels ( factor ( dt.candd ) )
                left.cat <- spl.left[[t]]
                right.cat <- Given.dt.cat  [ -match( spl.left[[t]], Given.dt.cat ) ]
                idxLeft <- which( dt.candd %in% left.cat )
                idxRight <- which( dt.candd %in% right.cat )
                rm( list=c("Given.dt.cat","left.cat","right.cat"))
           }
           #print( spl.left[[t]])
           
           Left.L <- data [ idxLeft, ] ; n.obs.L <- length ( idxLeft )
           Right.R <- data [ -idxLeft, ] ; n.obs.R <- ndata - n.obs.L
       
           if ( n.obs.L >= minsplit ){ L.resp <- table( Left.L [ , nx + 1 ] ) } else { L.resp <- 0 }
           if ( n.obs.R >= minsplit ) { R.resp <- table( Right.R [ , nx + 1 ] ) } else { R.resp <- 0 }


           if ( length( L.resp ) <= 1 ) {
                LL [t] <- NA 
           } else { 
                Ltmp <- five.sets.of.Z.COZMOS( Left.L, var.type, child=TRUE )
                LL [t] <- if ( all( is.na( Ltmp$z ) ) ) NA else max( Ltmp$z, na.rm=TRUE ) ;
                rm( Ltmp )
           }
                
           if ( length( R.resp ) <= 1 ) { 
                RR [t] <- NA
           } else { 
                Rtmp <- five.sets.of.Z.COZMOS( Right.R , var.type, child=TRUE )
                RR [t] <- if ( all( is.na( Rtmp$z ) ) ) NA else max( Rtmp$z, na.rm=TRUE ) 
                rm( Rtmp )
           }
          
       rm( list=c("dt.candd", "idxLeft", "idxRight", "Left.L", "Right.R", "n.obs.L", "n.obs.R", "L.resp", "R.resp"))
       }   

       if ( all( is.na ( c( LL[-1], RR[-1] ) ) )) {        
          V <- V.candd [ 1 ]
          idx12 <- 1 
          spl.pnt <- 0       #;  print( paste( "1d: ", V ) )                               
       } else {
          Z.candd <- apply ( cbind( LL, RR ), 1, min, na.rm=TRUE ) 
          Z.candd2 <- apply ( cbind( LL, RR ), 1, max, na.rm=TRUE ) 
          Z.candd3 <- apply ( cbind( LL, RR ), 1, mean, na.rm=TRUE )   #$$$$$$$$$$$ 수정했음 20210114      

          dvd.all <- as.numeric( apply( cbind( LL, RR), 1, function( x ) { ( 2 - sum( is.na(x)) ) } ) ) 
          M.sub <- which ( Z.candd3 [ 2 : 3 ] == max(  Z.candd3 [ 2 : 3 ], na.rm=TRUE ) )

          if ( length( M.sub ) > 1 ) { M.sub <- M.sub[ sample( 1 : 2, 1 ) ] } 
          if ( any(dvd.all<2) |   # dvd.all [ 1 ] < 2 | dvd.all [ M.sub + 1 ] < 2 | #$$$$$$$$$$$$$$$$ 수정했음 20210112
                 (( Z.candd[M.sub+1] - Z.candd2[1] ) < (control) ) ) {
            V <- V.candd [ 1 ]
            idx12 <- 1 
            spl.pnt <- 0  
          } else {    
                  print( "2d" ) 
                  V <- V.candd [ M.sub + 1 ]
                  idx12 <- 2 
                  spl.pnt <- spl.left[[ M.sub + 1 ]]
                  #print( V )
                  #print(spl.pnt)
          }
       }
       rm( list=c("V.candd", "n.candd","Z.candd", "LL","RR" , "dvd.all","M.sub","spl.left" ))
       }   
    return ( list( V=V, idx12 = idx12, spl.pnt = spl.pnt, skip=FALSE ) )
    }

COZMOS.var.selection.terminal <- function ( data, var.type ) {
    colnames ( data ) <- rownames ( data ) <- NULL
    vty <- c( "n", "c", "y" ) %in% var.type

    five.z.result <- five.sets.of.Z.COZMOS( data, var.type, child=TRUE ) 
    Z.chisq <- five.z.result$Z.chisq
    Z <- five.z.result$z
    CB <- five.z.result$cb

    if ( all( is.na( Z.chisq ) ) | all ( is.na( c( Z[1], Z[2] ) ) ) ) { 
        return ( list( V=NULL, idx12 = NULL, spl.pnt = NULL, skip=TRUE ) )
    } else {
        if ( any ( as.logical( lapply( CB, is.null ) )) ) {    # when all vars are num or cat
            not.na.idx <- which( lapply( CB, is.null )==FALSE )
            return ( list( V=as.numeric( CB [[ not.na.idx ]] ), idx12 = 1, spl.pnt = 0 , skip=FALSE) )
        } else {
        
            Z[ is.na ( Z ) ] <- -Inf
            if ( Z[1] >= Z[2] ) {
                V <- as.numeric( CB [[ 1 ]] )
            } else {
                V <- as.numeric( CB [[ 2 ]] )
            }
            return ( list( V=V, idx12 = 1, spl.pnt = 0, skip=FALSE ) )
        }
    }
    
    }

