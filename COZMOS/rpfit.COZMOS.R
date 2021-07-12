rpfit_COZMOS <- function ( OrigData, type, minsplit, maxdepth, pan.con ) {

  row.info.init <- sum( 2^( 0 : maxdepth ) ) + 5
    
  NodeInfo <- matrix ( 0, row.info.init, 10 )
    # 1: "Missclassify", 2: "SumWeight", 3: "Prediction", 
    # 4: "resp1", 5: "resp2", 6: "nodenumber", 7: "NumOfObs", 8: "NumOfObsFinal",
    # 9: "PathNum", 10 :"err.lmt"
    
  SplitInfo <- matrix ( 0, row.info.init, 5 )
    # 1: "nodenumber", 2: "SplitPoint", 3: "VariableNumber", 4: "Count", 5: "NumVSCat"
    
  mc <- if ( any ( type=="c" ) ) 
              max( apply ( as.matrix( OrigData [ , type=="c" ] ), 2, max ) ) else 0
  csplit <- matrix ( 0, 2, mc )
  NodeData <- rep ( 1 , nrow( OrigData ) )
  catsplit <- 0

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Top Down procedure <<<<<<<<<<<<<<<<<<<<<<<<<<<<# 

for ( i in 1 : maxdepth  ) {
  
  if ( max( NodeData ) < 2^( i - 1 ) ) next  
  
  for ( j in ( 2^( i - 1 ) ) : ( 2^i - 1 ) ) {

    #============ Node Information Update ========================#
    GivenData <- GivenData.forsplit <- OrigData[ NodeData==j, ]
    idxY <- ncol( OrigData )
    Given.index <- seq ( 1 : idxY )
    all.n <- length( as.numeric( GivenData ) )

    if ( all.n < idxY ) next
    if ( all.n == idxY ) {
      gdt.y <- GivenData[ length( GivenData ) ]
      n1 <- n2 <- 0
      if ( gdt.y == 1 ) n1 <- 1 else n2 <- 1
      Prediction <- if ( n1 > n2 ) 1 else 2
      weight <- 1
      miss <- 0
    } else {
      gdt.y <- GivenData[, idxY ]
      n1 <- length ( GivenData [ gdt.y == 1, 1 ] )
      n2 <- length ( GivenData [ gdt.y == 2, 1 ] )
      Prediction <- if ( n1 > n2 ) 1 else 2
      weight <- length ( gdt.y)
      miss <- length ( GivenData[ gdt.y != Prediction, 1 ]  )
    }

    NodeInfoTemp <- c( miss, weight, Prediction , n1, n2, j, weight,weight, i-1, 0 )
    row.start <- which( NodeInfo[ , 2 ] == 0) [3]
    NodeInfo [ row.start, ]<- NodeInfoTemp

    #============= Splitting Information Update ====================#
    if ( all.n < ( minsplit * idxY ) ) next
    if ( length ( levels ( factor ( gdt.y ) ) ) == 1  ) next
    if ( min( sum ( gdt.y == 1 ) , sum( gdt.y == 2 ) )  <= ( 0.2 * minsplit )  ) next

    try2dim.condi <-if ( nrow( GivenData ) < 200 ) FALSE else TRUE    
    if ( ( j < 2^( maxdepth - 2 ) ) & ( try2dim.condi == TRUE ) ) {
      var.select.result <- COZMOS.var.selection ( GivenData, type, minsplit ) 
    } else {
      var.select.result <- COZMOS.var.selection.terminal ( GivenData, type ) 
    }

    skip <- var.select.result$skip
    if ( skip==TRUE ) next
    sig.var <- Given.index [ var.select.result$V ]
    sig.var.idx12 <- var.select.result$idx12

    if ( length ( levels ( factor ( GivenData [ , sig.var ] ) ) ) == 1  ) next

    if ( type [ sig.var ] == "n" ) {
      spl.pnt <- if( sig.var.idx12 == 1 ) { 
                        split.L.ftn ( GivenData , sig.var, idxY )[[1]] 
                    } else { var.select.result$spl.pnt }
      idxLeft <- intersect( which( NodeData==j ), 
                            which( OrigData[, sig.var ] < spl.pnt ) )
      idxRight <- intersect( which( NodeData==j ), 
                             which( OrigData[, sig.var ] >= spl.pnt ) )
    } else {
      catsplit <- catsplit + 1
      spl.pnt <- catsplit
      spl.left <- if( sig.var.idx12 == 1 ) {
                        split.L.ftn.categorical ( GivenData , sig.var, idxY )[[1]] 
                    } else { var.select.result$spl.pnt } 
          
      Given.dt.cat <- levels ( factor ( GivenData [ , sig.var ] ) )
      numofcat <-  max( OrigData [ , sig.var ] )
      csplit.temp <- rep ( 1, numofcat )
      csplit.temp [ !is.element( 1 : numofcat, Given.dt.cat ) ] <- 0
      left.cat.orig <- spl.left
      right.cat.orig <- Given.dt.cat [ -match( spl.left, Given.dt.cat ) ]
      if ( any ( csplit.temp == 0 ) ) {
        zeroidx <- which( csplit.temp == 0 )            
        zero.tab <- table( OrigData[ OrigData [ , sig.var ] %in% zeroidx, idxY ] , 
                           OrigData[ OrigData [ , sig.var ] %in% zeroidx, sig.var] )
        if( nrow(zero.tab) >= 2 ){
          zero.cat <- rep( 2, length( zeroidx ) )
          zero.cat [ zero.tab[ 1, ] > zero.tab[ 2, ] ]<- 1

        } else {
          zero.cat <- rep( unique( OrigData[ OrigData [ , sig.var ] %in% zeroidx, idxY ]), length( zeroidx ) )
	      }
        left.cat.orig <- sort( c( left.cat.orig, zeroidx[ zero.cat == 1 ] ) )
        right.cat.orig <- sort( c( right.cat.orig, zeroidx[ zero.cat == 2 ] ) ) 
      }
      csplit.temp [ left.cat.orig ] <- -1
      csplit <- rbind ( csplit , 
                        c( csplit.temp, rep ( 0, mc - length( csplit.temp ) ) ) )
      idxLeft <- intersect( which( NodeData==j ), 
                            which( OrigData [ , sig.var ] %in%  left.cat.orig ) )
      idxRight <- intersect( which( NodeData==j ), 
                             which( OrigData [ , sig.var ] %in%  right.cat.orig ) )
    }
    if ( any( length(idxLeft ) == 0, length(idxRight)==0 ) ) next
    SplitInfoTemp <- c( j, spl.pnt, sig.var, length( gdt.y ),
                            if ( type [ sig.var ] == "n" ) { -1 } else { numofcat } )
    SplitInfo [ which( SplitInfo[ , 1 ] == 0 )[3], ] <- SplitInfoTemp
    NodeData [ idxLeft ] <- j * 2
    NodeData [ idxRight ] <- j * 2 + 1
    }
}
rm(i)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Bottom Up procedure <<<<<<<<<<<<<<<<<<<<<<<<<<<<#
input <- max( SplitInfo[,1] )
if ( input != 0 ) {
    inp <- 0
    bdr.pt <- c( 2^inp, 2^( inp + 1 ) )
    while ( !( input >= bdr.pt[1] & input < bdr.pt[2] ) ) {
      inp <- inp + 1
      bdr.pt <- c( 2^inp, 2^( inp + 1 ) )
    }
    i <- inp + 2
    #============ Final Node Information update ================#
    maxleftnode <- 2^( inp + 1 ) - 1
    NodeInfo <- NodeInfo [ NodeInfo[,6] <= maxleftnode,  ]
      for ( j.f in ( 2^( i - 1 ) ) : ( 2^i - 1 ) ) {
        GivenData.f <- OrigData[ NodeData==j.f , ]
        gdt.nr.f <- length( as.numeric( GivenData.f ) )
        if ( gdt.nr.f < idxY ) next
        if ( gdt.nr.f == idxY ) {
          n1.f <- n2.f <- 0
          if ( GivenData.f[ length( GivenData.f ) ] == 1 ) n1.f <- 1 else n2.f <- 1
          Prediction.f <- if ( n1.f > n2.f ) 1 else 2
          weight.f <- 1
          miss.f <- 0
        } else {
          gdt.y.f <- GivenData.f[, idxY ]
          n1.f <- length ( GivenData.f [ gdt.y.f == 1, 1] )
          n2.f <- length ( GivenData.f [ gdt.y.f == 2, 1] )
          Prediction.f <- if ( n1.f > n2.f ) 1 else 2
          weight.f <- length ( gdt.y.f  )
          miss.f <- length ( GivenData.f [ gdt.y.f != Prediction.f, 1]  )
        }
        NodeInfoTemp.f <- c( miss.f, weight.f, Prediction.f , n1.f, n2.f, j.f , 
                             weight.f, weight.f, i-1, 0 )
        row.restart <- which( NodeInfo[ , 2 ] == 0 ) [3]
        NodeInfo [ row.restart, ]<- NodeInfoTemp.f
      }
    #============ Pruning ================#
	  node.seq <- function ( j, max ) {
		  seq <- j
		  while ( max(seq)< max ){
			  seq <- c( seq,  ( seq * 2 ) , ( seq * 2 ) + 1 )
		  }
		  return ( sort( unique( seq ) ) )
  	}
    pan.con <- 0.11
    row.prune <- ( which( NodeInfo[,2] == 0) [3] - 1 )
    NodeInfo <- NodeInfo [ 1 : row.prune ,]
    SplitInfo <- SplitInfo [ 1 : ( which( SplitInfo[,1] == 0 ) [3] - 1 ) ,]

	  penal.err <- function ( x, data, Ndata, Orig, Ntmp) {
		  j <- x[1]; miss <- x[2]; errmv <- x[3]
      errlv <- Ntmp [ Ntmp[ , 1 ] == j, 2 ]
      Mnode <- max( data[ , 6 ] )
		  errsv <- length( data[data [ , 6 ] %in% node.seq( j, Mnode), 6 ] )
		  xv <- Ndata %in% ( node.seq ( j, Mnode )[-1])
		  NodeDatapred <- data[ match( Ndata, data[ , 6 ] ), 3 ]
		  Y <- Orig[ , ncol( Orig ) ]
		  errv <- sum( NodeDatapred[xv] != Y[xv]  )/ sum( xv )
		  errn <- ncol( Orig )-1;
		  errm <- nrow( Orig )
		  errdelta2 <- 0.95
      return( errv + pan.con * sqrt ( ( log( errn^errlv ) + log( errn^errsv ) + 
                                    log( errm/( 1 - errdelta2 ) ) ) / errmv )  )
	  }

    MNode <- max(NodeInfo[ , 6 ])
    for ( i.back in i : 2 ) {
      NIf.6 <- NodeInfo[ , 6 ]
      y.temp <- NodeInfo[ NIf.6 >= 2^( i.back - 1 ) & NIf.6 < 2^i.back, 6 ]

      for ( y in (y.temp[ y.temp %% 2 == 0 ] / 2 ) ) {
        err.sub <- NodeInfo [ NodeInfo[ , 6 ] == y , 10 ] <- 
            penal.err( x=NodeInfo[ NodeInfo[ , 6 ] == y, c( 6, 1, 7 ) ], data=NodeInfo, 
								 	     Ndata=NodeData, Orig=OrigData, Ntmp=NodeInfo[ ,c( 6, 9 ) ] )
        err.n <- NodeInfo [ NodeInfo[ , 6 ] == y , 7 ]
        err.zero<- NodeInfo [ NodeInfo[ , 6 ] == y , 1 ] / err.n
        nodeprune <- sort( unique( node.seq( y, MNode ) ) )
        if ( err.zero  < err.sub )   {
          NodeData [ NodeData %in% nodeprune[-1] ] <- y
          NTF1 <- !( NodeInfo[ , 6 ] %in% nodeprune[-1])
          NodeInfo <- NodeInfo [ NTF1 , ]
          
          Spl.rm <- !(SplitInfo[ , 1 ] %in% nodeprune )
          if( any( as.logical( SplitInfo[ !Spl.rm & 
                                  !( SplitInfo[ , 5 ] %in% c( 0, -1 ) ), 2 ] ) ) ) {
            csplit <- csplit[ -( SplitInfo[ !Spl.rm & 
                                  !( SplitInfo[ , 5 ] %in% c( 0, -1 ) ), 2 ] + 2), ]
            SplitInfo <- SplitInfo[ Spl.rm, ]
            catsplit <- length( SplitInfo[!( SplitInfo[ , 5] %in% c( 0, -1 ) ), 2 ] )
            SplitInfo[ !( SplitInfo[ , 5] %in% c( 0, -1 ) ), 2 ] <- 1 : catsplit
          } else {
            SplitInfo <- SplitInfo[ Spl.rm, ]
          }
          rm(Spl.rm)
          
        }
      }
    }

    #====== Deleting terminal pairs with the same class======#
    for ( i.back in i : 2 ) {
      NIf.6 <- NodeInfo[ , 6 ]
      y.temp <- NodeInfo[ NIf.6 >= 2^( i.back - 1 ) & NIf.6 < 2^i.back, 6 ]
      for ( y in y.temp [ y.temp %% 2 == 0 ] ) {
        tmn.seq <- c( node.seq( y, max( NIf.6 ) ), node.seq( y+1, max( NIf.6 ) ) )
        terminal <- length ( NodeInfo [ NodeInfo[ , 6 ] %in% tmn.seq, 6 ] )
        if ( terminal == 2  ) {
          if ( NodeInfo[ NodeInfo[ , 6 ]==y, 3 ] == 
               NodeInfo[ NodeInfo[ , 6 ]==y+1, 3 ] ) {
            NodeData [ NodeData==y |  NodeData==( y+1 ) ] <- y/2
            NTF2 <- NodeInfo[ , 6 ]  != y & NodeInfo[ , 6 ]  != ( y+1 )
            NodeInfo <- NodeInfo [ NTF2 , ]

            Spl.rm <- -which( SplitInfo[ , 1 ] == y/2 )
            if( any( as.logical( SplitInfo[ !Spl.rm & 
                                    !( SplitInfo[ , 5 ] %in% c( 0, -1 ) ), 2 ] ) ) ) {
              csplit <- csplit[ -( SplitInfo[ !Spl.rm & 
                                    !( SplitInfo[ , 5 ] %in% c( 0, -1 ) ), 2 ] + 2 ), ]
              SplitInfo <- SplitInfo[ Spl.rm, ]
              catsplit <- length( SplitInfo[!( SplitInfo[ , 5 ] %in% c( 0, -1 ) ), 2 ] )
              SplitInfo[ !( SplitInfo[ , 5 ] %in% c( 0, -1 ) ), 2 ] <- 1 : catsplit
            } else {
              SplitInfo <- SplitInfo[ Spl.rm, ]
            }
            rm(Spl.rm)
            
          }
        }
      }
    }

} else {
  row.end <- ( which( NodeInfo[ , 2 ] == 0 )[3] - 1 )
  NodeInfo <- NodeInfo [ 1 : row.end, ]
  SplitInfo <- SplitInfo [ 1 : row.end, ]
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Result Tables <<<<<<<<<<<<<<<<<<<<<<<<<<<<# 

nrspl <- nrow( SplitInfo )
dsplit <- cbind ( 0, SplitInfo[ , 2 ], 0 )
isplit <- cbind ( SplitInfo[ , 3 ], SplitInfo[ , 4 ], SplitInfo[ , 5 ] )
dnode <- cbind ( NodeInfo[ , 1 ], 0, NodeInfo [ , 2 : 5 ] )
mtemp <- match( NodeInfo[ , 6 ], SplitInfo[ , 1 ] ) - 2
mtemp[ which( is.na( mtemp ) ) ]<- 0
NodeInfo[ mtemp != 0 , 8 ] <- 0
inode <- cbind ( NodeInfo[ , 6 ], mtemp, 1, 0, NodeInfo[ , 7 ], NodeInfo[ , 8 ] )
ncat <- catsplit

dsplit = dsplit[ -( 1 : 2 ), ]
isplit = isplit[ -( 1 : 2 ), ]
csplit = csplit[ -( 1 : 2 ), ]
dnode = dnode[ -( 1 : 2 ), ]
inode = inode[ -( 1 : 2 ), ]

if ( ncat == 1 ) { csplit = t( as.matrix( csplit ) ) }

if ( is.null ( nrow( dsplit ) ) ) { dsplit <- t( as.matrix( dsplit ) ) }
if ( is.null ( nrow( isplit ) ) ) { isplit <- t( as.matrix( isplit ) ) }
if ( is.null ( nrow( csplit ) ) ) { csplit <- t( as.matrix( csplit ) ) }
if ( is.null ( nrow( dnode ) ) ) { dnode <- t( as.matrix( dnode ) ) }
if ( is.null ( nrow( inode ) ) ) { inode <- t( as.matrix( inode ) ) }

return ( list ( which=NodeData, cptable = matrix( 0, 5, 1 ),
        dsplit = dsplit,  isplit = isplit, csplit = csplit,
        dnode = dnode, inode = inode, ncat = ncat ) )
}
