rm( list=ls() )

# load R files
dir.R <-'E:\\COZMOS'   # directory containing codes and data
setwd( dir.R )

source( 'readydata.R' )
source( 'functions.R' )
source( 'write.formula.R' )
source( 'five.z.R' )
source( 'var.selection.R')   
source( 'split.select.R' )
source( 'rpfit.COZMOS.R' )  
source( 'COZMOS.R' ) 
source( 'pred.for.newdata.R' )
source( 'predict.COZMOS.R' )


# Preparing Data 

dat <- read.table( 'hea.dat', header=FALSE )
data0 <- dat[ , -ncol(dat) ]
out.dsc <- write.formula( 'hea.dsc', dir.R )
out.dsc$dsc$name [ out.dsc$dsc$type =='d' ] = 'y'  ## change the class var name to 'y'
colnames( data0 ) <- out.dsc$dsc$name[ out.dsc$dsc$col ]

for (j in 1 : ncol( data0 ) ) {  
    if ( out.dsc$dsc$type[j] == 'c' ) { data0[,j] <- factor( data0[ ,j ] ) }
    }
data0 <-na.omit( data0 ) 
nr <- nrow( data0 )
nc <- ncol( data0 )

type.num <- as.numeric( lapply( data0, is.factor ) )
type <- rep ( "n", length ( type.num ) )
type [ which( type.num==1 ) ] <- "c" 
type [ length ( type.num ) ] <- "y"

tr.ratio <- 0.7
tr.n <- round( nr * tr.ratio )
tr.idx <- sample( 1:nr, tr.n, replace=FALSE )
train <- data0[ tr.idx, ]
testX <- data0[ -tr.idx, -nc]
testY <- data0[ -tr.idx, nc ]


# COZMOS tree fitting

cozmos0fit <- COZMOS ( factor(y)~., data=train, minsplit=10, maxdepth=5, type=type )
cozmos0pred <- predict.COZMOS ( cozmos0fit , newdata=testX )

print( cozmos0fit )
print( table( testY, cozmos0pred ) )

setwd( dir.R )

	