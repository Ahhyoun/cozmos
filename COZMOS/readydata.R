
# - name.dsc: name of 'dsc' file
# - dir.dsc: directory containing 'dsc' file

ready.data <- function( dat, name.dsc, dir.dsc )
{
    # load description file
    setwd( dir.dsc )
    dsc <- read.table( name.dsc, skip=3, header=FALSE )
    colnames(dsc) <- c( 'col', 'name', 'type' )
    dsc$name <- as.character( dsc$name )
    dsc$name[dsc$type =='d']='y'

    # Data handling
    data0 <- dat[ , -ncol(dat) ]
    colnames(data0) <- dsc$name[ dsc$col ]

    #turn categorical variables into factor variable in data.frame
    for (j in 1:ncol(data0)) {   
        if (dsc$type[j]=='c') {data0[,j] <- factor(data0[,j])}
        }
    
    # return
    return( data0 )
}
