
# function for writing formula using 'dsc' file
# - name.dsc: name of 'dsc' file
# - dir.dsc: directory containing 'dsc' file
# e.g. out <- write.formula( 'aba.dsc', 'C:/' )

write.formula <- function( name.dsc, dir.dsc )
{
    # load description file
    
    setwd( dir.dsc )
    dsc <- read.table( name.dsc, skip=3, header=FALSE )
    colnames(dsc) <- c( 'col', 'name', 'type' )
    dsc$name <- as.character( dsc$name )
    
    
    # extract response, numerical & categorical variables
    
    resp <- dsc$name[ dsc$type=='d' ]
    num <- dsc$name[ dsc$type=='n' ]
    cate <- dsc$name[ dsc$type=='c' ]
    
    
    # write down formula: response (class)
    
    fo <- paste( "factor(",resp,") ~", sep="" )
    
    
    # write down formula: numerical variables
    
    if ( length(num) > 0 )
    {
        if ( length(num)==1 )
        {
            fo <- paste( fo, num[1] )
        } else
        {
            for ( i in 1:(length(num)-1) )
            { fo <- paste( fo, num[i], "+" ) }
            fo <- paste( fo, num[length(num)] )
        }
    }
    
    
    # write down formula: categorical variables
    
    if ( length(cate) > 0 )
    {
        if ( length(num) > 0 )
        { fo <- paste( fo, "+" ) }
        if ( length(cate)==1 )
        {
            fo <- paste( fo, " factor(",cate[1],")", sep="" )
        } else
        {
            for ( i in 1:(length(cate)-1) )
            { fo <- paste( fo, " factor(",cate[i],") +", sep="" ) }
            fo <- paste( fo, " factor(",cate[length(cate)],")", sep="" )
        }
    }
    
    
    # return formula
    
    obj <- list( formula=fo, dsc=dsc )
    return( obj )
}
