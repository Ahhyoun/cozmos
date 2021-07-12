pred.for.newdata <- function ( object, newdata ) {
	where <- rep ( 0, nrow ( newdata ) )#; print( dim( newdata))
for ( newdatanum in 1 : nrow ( newdata ) ){

	onedata <- newdata [ newdatanum , ]
	frame.fit <- object$frame
	split.fit <- object$split
	currentnode <- 1

	while ( frame.fit [ as.character( currentnode ), "var" ] != "<leaf>" ) {

	temp.var <- frame.fit [ as.character( currentnode ), "var" ]
	temp.point <- split.fit [ rownames(split.fit) == as.character( temp.var ) &
					  split.fit[ ,"count"]== frame.fit[ rownames( frame.fit ) == currentnode, "n" ], "index" ]
        if ( length ( which( names( attr( newdata,"column.levels") )==temp.var ) )== 1 ) type="c" else type="n"

	if ( type == "n" ) {
		if ( onedata [ which( colnames( newdata )==temp.var ) ] <= temp.point ) {
			currentnode <- currentnode * 2
		} else {
			currentnode <- currentnode * 2 + 1
		}
	} else {
		if ( onedata [ which( colnames( newdata )==temp.var ) ] %in% (which( object$csplit [temp.point, ] == 1 )) ) {
			currentnode <- currentnode * 2
		} else {
			currentnode <- currentnode * 2 + 1
		}

	}

	}
	where [ newdatanum ] <- which ( rownames( frame.fit )== currentnode )
}
return ( where )
}
