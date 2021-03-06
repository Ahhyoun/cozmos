#  I revised 'rpart.R' in rpart package.

COZMOS <- function(formula, data, weights, subset,
		   na.action=na.omit, method, model=FALSE, x=FALSE, y=TRUE,
		   parms, control, cost, minsplit, maxdepth, type, ...)
{
    Call <- match.call()
    if (is.data.frame(model)) {
	m <- model
	model <- FALSE
	}
    else {
	m <- match.call(expand.dots=FALSE)
	m$model <- m$method <- m$control<- NULL
	m$x <- m$y <- m$parms <- m$... <- NULL
	m$cost <- NULL
	m$minsplit <- m$maxdepth <- m$type <- NULL
	m$na.action <- na.action
	m[[1L]] <- as.name("model.frame")
	m <- eval(m, parent.frame())
	}
    Terms <- attr(m, "terms")
    if(any(attr(Terms, "order") > 1L))
	stop("Trees cannot handle interaction terms")

    Y <- model.response(m) 
    wt <- model.weights(m)
    if(length(wt)==0L) wt <- rep(1.0, nrow(m))
    offset <- attr(Terms, "offset")
    X <- rpart.matrix(m) 
    nobs <- nrow(X)
    nvar <- ncol(X)

    if (missing(method)) {
	if (is.factor(Y) || is.character(Y))      method <- "class"
        else if (inherits(Y, "Surv")) method <- "exp"
	else if (is.matrix(Y)) method <- "poisson"
	else                   method <- "anova"
    }

    if (is.list(method)) {
	## User written split methods
	mlist <- method
	method <- "user"

	init <- if (missing(parms)) mlist$init(Y, offset, wt=wt)
	else  mlist$init(Y, offset, parms, wt)

	method.int <- 4L     #the fourth entry in func_table.h

        ## assign this to avoid garbage collection
        keep <- rpartcallback(mlist, nobs, init)
    } else {
	method.int <- pmatch(method, c("anova", "poisson", "class", "exp"))
	if (is.na(method.int)) stop("Invalid method")
	method <- c("anova", "poisson", "class", "exp")[method.int]

	if (method.int == 4L) method.int <- 2L

	 init <-if (missing(parms))
             (get(paste("rpart", method, sep='.')))(Y, offset, ,wt)
	else
            (get(paste("rpart", method, sep='.')))(Y, offset, parms, wt)

        ns <- asNamespace("rpart")
        if(!is.null(init$print)) environment(init$print) <- ns
        if(!is.null(init$summary)) environment(init$summary) <- ns
        if(!is.null(init$text)) environment(init$text) <- ns
    
    }
    
    Y <- init$y
    xlevels <- attr(X, "column.levels")   #-> NULL
      
    cats <- rep(0L, ncol(X))
    if(!is.null(xlevels)) {
	cats[match(names(xlevels), dimnames(X)[[2L]])] <-
            unlist(lapply(xlevels, length))
    }

    ## We want to pass any ... args to rpart.control, but not pass things
    ##  like "dats=mydata" where someone just made a typo.  The use of ...
    ##  is just to allow things like "cp=.05" with easier typing
    extraArgs <- list(...)
    if (length(extraArgs)) {
	controlargs <- names(formals(rpart.control))  #legal arg names
	indx <- match(names(extraArgs), controlargs, nomatch=0L)
	if (any(indx==0L))
            stop(gettextf("Argument %s not matched", names(extraArgs)[indx==0L]),
                 domain = NA)
    }

    controls <- rpart.control(...)
    if (!missing(control)) controls[names(control)] <- control
    
    xval <- controls$xval
    if (is.null(xval) || (length(xval)==1L && xval==0) || method=="user") {
	xgroups <-0
	xval <- 0
    }
    else if (length(xval)==1L) {
	# make random groups
        xgroups <- sample(rep(1L:xval, length=nobs), nobs, replace=FALSE)
    }
    else if (length(xval) == nobs) {
	xgroups <- xval
	xval <- length(unique(xgroups))
    }
    else {
	# Check to see if observations were removed due to missing
	if (!is.null(attr(m, "na.action"))) {
	    # if na.rpart was used, then na.action will be a vector
	    temp <- as.integer(attr(m, "na.action"))
	    xval <- xval[-temp]
	    if (length(xval) == nobs) {
		xgroups <- xval
		xval <- length(unique(xgroups))
		}
	    else stop("Wrong length for 'xval'")
        }
	else stop("Wrong length for 'xval'")
    }

    
    ##
    ## Incorporate costs
    ##
    if (missing(cost)) cost <- rep(1.0, nvar)
    else {
	if (length(cost) != nvar)
            stop("Cost vector is the wrong length")
	if (any(cost <=0)) stop("Cost vector must be positive")
    }
    
    ##
    ## Have s_to_rp consider ordered categories as continuous
    ##  A right-hand side variable that is a matrix forms a special case
    ## for the code.
    ##
    #dyn.load("c:/Ctemp/s_to_rp.dll")
    tfun <- function(x) {
	if (is.matrix(x)) rep(is.ordered(x), ncol(x))
	else is.ordered(x)
    }
    labs <- sub("^`(.*)`$", "\\1", attr(Terms, 'term.labels'))
    isord <- unlist(lapply(m[labs], tfun))


    #===================================================================
	dep <- init$y
  	dep [ dep==0 ] <- 2 
	OrigData <- cbind( as.matrix(X), dep )
    	rp <- rpfit_COZMOS( OrigData, type, minsplit, maxdepth, pan.con  )
    #===================================================================

    # rpfit$newX[1:n] contains the final sorted order of the observations
    nodes <- nrow( rp$dnode )          # total number of nodes
    nsplit<- nrow( rp$dsplit )      # total number of splits, primary and surrogate
    numcp <- 1     # number of lines in cp table
    ncat  <- rp$ncat    # total number of categorical splits
    numresp<- 2 + 1    # length of the response vector

    if (nsplit == 0) xval <- 0
    cpcol <- if (xval>0 && nsplit>0) 5L else 3L
    if (ncat == 0) catmat <- 0
    else         catmat <- matrix(integer(1), ncat, max(cats))

    tname <- c("<leaf>", dimnames(X)[[2]])
    
    #if (cpcol==3) temp <- c("CP", "nsplit", "rel error")
    #else          temp <- c("CP", "nsplit", "rel error", "xerror", "xstd")
	temp <- c("CP", "nsplit", "rel error", "xerror", "xstd")

    dimnames(rp$cptable) <- list(temp, 1L:numcp)

    # R change for empty-vector calculations.
    dn1 <- if(nsplit == 0L) character(0L) else tname[rp$isplit[,1L]+1L]
    splits<- matrix(c(rp$isplit[,2L:3L], rp$dsplit), ncol=5L,
                    dimnames = list(dn1,
                    c("count", "ncat", "improve", "index", "adj")))
    if( !is.null ( nrow( rp$inode ) ) ) {
    index <- rp$inode[,2]  #points to the first split for each node 
	} else { index <- 1 }

    # Now, make ordered categories look like categories again (a printout
    #  choice)
    nadd <- sum(isord[rp$isplit[,1L]])
    if (nadd >0) {
	newc <- matrix(integer(1), nadd, max(cats))
	cvar <- rp$isplit[,1L]
	indx <- isord[cvar]		     # vector of T/F
	cdir <- splits[indx,2L]               # which direction splits went
	ccut <- floor(splits[indx,4L])        # cut point
	splits[indx,2L] <- cats[cvar[indx]]   #Now, # of categories instead
	splits[indx,4L] <- ncat + 1L:nadd      # rows to contain the splits

	# Next 4 lines can be done without a loop, but become indecipherable
	for (i in 1L:nadd) {
	    newc[i, 1L:(cats[(cvar[indx])[i]])] <- -1*as.integer(cdir[i])
	    newc[i, 1L:ccut[i]] <- as.integer(cdir[i])
	    }
	if (ncat==0) catmat <- newc
	else         catmat <- rbind(rp$csplit, newc)
	ncat <- ncat + nadd
	}
    else catmat <- rp$csplit
    
    if (nsplit==0) {  #tree with no splits
	frame <- data.frame(row.names=1,
			    var=  "<leaf>",
			    n =   rp$inode[5L],
			    wt=   rp$dnode[3L],
			    dev=  rp$dnode[1L],
			    yval= rp$dnode[4L],
			    complexity=rp$dnode[2L],
			    ncompete  = pmax(0L, rp$inode[3L] - 1L),
			    nsurrogate=rp$inode[4L])
	}
    else {
	temp <- ifelse(index==0, 1, index)

	svar <- ifelse(index==0, 0, rp$isplit[temp,1L]) #var number

	frame <- data.frame(row.names=rp$inode[,1L],
			    var=  factor(svar, 0:ncol(X), tname),
			    n =   rp$inode[,5L],
			    wt=   rp$dnode[,3L],
			    dev=  rp$dnode[,1L],
			    yval= rp$dnode[,4L],
			    complexity=rp$dnode[,2L],
			    ncompete  = pmax(0L, rp$inode[,3L] - 1L),
			    nsurrogate=rp$inode[,4L])
	}
	#print ( tname )
	#print ( frame )
	
    if (method.int == 3L) {
        numclass <- init$numresp -1L
        # Create the class probability vector from the class counts, and
        #   add it to the results
        # The "pmax" one line down is for the case of a factor y which has
        #   no one at all in one of its classes.  Both the prior and the
        #   count will be zero, which led to a 0/0.
        temp <- rp$dnode[,-(1L:4L), drop = FALSE] %*% diag(init$parms$prior*
                                           sum(init$counts)/pmax(1,init$counts))
        yprob <- temp /rowSums(temp)   #necessary with altered priors
        yval2 <- matrix(rp$dnode[, -(1L:3L)], ncol=numclass+1)
	frame$yval2 <- cbind(yval2, yprob)
	}
    else if (init$numresp >1L) frame$yval2 <- rp$dnode[,-(1L:3L), drop = FALSE]
 
    if (is.null(init$summary))
	    stop("Initialization routine is missing the 'summary' function")
    if (is.null(init$print))
	    functions <- list(summary=init$summary)
    else    functions <- list(summary=init$summary, print=init$print)
    if (!is.null(init$text)) functions <- c(functions, list(text=init$text))
    if (method=="user")	functions <- c(functions, mlist)

    where <- rep ( 0 , length ( rp$which ) )
    for ( whe in 1 : length ( rp$which ) ) {
	   where [ whe ] <- which ( as.numeric( rownames( frame ) )== as.numeric( rp$which )[ whe ] )
    }

    names(where) <- row.names(m)

    if (nsplit == 0L) {  # no 'splits' component
	ans <- list(frame = frame,
		    where = where,
		    call=Call, terms=Terms,
		    cptable =  t(rp$cptable),
		    method = method,
		    parms  = init$parms,
		    control= controls,
		    functions= functions)
    } else {
	ans <- list(frame = frame,
		    where = where,
		    call=Call, terms=Terms,
		    cptable =  t(rp$cptable),
		    splits = splits,
		    method = method,
		    parms  = init$parms,
		    control= controls,
		    functions= functions)
	}

    if (ncat > 0) ans$csplit <- catmat + 2L
    if (model) {
	ans$model <- m
	if (missing(y)) y <- FALSE
    }

    if (y) ans$y <- Y
    if (x) {
	ans$x <- X
	ans$wt<- wt
    }
    ans$ordered <- isord

    if(!is.null(attr(m, "na.action")))
        ans$na.action <- attr(m, "na.action")


    if (!is.null(xlevels)) attr(ans, "xlevels") <- xlevels


    if(method=="class") attr(ans, "ylevels") <- init$ylevels

    ans0 <- ans
    class(ans) <- "rpart"

#return( list( ans0=ans0, ans =ans  ) ) 
ans


}


