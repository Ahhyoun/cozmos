

rpart.class <- function(y, offset, parms, wt)
{
  if (!is.null(offset)) stop("No offset allowed in classification models")
  
  fy <- as.factor(y)
  y <- as.integer(fy)
  numclass <- max(y[!is.na(y)])
  counts <- tapply(wt, factor(y, levels=1:numclass), sum)
  counts <- ifelse(is.na(counts), 0, counts) #in case of an empty class
  if (missing(parms) || is.null(parms))
    parms <- list(prior=counts/sum(counts),
                  loss=matrix(rep(1,numclass^2)-diag(numclass),numclass),
                  split=1)
  else if (is.list(parms)) {
    if (is.null(names(parms))) stop("The parms list must have names")
    temp <- pmatch(names(parms), c("prior", "loss", "split"), nomatch=0L)
    if (any(temp==0L))
      stop(gettextf("'parms' component not matched: %s",
                    names(parms)[temp == 0L]), domain = NA)
    names(parms) <- c("prior", "loss", "split")[temp]
    
    if (is.null(parms$prior)) temp <- c(counts/sum(counts))
    else {
      temp <- parms$prior
      if (sum(temp) !=1) stop("Priors must sum to 1")
      if (any(temp<0)) stop("Priors must be >= 0")
      if (length(temp) != numclass) stop("Wrong length for priors")
    }
    
    if (is.null(parms$loss)) temp2<- 1 - diag(numclass)
    else {
      temp2 <- parms$loss
      if (length(temp2) != numclass^2)
        stop("Wrong length for loss matrix")
      temp2 <- matrix(temp2, ncol=numclass)
      if (any(diag(temp2)!=0))
        stop("Loss matrix must have zero on diagonals")
      if (any(temp2<0))
        stop("Loss matrix cannot have negative elements")
      if (any(rowSums(temp2)==0))
        stop("Loss matrix has a row of zeros")
    }
    
    if (is.null(parms$split)) temp3 <- 1L
    else {
      temp3 <- pmatch(parms$split, c("gini", "information"))
      if (is.null(temp3)) stop("Invalid splitting rule")
    }
    parms <- list(prior=temp, loss=matrix(temp2,numclass), split=temp3)
  }
  else stop("Parameter argument must be a list")
  
  list(y=y, parms=parms, numresp=numclass+1L, counts=counts,
       ylevels= levels(fy), numy=1L,
       print = function(yval, ylevel, digits) {
         if (is.null(ylevel))
           temp <- as.character(yval[,1L])
         else    temp <- ylevel[yval[,1L]]
         
         nclass <- (ncol(yval) -1L)/2
         if (nclass <5) {
           yprob <- format(yval[, 1L+nclass + 1L:nclass],
                           digits=digits,nsmall=digits)
         }
         else yprob <- formatg(yval[, 1L+nclass + 1L:nclass], digits=2)
         if(is.null(dim(yprob)))    # yprob is a vector
           yprob <- matrix(yprob, ncol=length(yprob))
         temp <- paste(temp, ' (', yprob[,1L], sep='')
         for(i in 2L:ncol(yprob))
           temp  <- paste(temp, yprob[, i], sep=' ')
         temp <- paste(temp, ")", sep="")
         temp
       },
       summary= function(yval, dev, wt, ylevel, digits) {
         nclass <- (ncol(yval)-1L) /2
         group <- yval[, 1L]
         counts <- yval[, 1L+ (1L:nclass)]
         yprob  <- yval[, 1L+nclass + 1L:nclass]
         if(!is.null(ylevel)) group <- ylevel[group]
         
         temp1 <- formatg(counts, format="%5g")
         temp2 <- formatg(yprob,  format="%5.3f")
         if (nclass >1) {
           temp1 <- apply(matrix(temp1, ncol=nclass), 1L,
                          paste, collapse=' ')
           temp2 <- apply(matrix(temp2, ncol=nclass), 1L,
                          paste, collapse=' ')
         }
         paste("  predicted class=", format(group, justify='left'),
               "  expected loss=", formatg(dev/wt, digits),"\n",
               "    class counts: ", temp1,"\n",
               "   probabilities: ", temp2,
               sep='')
       },
       text= function(yval, dev, wt, ylevel, digits, n, use.n) {
         
         nclass <- (ncol(yval)-1L) /2
         group <- yval[, 1L]
         counts <- yval[, 1L+ (1L:nclass)]
         if(!is.null(ylevel)) group <- ylevel[group]
         
         temp1 <- formatg(counts, digits)
         if (nclass >1) {
           temp1 <- apply(matrix(temp1, ncol=nclass), 1L,
                          paste, collapse='/')
         }
         
         if(use.n)  {out <- paste(format(group, justify='left'),"\n",
                                  temp1,sep="")}      else
                                  {out <- format(group,justify="left")}
         
         return(out)
       })
}


rpart.control <-
  function(minsplit=20, minbucket= round(minsplit/3), cp=.01,
           maxcompete=4, maxsurrogate=5, usesurrogate=2, xval=10,
           surrogatestyle =0, maxdepth=30, ... )
  {
    
    if (maxcompete<0) {
      warning("The value of 'maxcompete' supplied is < 0; the value 0 was used instead")
      maxcompete <-0
    }
    if (any(xval<0)) {
      warning("The value of 'xval' supplied is < 0; the value 0 was used instead")
      xval <-0
    }
    if (maxdepth > 30) stop("Maximum depth is 30")
    if (maxdepth < 1)  stop("Maximum depth must be at least 1")
    
    if (missing(minsplit) && !missing(minbucket)) minsplit <- minbucket*3
    
    if((usesurrogate < 0) || (usesurrogate > 2)) {
      warning("The value of 'usesurrogate' supplied was out of range, the default value of 2 is used instead.")
      usesurrogate <- 2
    }
    if((surrogatestyle < 0) || (surrogatestyle > 1)) {
      warning("The value of 'surrogatestyle' supplied was out of range, the default value of 0 is used instead.")
      surrogatestyle <- 0
    }
    
    # Because xval can be of length either 1 or n, and the C code
    #   refers to parameters by number, i.e., "opt[5]" in rpart.c,
    #   the xval parameter should always be last on the list.
    list(minsplit=minsplit, minbucket=minbucket, cp=cp,
         maxcompete=maxcompete, maxsurrogate=maxsurrogate,
         usesurrogate=usesurrogate,
         surrogatestyle=surrogatestyle, maxdepth=maxdepth, xval=xval )
  }

rpart.matrix <- function(frame)
{
  if(!inherits(frame, "data.frame")) return(as.matrix(frame))
  frame$"(weights)" <- NULL
  terms <- attr(frame, "terms")
  if(is.null(terms)) predictors <- names(frame)
  else {
    a <- attributes(terms)
    predictors <- as.character(a$variables)[-1L] # R change
    ## and this might include backquotes
    predictors <- sub("^`(.*)`$", "\\1", predictors)
    removals <- NULL
    if((TT <- a$response) > 0L) {
      removals <- TT
      frame[[predictors[TT]]] <- NULL
    }
    if(!is.null(TT <- a$offset)) {
      removals <- c(removals, TT)
      frame[[predictors[TT]]] <- NULL
    }
    if(!is.null(removals)) predictors <- predictors[ - removals]
    labels <- a$term.labels
    if(abs(length(labels)-length(predictors))>0L)
      predictors <- predictors[match(labels,predictors)]
  }
  
  factors <- sapply(frame, function(x) !is.null(levels(x)))
  characters <- sapply(frame, is.character)
  if(any(factors | characters)) {
    ## change characters to factors
    for (preds in predictors[characters])
      frame[[preds]] <- as.factor(frame[[preds]])
    factors <- factors | characters
    column.levels <- lapply(frame[factors], levels)
    
    ## Now make them numeric
    for (preds in predictors[factors])
      frame[[preds]] <- as.numeric(frame[[preds]])
    x <- as.matrix(frame)
    attr(x, "column.levels") <- column.levels
  }
  else x <- as.matrix(frame[predictors])
  class(x) <- "rpart.matrix"
  x
}


