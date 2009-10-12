########################################################################
#
#
# The functions are modified from the package marray (Yang Y.H et al.).
#
#
# cbind method
########################################################################
cbind2laneRaw <- function(newx, x){
   if(length(expVals(x))!=0)
      expVals(newx) <- cbind(expVals(newx), expVals(x))
   if(length(laGenes(x))!=0)
      laGenes(newx) <- cbind(laGenes(newx), laGenes(x))
   if(length(laNotes(x)) != 0)
      slot(newx,"laNotes")<- as.character(rbind(as.matrix(laNotes(newx)), as.matrix(laNotes(x))))
   dimnames(slot(newx,"expVals")) <- list(as.character(slot(newx,"laGenes")[,1]),slot(newx,"laNotes"))
   return(newx)
}

cbind.laneRaw <- function(..., deparse.level = 1){
    data <- list(...)
    newx <- data[[1]]
    for(x in data[2:length(data)]){
        if(length(expVals(x))!=0)
          expVals(newx) <- cbind(expVals(newx), expVals(x))
        if(length(laGenes(x))!=0)
          laGenes(newx) <- cbind(laGenes(newx), laGenes(x))
        if(length(laNotes(x)) != 0)
          slot(newx,"laNotes")<- as.character(rbind(as.matrix(laNotes(newx)), as.matrix(laNotes(x))))
    }
    dimnames(slot(newx,"expVals")) <- list(as.character(slot(newx,"laGenes")[,1]),slot(newx,"laNotes"))
    return(newx)
}

cbind.lanePair <- function(..., deparse.level = 1){
    data <- list(...)
    newx <- data[[1]]
    for(x in data[2:length(data)]){
        if(length(expVals1(x))!=0)
          expVals1(newx) <- cbind(expVals1(newx), expVals1(x))
        if(length(expVals2(x))!=0)
          expVals2(newx) <- cbind(expVals2(newx), expVals2(x))
        if(length(PairGenes(x))!=0)
          PairGenes(newx) <- cbind(PairGenes(newx), PairGenes(x))
        if(length(PairNotes(x)) != 0)
          PairNotes(newx) <- as.character(rbind(as.matrix(PairNotes(newx)), as.matrix(PairNotes(x))))
    }
    dimnames(slot(newx,"expVals1")) <- list(as.character(slot(newx,"PairGenes")[,1]),slot(newx,"PairNotes"))
    dimnames(slot(newx,"expVals2")) <- list(as.character(slot(newx,"PairGenes")[,1]),slot(newx,"PairNotes"))
    return(newx)
}

cbind.PairNorm <-  function(..., deparse.level = 1){
    data <- list(...)
    newx <- data[[1]]
    for(x in data[2:length(data)]){
        if(length(AVal(x))!=0)
          AVal(newx) <- cbind(AVal(newx), AVal(x))
        if(length(MVal(x))!=0)
          MVal(newx) <- cbind(MVal(newx), MVal(x))
        if(length(pMloc(x))!=0)
          pMloc(newx) <- cbind(newx@pMloc, x@pMloc)
        if(length(pMscale(x))!=0)
          pMscale(newx) <- cbind(pMscale(newx), pMscale(x))
        if(length(PairGenes(x))!=0)
          PairGenes(newx) <- cbind(PairGenes(newx), PairGenes(x))
        if(length(PairNotes(x)) != 0)
          PairNotes(newx) <- as.character(rbind(as.matrix(PairNotes(newx)), as.matrix(PairNotes(x))))
    }
    newx@PairNormCall=NULL
    dimnames(slot(newx,"AVal")) <- list(as.character(slot(newx,"PairGenes")[,1]),slot(newx,"PairNotes"))
    dimnames(slot(newx,"MVal")) <- list(as.character(slot(newx,"PairGenes")[,1]),slot(newx,"PairNotes"))
    return(newx)
}


getPairs <- function(lanebatch1, lanebatch2, smallest_value=0, method="LRT"){
     newPairs <- new("lanePair")
     lane1_i = lanebatch1
     lane2_j = lanebatch2
     one_pair <- new("lanePair")
     expVals1(one_pair) <- expVals(lane1_i)
     expVals2(one_pair) <- expVals(lane2_j)
     if(smallest_value==0){
        smallest_value1 <- min(expVals1(one_pair)[expVals1(one_pair)!=0])
        smallest_value2 <- min(expVals2(one_pair)[expVals2(one_pair)!=0])
        smallest_value <- min(smallest_value1,smallest_value2)
     }
     #expVals1(one_pair)[expVals1(one_pair)==0] <- smallest_value
     #expVals2(one_pair)[expVals2(one_pair)==0] <- smallest_value
            
     #expVals1(one_pair)[expVals1(one_pair)==0] <- NA
     #expVals2(one_pair)[expVals2(one_pair)==0] <- NA
            
     #cat("count 0 in expVals1:", length(expVals1(one_pair)[expVals1(one_pair)==0]),"\n")
     #cat("count 0 in expVals2:", length(expVals2(one_pair)[expVals2(one_pair)==0]),"\n")
     #cat("count 0 in 1 and 2",length(expVals1(one_pair)[(expVals(lane1_i)==0)&(expVals(lane2_j)==0)]),"\n")
            
     expVals1(one_pair)[(expVals(lane1_i)==0)&(expVals(lane2_j)==0)] <- NA
     expVals2(one_pair)[(expVals(lane1_i)==0)&(expVals(lane2_j)==0)] <- NA

     if((method == "MARS")||(method == "MATR")){
       expVals1(one_pair)[(expVals(lane1_i)==0)&(expVals(lane2_j)>=5)] <- 0.5
       expVals2(one_pair)[(expVals(lane1_i)>=5)&(expVals(lane2_j)==0)] <- 0.5
     }
     	         
     PairGenes(one_pair) <- laGenes(lane1_i)
     PairNotes(one_pair) <- paste(laNotes(lane1_i),"VS",laNotes(lane2_j))

     newPairs <- one_pair
     dimnames(slot(newPairs,"expVals1")) <- list(as.character(slot(newPairs,"PairGenes")[,1]),slot(newPairs,"PairNotes"))
     dimnames(slot(newPairs,"expVals2")) <- list(as.character(slot(newPairs,"PairGenes")[,1]),slot(newPairs,"PairNotes"))
     newPairs
}
