######################################
#### functions.R
####
#### other functions in DEGseq
####################################### 
## funtions for laneRaw

#######################################

ReadLane <- function(file1, geneCol=1, valCol=2, label = NULL, header=TRUE, sep=""){
   rt1 <- read.table(file1,header=header, sep=sep)
   
   exp_values <- as(rt1[valCol], "matrix")
   if(!is.numeric(exp_values)){
     cat("Some invalid values in the matrix!\n")
     cat("Please check the input file!\n")
   }
   exp_values[is.na(exp_values)] <- 0
   exp_values_sum <- exp_values[,1]

   if(length(valCol) >= 2){
      for(i in (2:length(valCol))){
          exp_values_sum <- exp_values_sum + exp_values[,i]
      }
   }
   Gene_names <- as(rt1[geneCol], "matrix")
  
   laneRawNew<-new("laneRaw")
   expVals(laneRawNew) <- as(exp_values_sum, "matrix")
   laGenes(laneRawNew) <- Gene_names
   laNotes(laneRawNew) <- label
   dimnames(slot(laneRawNew,"expVals")) <- list(as.character(slot(laneRawNew,"laGenes")[,1]),slot(laneRawNew,"laNotes"))
   laneRawNew
}
###########################################################################

## funtions for class Pair2Norm
#############################                  
Pair2Norm <-function(from){
    lnorm<-new("PairNorm", AVal=AVal(from), MVal=MVal(from),PairGenes=PairGenes(from), PairNotes=PairNotes(from))
    dimnames(slot(lnorm,"AVal")) <- list(as.character(slot(lnorm,"PairGenes")[,1]),slot(lnorm,"PairNotes"))
    dimnames(slot(lnorm,"MVal")) <- list(as.character(slot(lnorm,"PairGenes")[,1]),slot(lnorm,"PairNotes"))
    lnorm
}

#######################
#######################################################
## funtions for class PairResult
##########################
Norm2Result <- function(from){
    lnorm<-new("PairResult", AVal=AVal(from), MVal=MVal(from),PairGenes=PairGenes(from), PairNotes=PairNotes(from),
                pMloc=pMloc(from),pMscale=pMscale(from),PairNormCall=slot(from,"PairNormCall"))
    lnorm
}

#######################
###########################################################################
