######################################
#### functions.R
####
#### other functions in DEGseq
####################################### 
## funtions for laneRaw

#######################################

ReadLane <- function(file1, geneCol=1, valCol=2, label = NULL, header=TRUE, sep=""){
   rt1 <- read.table(file1,header=header, sep=sep, row.names=NULL)
   
   exp_values <- as(rt1[valCol], "matrix")
   exp_values <- matrix(as.numeric(exp_values), ncol=ncol(exp_values))
   if(sum(is.na(exp_values)) > 0){
     cat("Some invalid values in the matrix!\n")
     cat("Please check the input file!\n")
     stop("");
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
#######################################

ReadLane2 <- function(geneExpMatrix, geneCol=1, valCol=2, label = NULL){
   
   exp_values <- as(geneExpMatrix[, valCol], "matrix")
   exp_values <- matrix(as.numeric(exp_values), ncol=ncol(exp_values))
  
   if(sum(is.na(exp_values)) > 0){
     cat("Some invalid values in the matrix!\n")
     cat("Please check the input file!\n")
     stop("");
   }
   exp_values[is.na(exp_values)] <- 0
   exp_values_sum <- exp_values[,1]

   if(length(valCol) >= 2){
      for(i in (2:length(valCol))){
          exp_values_sum <- exp_values_sum + exp_values[,i]
      }
   }
   Gene_names <- as(geneExpMatrix[, geneCol], "matrix")
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
#################################################################
exportAlnAsBed <- function(aln, output){
  if(!is(aln, "AlignedRead")){
     stop("aln is not an object of class AlignedRead\n")
  }
  aln <- aln[!is.na(position(aln))]
  string <- as.character(sread(aln))
  chromosome <- as.character(chromosome(aln))
  start <- position(aln)-1
  end <- start + nchar(string)
  start <- formatC(start, format="f", digits=11, drop0trailing=FALSE)
  end <- formatC(end, format="f", digits=11, drop0trailing=FALSE)
  id <- as.character(id(aln))
  score <- as.numeric(quality(alignQuality(aln)))
  score[is.na(score)] <- 0
  strand <- as.character(strand(aln))
  strand[(strand != "+")&(strand != "-")] <- "+"
  result <- cbind(chromosome, start, end, id, score, strand)
  write.table(result, col.name=FALSE, row.name=FALSE, file=output,quote=FALSE, sep="\t")
}
#######################################
DEV_OFF <- function(off=TRUE, dev_cur=2){    
 for(i in dev.list())
  if((i != 1) && (off) && (i > dev_cur)){
    dev.off(i)
  }
}
open_html <- function(path){
  if(.Platform$OS.type == "windows"){
     if(substr(path,nchar(path),nchar(path)) == "/"){
        path <- substr(path, 1, nchar(path)-1)
     }
     if(file.access("C:/Program Files/Internet Explorer/iexplore.exe", mode = 0)==0){
        ie_path <- "\"C:/Program Files/Internet Explorer/iexplore.exe\" "
        quotation_str <- "\"";
        out_html_path <- paste(quotation_str,path,"/output.html",quotation_str,sep="");
        system(paste(ie_path,out_html_path,sep=""),invisible=FALSE,wait=FALSE,show.output.on.console = FALSE)
     }
  }
}

######################################################