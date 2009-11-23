##################################################################
#########   MainFunctionWrap.R
#########   functions:
#########         getGeneExp
#########         DEGseq
#################################################################
	
getGeneExp <- function(mapResultBatch, fileFormat="bed", readLength=32, strandInfo=FALSE,
                       refFlat, output=paste(mapResultBatch[1],".exp",sep=""), min.overlapPercent=1){
  cat("Please wait...\n")
  needStrand <- 0 # need not
  if(strandInfo == FALSE){
    needStrand <- 0  # need not
  }else{
    needStrand <- 1  # need
  }
  if(output == "none"){
     output <- paste(mapResultBatch, "exp", sep="")
  }
  cat("\n")
  kk <- .C(".getGeneExp",as.character(refFlat),as.character(mapResultBatch),as.integer(length(mapResultBatch)),
           as.character(output),as.character(fileFormat),as.integer(readLength),as.integer(needStrand),
           as.numeric(min.overlapPercent))
}

getExonAnnotationFile <- function(refFlat, output){
  kk <- .C(".getExonAnnotationFile",as.character(refFlat),as.character(output))
}

readGeneExp <- function(file, geneCol=1, valCol=2, label = NULL, header=TRUE, sep=""){
  rt1 <- read.table(file, header=header, sep=sep, row.names=NULL)
  exp_values <- as(rt1[valCol], "matrix")
  exp_values[is.na(exp_values)] <- 0
  Gene_names <- as(rt1[geneCol], "matrix")
  if(!is.null(label)){
     dimnames(exp_values) <- list(as.character(Gene_names), label)
  }else{
     dimnames(exp_values) <- list(as.character(Gene_names), dimnames(exp_values)[[2]])
  }
  if(!is.numeric(exp_values)){
     cat("Some invalid values in the matrix!\n")
     cat("Please check the input file!\n")
  }
  exp_values
}

Getlines <- function(file){
  rt <- read.table(file,header = FALSE)
  nrow(rt)
}

DEGseq <- function(mapResultBatch1, mapResultBatch2, fileFormat="bed", readLength=32,
                   strandInfo=FALSE, refFlat, groupLabel1="group1", groupLabel2="group2",
                   method=c("LRT", "CTR", "FET", "MARS", "MATR", "FC"), 
                   pValue=1e-3, zScore=4, qValue=1e-3, foldChange=4, thresholdKind=1,
                   outputDir="none", normalMethod=c("none", "loess", "median"),
                   depthKind=1, replicate1="none", replicate2="none",
                   replicateLabel1="replicate1", replicateLabel2="replicate2"){
  
  method <- match.arg(method)
  normalMethod <- match.arg(normalMethod)
  cat("Please wait...\n")
  if(outputDir == "none"){
     cat("\noutputDir is not definite!")
     stop("\n")
     cat("Only generate statistic summary report graphs!\n")
  }
  if(outputDir != "none"){
     dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
     dir.create(paste(outputDir,"/output",sep=""), showWarnings = FALSE, recursive = TRUE)
  }else{
     ## outputDir <- tempdir()
  }
  if((outputDir != "none")&&(file.access(outputDir, mode = 0) != 0)){
    cat("Can not creat ",outputDir, "\n")
    stop("\n")
  }
      
  cat("\n")
  cat("mapResultBatch1: \n")
  for(i in (1:length(mapResultBatch1))){
     cat("\t", mapResultBatch1[i], "\n")
  }
  cat("mapResultBatch2: \n")
  for(i in (1:length(mapResultBatch2))){
     cat("\t", mapResultBatch2[i], "\n")
  }
  cat("file format:",fileFormat,"\n")
  cat("refFlat: ",refFlat,"\n")
  if(method == "MATR"){
     if((replicate1 == "none")||(replicate2 == "none")){
        cat("You must provide two replicates for method MATR!")
        cat("\n")
        stop();
     }
     cat("replicate1: ",replicate1,"\n")
     cat("replicate2: ",replicate2,"\n")
  }

  
  if(strandInfo == TRUE){
     cat("Consider the strand information when count the reads mapped to genes!\n")
  }else{
     cat("Ignore the strand information when count the reads mapped to genes!\n")
  }

  flush.console();
  
  GeneExp1 <- paste(outputDir,"/",groupLabel1,".exp",sep="");
  GeneExp2 <- paste(outputDir,"/",groupLabel2,".exp",sep="");
  unlink(GeneExp1)
  unlink(GeneExp2)
  cat("Count the number of reads mapped to each gene ... \n");
  cat("This will take several minutes, please wait patiently!\n")
  flush.console();
  getGeneExp(mapResultBatch1, fileFormat, readLength, strandInfo, refFlat, GeneExp1)
  getGeneExp(mapResultBatch2, fileFormat, readLength, strandInfo, refFlat, GeneExp2)
  ControlExp1 <- paste(outputDir,"/", replicateLabel1, ".exp",sep="");
  ControlExp2 <- paste(outputDir,"/", replicateLabel2, ".exp",sep="");

  
  if(method == "MATR"){
     unlink(ControlExp1)
     unlink(ControlExp2)
     getGeneExp(replicate1, fileFormat, readLength, strandInfo, refFlat, ControlExp1)
     getGeneExp(replicate2, fileFormat, readLength, strandInfo, refFlat, ControlExp2)
  }else{
     ControlExp1 <- "none"
     ControlExp2 <- "none"
  }

  
  expCol1 <- seq(2, by=3, length=length(mapResultBatch1))
  expCol2 <- seq(2, by=3, length=length(mapResultBatch2))
  if(depthKind == 0){
     depth1 <- rep(0, length(mapResultBatch1))
     depth2 <- rep(0, length(mapResultBatch2))
     cDepth1 <- rep(0, length(replicate1))
     cDepth2 <- rep(0, length(replicate2))
  }else{
     depth1 <- rep(-1, length(mapResultBatch1))
     depth2 <- rep(-1, length(mapResultBatch2))
     cDepth1 <- rep(-1, length(replicate1))
     cDepth2 <- rep(-1, length(replicate2))
  }

  DEGexp(GeneExp1, 1, expCol1, depth1, groupLabel1,
         GeneExp2, 1, expCol2, depth2, groupLabel2,
	       method=method, pValue=pValue, zScore=zScore, foldChange=foldChange, qValue=qValue, thresholdKind=thresholdKind,
         replicate1=ControlExp1, geneColR1=1, expColR1=2, depthR1=cDepth1,
         replicate2=ControlExp2, geneColR2=1, expColR2=2, depthR2=cDepth2,
	       replicateLabel1=replicateLabel1, replicateLabel2=replicateLabel2,
         outputDir=outputDir, normalMethod=normalMethod)
}                  
