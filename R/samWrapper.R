keep <- function(groups, measure1, measure2){
   keepR <- rep(TRUE, nrow(groups))
   for(i in (1:nrow(groups))){
      geneRow <- as.vector(groups[i,])
      row1 <- geneRow[1:length(measure1)]
      cat(i,":",row1,"\t")
      row2 <- geneRow[(length(measure1)+1):(length(measure1)+length(measure2))]
      cat(row2,"\t")
      zero_count1 <- length(row1[row1==0])
      zero_count2 <- length(row2[row2==0])
      cat(zero_count1,":",zero_count2,"\n")
      # if((zero_count1 >= length(measure1)/2)&&(zero_count2 >= length(measure2)/2)){
       #   keepR[i] <- FALSE
      # }
      if((zero_count1 >= 1)||(zero_count2 >= 1)){
          keepR[i] <- FALSE
      }
   }
   keepR 
}


####################################################

samWrapper <- function(geneExpFile1, geneCol1=1, expCol1=2, measure1=rep(1, length(expCol1)),
                    geneExpFile2, geneCol2=1, expCol2=2, measure2=rep(2, length(expCol2)),
                    header=TRUE, sep="", paired=FALSE, s0=NULL, s0.perc=NULL, nperms=100, 
                    testStatistic=c("standard","wilcoxon"), max.qValue=1e-3, min.foldchange=0,
                    logged2=FALSE, output){
   library(samr)
   cat("geneExpFile1:", geneExpFile1, "\n")
   cat("geneCol1:", geneCol1, "\n")
   cat("expCol11:", expCol1, "\n")
   cat("measure1:", measure1, "\n")
   
   cat("geneExpFile2:", geneExpFile2, "\n")
   cat("geneCol2:", geneCol2, "\n")
   cat("expCol2:", expCol2, "\n")
   cat("measure2:", measure2, "\n")
   
   cat("paired:", paired, "\n")
   cat("nperms:", nperms, "\n")
   
   testStatistic <- match.arg(testStatistic)
   if(paired==FALSE){
      cat("testStatistic", testStatistic, "\n")
   }
    
   cat("max.qValue:", max.qValue, "\n")
   
   if((min.foldchange < 1)&&(min.foldchange > 0)){
      min.foldchange <- 1/min.foldchange
   }
   cat("min.foldchange:", min.foldchange, "\n")
   cat("output:", output, "\n")
   cat("Please wait...\n")
   
   group1 <- readGeneExp(geneExpFile1, geneCol=geneCol1, valCol=expCol1, header=TRUE, sep=sep)
   d1 <- as.character(dimnames(group1)[[1]])
   d2 <- as.character(dimnames(group1)[[2]])
   d2 <- d2[2:ncol(group1)]
   group1 <- as(group1[,2:ncol(group1)], "matrix")
   group1 <- matrix(as.numeric(group1), ncol=ncol(group1))
   dimnames(group1) <- list(d1, d2)
   if(sum(is.na(group1)) > 0){
     cat("Some invalid values in the matrix!\n")
     cat("Please check the input file!\n")
     stop("");
   }             
   group2 <- readGeneExp(geneExpFile2, geneCol=geneCol2, valCol=expCol2, header=TRUE, sep=sep)
   d1 <- as.character(dimnames(group2)[[1]])
   d2 <- as.character(dimnames(group2)[[2]])
   d2 <- d2[2:ncol(group2)]
   group2 <- as(group2[,2:ncol(group2)], "matrix")
   group2 <- matrix(as.numeric(group2), ncol=ncol(group2))
   dimnames(group2) <- list(d1, d2)
   if(sum(is.na(group2)) > 0){
     cat("Some invalid values in the matrix!\n")
     cat("Please check the input file!\n")
     stop("");
   } 
   ########################################################
   ########################################################
   min_value <- min(min(group1[group1!=0]), min(group2[group2!=0]))
   groups <- cbind(group1, group2)
   len <- length(groups[groups==0])
   # groups <- groups[keep(groups,measure1,measure2),]
   if(len > 0){
     groups[groups==0] <- sample(seq(0,min_value,by=min_value/(len*100)),length(groups[groups==0]))
   }
   Measure <- as.vector(cbind(matrix(measure1, nrow=1), matrix(measure2, nrow=1)))
   if(paired==FALSE){
      check.format(Measure, "Two class unpaired")
   }else{
      check.format(Measure, "Two class paired")
   }
   geneNames <- as.vector(dimnames(groups)[[1]])   
   data <- list(x=groups, y=Measure, geneid=geneNames, genenames=geneNames, logged2=logged2)
   if(paired==FALSE){        
      samr.obj<-samr(data, resp.type="Two class unpaired", nperms=nperms, testStatistic=testStatistic, s0=s0, s0.perc=s0.perc)
   }else{
      samr.obj<-samr(data, resp.type="Two class paired", nperms=nperms, testStatistic=testStatistic, s0=s0, s0.perc=s0.perc)
   }
   
   delta.table <- samr.compute.delta.table(samr.obj, nvals=1000)
   siggenes.table <- samr.compute.siggenes.table(samr.obj, 0.4, data, delta.table, 
                                                 min.foldchange, all.genes=TRUE) 
   ##################################################################################
   ####### 	plot using the method samr.plot2 in DEGseq
   # p=length(samr.obj$tt)
   # pup=(1:p)[samr.obj$tt>=0]
   # plo=(1:p)[samr.obj$tt<0]
   # sig=list(pup=pup, plo=plo)
   
   # qvalues=qvalue.func(samr.obj,sig, delta.table)  ## can not be used. samr does not export the method.
   # samrWrapper.plot(samr.obj, sig, qvalues, max.qValue, min.foldchange)
   ##################################################################################
   all <- rbind(siggenes.table$genes.up, siggenes.table$genes.lo)
   # sort by gene names
   all[rank(all[,2]),] <- all  
   # sort by qvalue
   all[, ncol(all)] <- as.numeric(all[, ncol(all)])/100
   all[rank(as.numeric(all[, ncol(all)]), ties.method="first"),] <- all 
   all <- all[, 3:ncol(all)]
   sig1 <- as.numeric(all[, ncol(all)]) <= max.qValue
   if(min.foldchange>0){
       sig2 <- as.numeric(all[,ncol(all)-1]) >= min.foldchange
       sig3 <- as.numeric(all[,ncol(all)-1]) <= 1/min.foldchange
       sig1 <- sig1&(sig2|sig3)
   }
   all <- cbind(all, sig1)
   cat(file=output,"\"Gene Name\"\t\"Score(d)\"\t\"Numerator(r)\"\t\"Denominator(s+s0)\"\t\"Fold Change\"\t\"q-value\"\t\"Signature")
   cat(file=output,"(foldchange >=", min.foldchange, " && qValue <=", max.qValue,")\"\n", append=TRUE, sep="")
   write.table(all, file=output, quote=FALSE, sep="\t", row.names=FALSE,col.names=FALSE, append=TRUE)     
}

########################################
####################################################
samrWrapper.plot <- function(samr.obj, sig, qvalues, max.qValue, min.foldchange){

## make observed-expected plot
## takes foldchange into account too
  
  pup <- sig$pup
  plo <- sig$plo
  
  pupQ <- pup[qvalues$qvalue.up <= max.qValue]
  ploQ <- plo[qvalues$qvalue.lo <= max.qValue]
  
  pupT <- pup[(samr.obj$foldchange[pup] < min.foldchange)&(qvalues$qvalue.up <= max.qValue)]
  ploT <- plo[(samr.obj$foldchange[plo] > 1/min.foldchange)&(qvalues$qvalue.lo <= max.qValue)]
  
  #cat("pup:", pup, "\n")
  #cat("plo:", plo, "\n")
  
  col=rep(1,length(samr.obj$evo))
  pch=rep(21, length(samr.obj$evo))
  col[ploQ]=3
  pch[ploQ]=22
  col[pupQ]=2
  pch[pupQ]=22
  #col[ploT]=4
  #col[pupT]=6
  col.ordered=col[order(samr.obj$tt)]
  pch.ordered=pch[order(samr.obj$tt)]
  ylims <- range(samr.obj$tt)
  xlims <- range(samr.obj$evo)
  plot(samr.obj$evo,sort(samr.obj$tt),xlab="expected score", ylab="observed score",ylim=ylims, 
         xlim=xlims, type="n")
  points(samr.obj$evo,sort(samr.obj$tt),col=col.ordered, pch=pch.ordered)
  abline(0,1)
}
