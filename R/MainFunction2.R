#####################################################
########### MainFunction.R
########### functions:
###########            DEGexp2: the old version of DEGexp
#####################################################
DEGexp2 <- function(geneExpFile1, geneCol1=1, expCol1=2, depth1=rep(0, length(expCol1)), groupLabel1="group1",
                   geneExpFile2, geneCol2=1, expCol2=2, depth2=rep(0, length(expCol2)), groupLabel2="group2",
                   header=TRUE, sep="", method=c("LRT", "CTR", "FET", "MARS", "MATR", "FC"),
                   pValue=1e-3, zScore=4, qValue=1e-3, foldChange=4,
                   thresholdKind=1, outputDir="none", normalMethod=c("none", "loess", "median"),
                   replicate1="none", geneColR1=1, expColR1=2, depthR1=rep(0, length(expColR1)), replicateLabel1="replicate1",
                   replicate2="none", geneColR2=1, expColR2=2, depthR2=rep(0, length(expColR2)), replicateLabel2="replicate2", rawCount=TRUE){

 dev_cur <- dev.cur();
 cat("Please wait...\n")
 method <- match.arg(method)
 normalMethod <- match.arg(normalMethod)
 for(i in (1:length(expCol1))){
   if(depth1[i] == -1){
      rt1 <- read.table(geneExpFile1,header=header,sep=sep,row.names=NULL)
      exp_values <- as(rt1[expCol1[i]+2], "matrix")
      exp_values[is.na(exp_values)] <- 0
      depth1[i] <- as.numeric(exp_values[1])
   }
   if(depth1[i] == 0){
      rt1 <- read.table(geneExpFile1,header=header,sep=sep,row.names=NULL)
      exp_values <- as(rt1[expCol1[i]], "matrix")
      exp_values[is.na(exp_values)] <- 0
      depth1[i] <- sum(exp_values)
      warning_msg <- 1
   }
 }

 for(i in (1:length(expCol2))){
   if(depth2[i] == -1){
      rt1 <- read.table(geneExpFile2,header=header,sep=sep,row.names=NULL)
      exp_values <- as(rt1[expCol2[i]+2], "matrix")
      exp_values[is.na(exp_values)] <- 0
      depth2[i] <- as.numeric(exp_values[1])
   }
   if(depth2[i] == 0){
      rt1 <- read.table(geneExpFile2,header=header,sep=sep,row.names=NULL)
      exp_values <- as(rt1[expCol2[i]], "matrix")
      exp_values[is.na(exp_values)] <- 0
      depth2[i] <- sum(exp_values)
   }
 }
 
 if(method == "MATR"){
   for(i in (1:length(expColR1))){
       if(depthR1[i] == -1){
           rt1 <- read.table(replicate1,header=header,sep=sep,row.names=NULL)
           exp_values <- as(rt1[expColR1[i]+2], "matrix")
	         exp_values[is.na(exp_values)] <- 0
           depthR1[i] <- as.numeric(exp_values[1])
       }
       if(depthR1[i] == 0){
           rt1 <- read.table(replicate1,header=header,sep=sep,row.names=NULL)
           exp_values <- as(rt1[expColR1[i]], "matrix")
	         exp_values[is.na(exp_values)] <- 0
           depthR1[i] <- sum(exp_values)
       }
   }
   for(i in (1:length(expColR2))){
       if(depthR2[i] == -1){
           rt1 <- read.table(replicate2,header=header,sep=sep,row.names=NULL)
           exp_values <- as(rt1[expColR2[i]+2], "matrix")
	         exp_values[is.na(exp_values)] <- 0
           depthR2[i] <- as.numeric(exp_values[1])
       }
       if(depthR2[i] == 0){
           rt1 <- read.table(replicate2,header=header,sep=sep,row.names=NULL)
           exp_values <- as(rt1[expColR2[i]], "matrix")
	         exp_values[is.na(exp_values)] <- 0
           depthR2[i] <- sum(exp_values)      
       }
   }
 }

 cat("\ngeneExpFile1: ",geneExpFile1,"\n")
 cat("gene id column in geneExpFile1: ",geneCol1,"\n")
 cat("expression value column(s) in geneExpFile1:",expCol1,"\n")
 cat("total number of reads uniquely mapped to genome obtained from sample1:",depth1,"\n")

 cat("\ngeneExpFile2: ",geneExpFile2,"\n")
 cat("gene id column in geneExpFile2: ",geneCol2,"\n")
 cat("expression value column(s) in geneExpFile2:",expCol2,"\n")
 cat("total number of reads uniquely mapped to genome obtained from sample2:",depth2,"\n\n")

 if(replicate1 != "none"){
    cat("replicate1: ",replicate1,"\n")
    cat("gene id column in the expression file for replicate1: ",geneColR1,"\n")
    cat("expression value column(s) in the expression file for replicate1:",expColR1,"\n")
    cat("total number of reads uniquely mapped to genome obtained from replicate1:",depthR1,"\n\n")
    cat("replicate2: ",replicate2,"\n")
    cat("gene id column in the expression file for replicate2: ",geneColR2,"\n")
    cat("expression value column(s) in the expression file for replicate2:",expColR2,"\n")
    cat("total number of reads uniquely mapped to genome obtained from replicate2:",depthR2,"\n\n")
 }

 cat("method to identify differentially expressed genes: ",method,"\n")
 
 pValue_threshold <- pValue
 qValue_threshold <- qValue
 threshold <- zScore
 
 if((thresholdKind != 1)&&(thresholdKind != 2)&&(thresholdKind != 3)&&(thresholdKind != 4)){
     cat("Wrong value for thresholdKind!\n")
     cat("Use default value for thresholdKind.\n")
     thresholdKind <- 1
 }
 if(method != "FC"){
   if(thresholdKind == 2){
      cat("zScore threshold:",threshold,"\n")
      pValue_threshold <- 2*pnorm(-abs(threshold))
   }
   if(thresholdKind == 1){
      cat("pValue threshold:",pValue_threshold,"\n")
      threshold <- abs(qnorm(pValue_threshold/2))
   }
   if(thresholdKind == 3){
      pValue_threshold <- 0
      threshold <- 4
      cat("qValue threshold (Benjamini et al. 1995):",qValue_threshold,"\n")
   }
   if(thresholdKind == 4){
      pValue_threshold <- 0
      threshold <- 4
      cat("qValue threshold (Storey et al. 2003):",qValue_threshold,"\n")
   }
 }
 
 if(method == "FC"){
    cat("log2 fold change:", foldChange,"\n")
 }

 cat("output directory:",outputDir,"\n")
 if(outputDir=="none"){
    cat("\nThe outputDir is not specified! ")
    # cat("TheOnly generate statistic summary report graphs!\n")
 }
 if((method != "FET")&&(method != "LRT")&&(method != "MATR")
   &&(method != "CTR")&&(method != "FC")&&(method != "MARS")){
     cat(method,"\tis not a valid method name!\n")
     stop("\n")
 }
 if((normalMethod != "none")&&(normalMethod != "loess")&&(normalMethod != "median")){
     cat(normalMethod,"\tis not a valid normalization method name!\n")
     stop("\n")
 }
 if(normalMethod != "none"){
     cat("normalMethod:",normalMethod,"\n")
 }
 if(outputDir!="none"){
    dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
    dir.create(paste(outputDir,"/output",sep=""), showWarnings = FALSE, recursive = TRUE)
    if(file.access(outputDir, mode = 0) != 0){
       cat("Can not creat ",outputDir, "\n")
       stop("\n")
    }
 }
 cat("\n")
 cat("Please wait ...\n")
 flush.console();
 library(methods)

 label1 <- groupLabel1
 label1 <- substr(label1,1,14)
 label2 <- groupLabel2
 label2 <- substr(label2,1,14)
 if(label1 == label2){
    label1 <- "sample1"
    label2 <- "sample2"
 }
 Sample1 <- ReadLane(geneExpFile1,geneCol=geneCol1,valCol=expCol1,label=label1,header=header,sep=sep)
 Sample2 <- ReadLane(geneExpFile2,geneCol=geneCol2,valCol=expCol2,label=label2,header=header,sep=sep)
 png_new("/output/Sample1_hist.png",outputDir)
 scatterMain <- paste(label1," VS ",label2,sep="")
 #label1 <- paste("log2(",label1,")",sep="")
 #label2 <- paste("log2(",label2,")",sep="")
 if(outputDir != "none"){
    par(cex.axis = 1.5,cex.main=2,cex.lab=1.5,font.axis=2,font.lab=2,lwd=3)
 }else{
    #par(def.par)
 }
 op <- par(lwd=1)
 hist(LogVal(Sample1),main=label1,xlab="log2(Number of reads mapped to a gene)",col=4,breaks=100,freq=FALSE,ylim=c(0,0.5))
 par(op)
 logVal1 <- LogVal(Sample1);
 logVal1 <- logVal1[!is.na(logVal1)]
 lines(density(logVal1),col="red")
 png_new("/output/Sample2_hist.png",outputDir)
 op <- par(lwd=1)
 hist(LogVal(Sample2),main=label2,xlab="log2(Number of reads mapped to a gene)",col=4,breaks=100,freq=FALSE,ylim=c(0,0.5),lwd=1)
 par(op)
 logVal2 <- LogVal(Sample2);
 logVal2 <- logVal2[!is.na(logVal2)]
 lines(density(logVal2),col="red")
 DEV_OFF((outputDir != "none"), dev_cur=dev_cur)
 png_new("/output/SampleS_box.png",outputDir)
 op <- par(lwd=1)
 boxplot(cbind2laneRaw(Sample1,Sample2))
 par(op)
 DEV_OFF((outputDir != "none"), dev_cur=dev_cur)
 ylab <- paste("log2(read counts for each gene) in ",label1,sep="")
 xlab <- paste("log2(read counts for each gene) in ",label2,sep="")
 png_new("/output/Sample1_Sample2_compare.png",outputDir)
 LogVal1 <- log(expVals(Sample1),2)
 LogVal2 <- log(expVals(Sample2),2)
 min_value <- min(min(LogVal1[is.finite(LogVal1)]), min(LogVal2[is.finite(LogVal2)]))
 max_value <- max(max(LogVal1[is.finite(LogVal1)]), max(LogVal2[is.finite(LogVal2)]))
 xy_lim <- c(min_value, max_value)
 plot(LogVal2, LogVal1, pch=".",xlab=xlab,ylab=ylab,main=scatterMain,xlim=xy_lim,ylim=xy_lim)
 abline(0,1,lty=5,col="blue")
 DEV_OFF((outputDir != "none"), dev_cur=dev_cur)

 Criterion1 <- Sample1
 Criterion2 <- Sample2
 if(replicate1 != "none"){
    c_label1 <- replicateLabel1
    c_label1 <- substr(c_label1,1,14)
    c_label2 <- replicateLabel2
    c_label2 <- substr(c_label2,1,14)
    if(c_label1 == c_label2){
       c_label1 <- "replicate1"
       c_label2 <- "replicate2"
    }
    Criterion1 <- ReadLane(replicate1,geneCol=geneColR1,valCol=expColR1,label=c_label1,header=header,sep=sep)
    Criterion2 <- ReadLane(replicate2,geneCol=geneColR2,valCol=expColR2,label=c_label2,header=header,sep=sep)
    scatterMain <- paste(c_label1," VS ",c_label2,sep="")
    ylab <- paste("log2(read counts for each gene) in ",c_label1,sep="")
    xlab <- paste("log2(read counts for each gene) in ",c_label2,sep="")
    png_new("/output/control1_control2_compare.png",outputDir)
    LogVal1 <- log(expVals(Criterion1),2)
    LogVal2 <- log(expVals(Criterion2),2)
    min_value <- min(min(LogVal1[is.finite(LogVal1)]), min(LogVal2[is.finite(LogVal2)]))
    max_value <- max(max(LogVal1[is.finite(LogVal1)]), max(LogVal2[is.finite(LogVal2)]))
    xy_lim <- c(min_value, max_value)
    plot(LogVal2, LogVal1, pch=".",xlab=xlab,ylab=ylab,main=scatterMain,xlim=xy_lim,ylim=xy_lim)
    abline(0,1,lty=5,col="blue")
    DEV_OFF((outputDir != "none"), dev_cur=dev_cur);
 }

 PairS1S2 <- getPairs(Sample1,Sample2,method=method)
 PairC1C2 <- getPairs(Criterion1,Criterion2)


 count1 <- sum(depth1)
 count2 <- sum(depth2)

 switch(normalMethod,
         "none" = {
                    #png_new("/output/pre_norm.png",outputDir);
                    #plot(PairS1S2);
                    #LineOnPlot(PairS1S2,c("median"))
                    #LineOnPlot(PairS1S2,c("median","loess"))
                    PairS1S2_norm<-Pair2Norm(PairS1S2)
                    #png_new("/output/after_norm.png",outputDir)
                    #plot(PairS1S2_norm)
                    #LineOnPlot(PairS1S2_norm,c("median"))
                    #LineOnPlot(PairS1S2_norm,c("median","loess"))
                  },

         "median" = {
                     png_new("/output/pre_norm.png",outputDir)
                     plot(PairS1S2);
                     LineOnPlot(PairS1S2,"median");
                     PairS1S2_norm <- maNorm(PairS1S2,"median")
                     png_new("/output/after_norm.png",outputDir)
                     plot(PairS1S2_norm)
                     LineOnPlot(PairS1S2_norm,"median")
                    },

         "loess" = {
                     png_new("/output/pre_norm.png",outputDir)
                     plot(PairS1S2);
                     LineOnPlot(PairS1S2,"loess");
                     PairS1S2_norm <- maNorm(PairS1S2,"loess")
                     #points(AVal(PairS1S2_norm),pMloc(PairS1S2_norm),col=2)
                     png_new("/output/after_norm.png",outputDir)
                     plot(PairS1S2_norm)
                     LineOnPlot(PairS1S2_norm,"loess")
                   }
      )
      
 DEV_OFF((outputDir != "none"), dev_cur=dev_cur);

 if(method == "MATR"){
    switch(normalMethod,
         "none" = {
                    #png_new("/output/pre_norm_c.png",outputDir);
                    #plot(PairC1C2);
                    PairC1C2_norm<-Pair2Norm(PairC1C2)
                    if(rawCount==TRUE){
                       MVal(PairC1C2_norm) <- MVal(PairC1C2_norm)-log((depthR1/depthR2),2)+log((count1/count2),2)
                    }else{
                    }
                    #png_new("/output/after_norm_c.png",outputDir)
                    #plot(PairC1C2_norm)
                  },

         "median" = {
                     png_new("/output/pre_norm_c.png",outputDir)
                     plot(PairC1C2);
                     LineOnPlot(PairC1C2,"median");
                     PairC1C2_norm <- maNorm(PairC1C2,"median")
                     png_new("/output/after_norm_c.png",outputDir)
                     plot(PairC1C2_norm)
                     LineOnPlot(PairC1C2_norm,"median")
                    },

         "loess" = {
                     png_new("/output/pre_norm_c.png",outputDir)
                     plot(PairC1C2);
                     LineOnPlot(PairC1C2,"loess");
                     PairC1C2_norm <- maNorm(PairC1C2,"loess")
                     #points(AVal(PairS1S2_norm),pMloc(PairC1C2_norm),col=2)
                     png_new("/output/after_norm_c.png",outputDir)
                     plot(PairC1C2_norm)
                     LineOnPlot(PairC1C2_norm,"loess")
                   }
      )
      
   DEV_OFF((outputDir != "none"), dev_cur=dev_cur);
 }

 if(normalMethod != "none"){
    cat("Normalization done.\n")
 }
 cat("Identifying differentially expressed genes ...\nPlease wait patiently ...\n")
 flush.console();
 result_pic <- "none"
 result_score <- "none"
 if(outputDir != "none"){
    result_pic <- paste(outputDir,"/output/result.png",sep="")
    result_score <- paste(outputDir,"/output_score.txt",sep="")
 }
 switch(method,
        "FC" = {  
                  PairS1S2_result <- FindDiff_FC(PairS1S2_norm,threshold=foldChange,picfile=result_pic)
                  output_FC(PairS1S2_result,PairS1S2,result_score,threshold=foldChange)
                        },
        
        "MATR" = {
                  PairS1S2_result <- FindDiff_SD(PairS1S2_norm,PairC1C2_norm,,threshold=threshold,
                                                 pvalue_t=pValue_threshold,qvalue_t=qValue_threshold,
                                                 thresholdKind=thresholdKind,picfile=result_pic,control=1)
                  output_SD(PairS1S2_result,PairS1S2,result_score,pvalue_t=pValue_threshold,qvalue_t=qValue_threshold,thresholdKind=thresholdKind)
                        },
       
        "MARS" = {
                  PairS1S2_result <- FindDiff_MARS(PairS1S2_norm,count1=count1,count2=count2,threshold=threshold,
                                                   pvalue_t=pValue_threshold,qvalue_t=qValue_threshold,
                                                   thresholdKind=thresholdKind,picfile=result_pic)
                  output_SD(PairS1S2_result,PairS1S2,result_score,pvalue_t=pValue_threshold,qvalue_t=qValue_threshold,thresholdKind=thresholdKind)
                        },
        "LRT" = {
                  PairS1S2_result <- FindDiff_LRT(PairS1S2,count1=count1,count2=count2,picfile=result_pic,
                                                  pvalue_t=pValue_threshold,qvalue_t=qValue_threshold,thresholdKind=thresholdKind)
                  output_LRT(PairS1S2_result,PairS1S2,result_score,pvalue_t=pValue_threshold,qvalue_t=qValue_threshold,thresholdKind=thresholdKind)
                        },
        "FET" = {
                  PairS1S2_result <- FindDiff_FET(PairS1S2,count1=count1,count2=count2,picfile=result_pic,
                                                  pvalue_t=pValue_threshold,qvalue_t=qValue_threshold,thresholdKind=thresholdKind)
                  output_FET(PairS1S2_result,PairS1S2,result_score,pvalue_t=pValue_threshold,qvalue_t=qValue_threshold,thresholdKind=thresholdKind)
                        },
        "CTR" = {
                  Check_TR(PairS1S2_norm,count1=count1,count2=count2,picfile=result_pic,pvalue_t=pValue_threshold)
                    }

      )
 DEV_OFF((outputDir != "none"), dev_cur=dev_cur)
 cat("output ...\n")
 flush.console();



  
#######################################################################
#########generate the XHTML file ######################################
#######################################################################
 if(outputDir != "none"){
    output_html <- paste(outputDir,"/output.html",sep="")
 }else{
   output_html <- "none"
 }


output2html <- function(content,file=output_html,append=TRUE){
   if(file != "none"){
      cat(content,"\n",file=output_html,append=append)
   }
}


output_v <- function(parray){
   tmp <- parray
   for(i in (1:length(tmp))){
     #cat(paste("\"", tmp[i], "\""))
     output2html(paste("\"", tmp[i], "\"",sep=""));
     output2html("&nbsp");
   }
}

   type <- getOption("pkgType")
   if(type=="win.binary"){
      width <- 480
      heght <- 480
   }else{
      width <- 576
      heght <- 576
   }

output2html(" <!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\"> ",append=FALSE)
output2html(" <html xmlns=\"http://www.w3.org/1999/xhtml\" lang=\"en\" xml:lang=\"en\"> ")
output2html("  <head>                                                ")
output2html(" 	<meta http-equiv=\"Content-Type\" content=\"text/html; charset=ISO-8859-1\" /> ")
output2html("    <title> DEGseq RESULT </title>                    ")		 
output2html("       <style type=\"text/css\">                        ")
output2html("            body{                                       ")
output2html("                 background-color: #FFFFFF;             ")
output2html("                 margin-left: 5%;                       ")
output2html("                 margin-right: 45%;                      ")
output2html("                 border: 0px dotted gray;               ")
output2html("                 padding: 10px 10px 10px 10px;          ")
output2html("                 font-family: sans-serif;               ")
output2html("                 font-size: 130%;                       ")
output2html("            }                                           ")
output2html("            img{                                        ")
output2html(paste("              width: ",width,"px;",sep=""))
output2html(paste("              heght: ",heght,"px;",sep=""))
output2html("            }                                           ")
output2html("            td{text-align:center;}                      ")
output2html("            p.tips{                                     ")
output2html("                color: green;                           ")
output2html("                font-size: 130%;                        ")
output2html("                font-family: Trebuchet MS;              ")
output2html("                width: 200%;                            ")
output2html("            }                                           ")
output2html("            p.output{                                     ")
output2html("                font-size: 100%;                        ")
output2html("                width: 200%;                            ")
output2html("            }                                           ")
output2html("      </style>                                          ")
output2html("  </head>                                               ")
output2html("  <body>                                                      ")
output2html("      <p class = \"tips\"> Tips: change zoom level with Ctrl + or Ctrl - </p>")
output2html("      <table cellpadding=\"10\">                                ")       
output2html("      <tr>                                                      ")
output2html("        <td> <img src = \"./output/Sample1_hist.png\" alt=\"Sample1_hist.png\" /></td>    ")            
output2html("        <td> <img src = \"./output/Sample2_hist.png\" alt=\"Sample2_hist.png\" /></td>    ")          
output2html("      </tr>                                                      ")                                   
output2html("      <tr>                                                       ") 
output2html("         <td>&nbsp;&nbsp;&nbsp;&nbsp; Histogram of the number of reads for genes</td>    ")   
output2html("         <td>&nbsp;&nbsp;&nbsp;&nbsp; Histogram of the number of reads for genes</td>            ")
output2html("      </tr>                                                      ") 
output2html("      </table>                                                   ") 
output2html("      <p> <br/> </p>                                             ") 
output2html("      <table cellpadding=\"10\">                                   ")
output2html("      <tr>                                                       ")                   
output2html("        <td> <img src = \"./output/SampleS_box.png\" alt=\"SampleS_box.png\" /> </td>")           
output2html("        <td> <img src = \"./output/Sample1_Sample2_compare.png\" alt=\"Sample1_Sample2_compare.png\" /> </td>")
output2html("      </tr>                                                      ")    
output2html("      <tr>                                                       ") 
output2html("      <td>&nbsp;&nbsp;&nbsp;&nbsp; Boxplot of read counts for each gene</td>                  ")      
output2html(paste("<td>Scatterplot comparing the number of reads for each gene for ",label1,"and",label2,"</td>"))  
output2html("      </tr>                                                      ") 
output2html("      </table>                                                   ") 
output2html("      <p> <br/> </p>                                             ") 

if(normalMethod!="none"){      
output2html("      <table cellpadding=\"10\">                                   ")
output2html("      <tr>                                                       ")       
output2html("        <td> <img src = \"./output/pre_norm.png\" alt=\"pre_norm.png\" /> </td>")
output2html("        <td> <img src = \"./output/after_norm.png\" alt=\"after_norm.png\" /> </td> ")
output2html("      </tr>                                                      ")                 
output2html("      <tr>                                                       ")         
output2html("         <td>&nbsp;&nbsp;&nbsp;&nbsp; MA-plot before normalization</td>              ") 
output2html("         <td>&nbsp;&nbsp;&nbsp;&nbsp; MA-plot after normalization</td>               ")  
output2html("      </tr>                                                      ") 
output2html("      </table>                                                   ")   
output2html("      <p> <br/> </p>                                             ")
}
output2html("      <table cellpadding=\"10\">                                   ")                            
output2html("      <tr>                                                       ")
output2html("        <td> <img src = \"./output/result.png\" alt=\"result.png\"  /> </td>")
if(method == "MATR"){
output2html("        <td> <img src = \"./output/control1_control2_compare.png\" alt=\"control1_control2_compare.png\" /> </td>")
}
output2html("      </tr>                                                      ")              
output2html("      <tr>                                                       ")    
if(method != "CTR"){
output2html("       <td>&nbsp;&nbsp;&nbsp;&nbsp; Differentially expressed genes on the MA-plot</td>       ")
}else{
output2html("       <td>&nbsp;&nbsp;&nbsp;&nbsp; The variation between technical replicates</td>       ")
}
if(method == "MATR"){
output2html(paste("<td>Scatterplot comparing the number of reads for each gene for ",c_label1,"and",c_label2,"</td>"))
}
output2html("      </tr>                                                      ")
output2html("      </table>                                                   ")


output2html("<p class = \"output\">")

output2html("<br/>")
output2html("<br/>")

output2html(paste("geneExpFile1:",geneExpFile1))
output2html("<br/>")

output2html(paste("geneExpFile2:",geneExpFile2))
output2html("<br/>")

if(replicate1!="none"){
   output2html(paste("replicate1:",replicate1))
   output2html("<br/>")
   output2html(paste("replicate2:",replicate2))
   output2html("<br/>")
}


output2html(paste("method to identify differentially expressed genes:",method))


if(normalMethod !="none"){
  output2html("<br/>")
  output2html(paste("normalization method:", normalMethod))
}
output2html("<br/>")
output_str <- paste("output directory: ", "<a href = \"file:///", outputDir, "\">",outputDir, "</a>", sep="")
output2html(output_str, outputDir)

output2html("</p>")
output2html("  </body>                                                        ")
output2html("</html>                                                  ")                            
if(outputDir != "none"){
   options(warn=-1)
   # tmp <- try(open_html(path=outputDir),silent = TRUE)
   options(warn=0)
   cat("\nDone ...\n")
}
cat("The results can be observed in directory: ", outputDir, "\n")
} # main function end
