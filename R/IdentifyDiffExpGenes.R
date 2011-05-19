###############################################################
########### IdentifyDiffExpGenes.R
########### functions:
###########          FindDiff_SD      for MART and CTR
###########          FindDiff_MARS    for MARS
###########          FindDiff_LRT     for LRT
###########          FindDiff_FET     for FET
###########          FindDiff_FC      for FC
###########
########### see document for more detail 
###############################################################

getPQvalue <- function(rbatch){
   Zscore_v <- as.vector(Zscore(rbatch))
   Pvalue_v <- 2*pnorm(-abs(Zscore_v))
   n_count <- length(Pvalue_v[!is.na(Pvalue_v)])
   Qvalue_v <- Pvalue_v*n_count/rank(Pvalue_v)
   Qvalue_v[Qvalue_v > 1] <- 1
   Pvalue(rbatch) <- cbind(Pvalue_v)
   Qvalue(rbatch) <- cbind(Qvalue_v)
   rbatch
}

getQvalue1 <- function(Pvalue_v){
   n_count <- length(Pvalue_v[!is.na(Pvalue_v)])
   Qvalue_v <- Pvalue_v*n_count/rank(Pvalue_v)
   Qvalue_v[Qvalue_v > 1] <- 1
   Qvalue_v[Qvalue_v < 0] <- 0
   Qvalue_v
}

getQvalue2 <- function(Pvalue_v){
   Qvalue_v <- rep(NA,length(Pvalue_v))
   tmp <- try(library(qvalue))
   if(class(tmp) !="try-error") {
      index_v <- (1:length(Pvalue_v))
      index_v <- index_v[!is.na(Pvalue_v)]
      Pvalue_tmp <- Pvalue_v[!is.na(Pvalue_v)]
      Pvalue_tmp[Pvalue_tmp > 1] <- 1
      Pvalue_tmp[Pvalue_tmp < 0] <- 0
      tmp2 <- try({res <- qvalue(Pvalue_tmp);Qvalue_v[index_v] <- res$qvalues})
      if(class(tmp2) =="try-error"){
         cat("Warning: The default values for qvalue function are not suitable!!!\n")
         cat("Ouput qvalues(Benjamini et al. 1995) instead.\n")
         Qvalue_v  <- getQvalue1(Pvalue_v)
      }
   }else{
      cat("Warning:fail to load library qvalue(Storey et al. 2003)!\n")
      cat("Only ouput qvalue(Benjamini et al. 1995).\n")   
   }
   Qvalue_v
}


output_SD <- function(rbatch, rpairraw, file, pvalue_t=0.0001, qvalue_t=0.0001, thresholdKind=1, count1, count2){
    
  expVals1_v <- expVals1(rpairraw)
  expVals2_v <- expVals2(rpairraw)

  expVals1_v[(expVals1_v > 0) & (expVals1_v < 1)] <- 0
  expVals2_v[(expVals2_v > 0) & (expVals2_v < 1)] <- 0
  
  PairGenes_v <- PairGenes(rbatch)
  AVal_v <- AVal(rbatch)
  MVal_v <- MVal(rbatch)
  Zscore_v <- Zscore(rbatch)
  Pvalue_v <- Pvalue(rbatch)
  Qvalue_v1 <- Qvalue(rbatch)
  Qvalue_v2 <- Qvalue2(rbatch)
  
  #################################### sort
  tmp <- Zscore_v
  tmp[(!is.na(tmp))&(tmp < 0)] <- tmp[(!is.na(tmp))&(tmp < 0)]*(-1);
  order_matrix <- order(tmp, decreasing=TRUE)
  expVals1_v <- expVals1_v[order_matrix]
  expVals2_v <- expVals2_v[order_matrix]
  PairGenes_v <- PairGenes_v[order_matrix]
  MVal_v <- MVal_v[order_matrix]
  DfGenes_v <- DfGenes(rbatch)[order_matrix]
  Qvalue_v1 <- Qvalue_v1[order_matrix]
  Pvalue_v <- Pvalue_v[order_matrix]
  Qvalue_v2 <- Qvalue_v2[order_matrix]
  Zscore_v <- Zscore_v[order_matrix]
  ####################################
  order_matrix <- order(DfGenes_v, decreasing=TRUE)
  ####################################
  expVals1_v <- expVals1_v[order_matrix]
  expVals2_v <- expVals2_v[order_matrix]
  PairGenes_v <- PairGenes_v[order_matrix]
  MVal_v <- MVal_v[order_matrix]
  DfGenes_v <- DfGenes_v[order_matrix]
  Qvalue_v1 <- Qvalue_v1[order_matrix]
  Pvalue_v <- Pvalue_v[order_matrix]
  Qvalue_v2 <- Qvalue_v2[order_matrix]
  Zscore_v <- Zscore_v[order_matrix]
  #################################### sort 
  
  tmp_string <- paste("\"Signature(p-value < ",pvalue_t,")\"",sep="")
  if((thresholdKind ==1)||(thresholdKind ==2)){
      tmp_string <- paste("\"Signature(p-value < ",pvalue_t,")\"",sep="")
  }
  if(thresholdKind ==3){
      tmp_string <- paste("\"Signature(q-value(Benjamini et al. 1995) < ",qvalue_t,")\"",sep="")
  }
  if(thresholdKind ==4){
      tmp_string <- paste("\"Signature(q-value(Storey et al. 2003) < ",qvalue_t,")\"",sep="")
  }
  if(thresholdKind ==5){
     tmp_string <- paste("\"Signature(q-value(Storey et al. 2003) < ",qvalue_t,")\"",sep="")
  }
  MVal_v_n <- MVal_v+ log(count2/count1,2)
  ##################################### sorting down 
  
  tmp <- cbind(PairGenes_v, expVals1_v, expVals2_v, MVal_v,MVal_v_n, Zscore_v, Pvalue_v, Qvalue_v1, Qvalue_v2, DfGenes_v)
  dimnames(tmp) <- list(c(),c("\"GeneNames\"","\"value1\"","\"value2\"","\"log2(Fold_change)\"","\"log2(Fold_change) normalized\"",
                        "\"z-score\"","\"p-value\"","\"q-value(Benjamini et al. 1995)\"","\"q-value(Storey et al. 2003)\"",tmp_string))
  if(file != "none"){
     write.table(tmp,file=file,append=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  }
}

output_LRT <- function(rbatch, rpairraw, file, pvalue_t=0.0001, qvalue_t=0.0001, thresholdKind=1,count1, count2){
  
  expVals1_v <- expVals1(rpairraw)
  expVals2_v <- expVals2(rpairraw)

  PairGenes_v <- PairGenes(rbatch)
  AVal_v <- AVal(rbatch)
  MVal_v <- MVal(rbatch)
  #Zscore_v <- Zscore(rbatch)
  Pvalue_v <- Pvalue(rbatch)
  Qvalue_v1 <- Qvalue(rbatch)
  Qvalue_v2 <- Qvalue2(rbatch)
  
  #################################### sort 
  tmp_string <- paste("\"Signature(p-value < ",pvalue_t,")\"",sep="")
  if((thresholdKind ==1)||(thresholdKind ==2)){
      expVals1_v <- expVals1_v[order(Pvalue_v)]
      expVals2_v <- expVals2_v[order(Pvalue_v)]
      PairGenes_v <- PairGenes_v[order(Pvalue_v)]
      MVal_v <- MVal_v[order(Pvalue_v)]
      #Zscore_v <- Zscore_v[order(Pvalue_v)]
      DfGenes_v <- DfGenes(rbatch)[order(Pvalue_v)]
      Qvalue_v1 <- Qvalue_v1[order(Pvalue_v)]
      Qvalue_v2 <- Qvalue_v2[order(Pvalue_v)]
      Pvalue_v <- Pvalue_v[order(Pvalue_v)]
  }
  if(thresholdKind ==3){
      expVals1_v <- expVals1_v[order(Qvalue_v1)]
      expVals2_v <- expVals2_v[order(Qvalue_v1)]
      PairGenes_v <- PairGenes_v[order(Qvalue_v1)]
      MVal_v <- MVal_v[order(Qvalue_v1)]
      #Zscore_v <- Zscore_v[order(Qvalue_v1)]
      DfGenes_v <- DfGenes(rbatch)[order(Qvalue_v1)]
      Qvalue_v2 <- Qvalue_v2[order(Qvalue_v1)]
      Pvalue_v <- Pvalue_v[order(Qvalue_v1)]
      Qvalue_v1 <- Qvalue_v1[order(Qvalue_v1)]
      tmp_string <- paste("\"Signature(q-value(Benjamini et al. 1995) < ",qvalue_t,")\"",sep="")
  }
  if(thresholdKind ==4){
      expVals1_v <- expVals1_v[order(Qvalue_v2)]
      expVals2_v <- expVals2_v[order(Qvalue_v2)]
      PairGenes_v <- PairGenes_v[order(Qvalue_v2)]
      MVal_v <- MVal_v[order(Qvalue_v2)]
      #Zscore_v <- Zscore_v[order(Qvalue_v2)]
      DfGenes_v <- DfGenes(rbatch)[order(Qvalue_v2)]
      Qvalue_v1 <- Qvalue_v1[order(Qvalue_v2)]
      Pvalue_v <- Pvalue_v[order(Qvalue_v2)]
      Qvalue_v2 <- Qvalue_v2[order(Qvalue_v2)]
      tmp_string <- paste("\"Signature(q-value(Storey et al. 2003) < ",qvalue_t,")\"",sep="")
  }
  MVal_v_n <- MVal_v+ log(count2/count1,2)
  ##################################### sorting down 
  
  tmp <- cbind(PairGenes_v, expVals1_v, expVals2_v, MVal_v,MVal_v_n, Pvalue_v, Qvalue_v1, Qvalue_v2, DfGenes_v)
  dimnames(tmp) <- list(c(),c("\"GeneNames\"","\"value1\"","\"value2\"","\"log2(Fold_change)\"","\"log2(Fold_change) normalized\"",
                        "\"p-value\"","\"q-value(Benjamini et al. 1995)\"","\"q-value(Storey et al. 2003)\"",tmp_string))
  if(file != "none"){
     write.table(tmp,file=file,append=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  }
}

output_FET <- function(rbatch, rpairraw, file, pvalue_t=0.0001, qvalue_t=0.0001, thresholdKind=1,count1, count2){
  
  expVals1_v <- expVals1(rpairraw)
  expVals2_v <- expVals2(rpairraw)

  PairGenes_v <- PairGenes(rbatch)
  AVal_v <- AVal(rbatch)
  MVal_v <- MVal(rbatch)
  #Zscore_v <- Zscore(rbatch)
  Pvalue_v <- Pvalue(rbatch)
  Qvalue_v1 <- Qvalue(rbatch)
  Qvalue_v2 <- Qvalue2(rbatch)
  

  #################################### sort 
  tmp_string <- paste("\"Signature(p-value < ",pvalue_t,")\"",sep="")
  if((thresholdKind ==1)||(thresholdKind ==2)){
      expVals1_v <- expVals1_v[order(Pvalue_v)]
      expVals2_v <- expVals2_v[order(Pvalue_v)]
      PairGenes_v <- PairGenes_v[order(Pvalue_v)]
      MVal_v <- MVal_v[order(Pvalue_v)]
      #Zscore_v <- Zscore_v[order(Pvalue_v)]
      DfGenes_v <- DfGenes(rbatch)[order(Pvalue_v)]
      Qvalue_v1 <- Qvalue_v1[order(Pvalue_v)]
      Qvalue_v2 <- Qvalue_v2[order(Pvalue_v)]
      Pvalue_v <- Pvalue_v[order(Pvalue_v)]
  }
  if(thresholdKind ==3){
      expVals1_v <- expVals1_v[order(Qvalue_v1)]
      expVals2_v <- expVals2_v[order(Qvalue_v1)]
      PairGenes_v <- PairGenes_v[order(Qvalue_v1)]
      MVal_v <- MVal_v[order(Qvalue_v1)]
      #Zscore_v <- Zscore_v[order(Qvalue_v1)]
      DfGenes_v <- DfGenes(rbatch)[order(Qvalue_v1)]
      Qvalue_v2 <- Qvalue_v2[order(Qvalue_v1)]
      Pvalue_v <- Pvalue_v[order(Qvalue_v1)]
      Qvalue_v1 <- Qvalue_v1[order(Qvalue_v1)]
      tmp_string <- paste("\"Signature(q-value(Benjamini et al. 1995) < ",qvalue_t,")\"",sep="")
  }
  if(thresholdKind ==4){
      expVals1_v <- expVals1_v[order(Qvalue_v2)]
      expVals2_v <- expVals2_v[order(Qvalue_v2)]
      PairGenes_v <- PairGenes_v[order(Qvalue_v2)]
      MVal_v <- MVal_v[order(Qvalue_v2)]
      #Zscore_v <- Zscore_v[order(Qvalue_v2)]
      DfGenes_v <- DfGenes(rbatch)[order(Qvalue_v2)]
      Qvalue_v1 <- Qvalue_v1[order(Qvalue_v2)]
      Pvalue_v <- Pvalue_v[order(Qvalue_v2)]
      Qvalue_v2 <- Qvalue_v2[order(Qvalue_v2)]
      tmp_string <- paste("\"Signature(q-value(Storey et al. 2003) < ",qvalue_t,")\"",sep="")
  }  
  MVal_v_n <- MVal_v+ log(count2/count1,2)
  ##################################### sorting down 

  tmp <- cbind(PairGenes_v, expVals1_v, expVals2_v, MVal_v,MVal_v_n, Pvalue_v, Qvalue_v1, Qvalue_v2, DfGenes_v)
  dimnames(tmp) <- list(c(),c("\"GeneNames\"","\"value1\"","\"value2\"","\"log2(Fold_change)\"","\"log2(Fold_change) normalized\"",
                        "\"p-value\"","\"q-value(Benjamini et al. 1995)\"","\"q-value(Storey et al. 2003)\"",tmp_string))
  if(file != "none"){
     write.table(tmp,file=file,append=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  }
}


output_FC <- function(rbatch, rpairraw, file, threshold=3, count1, count2){
  expVals1_v <- expVals1(rpairraw)
  expVals2_v <- expVals2(rpairraw)
  PairGenes_v <- PairGenes(rbatch)
  AVal_v <- AVal(rbatch)
  MVal_v <- MVal(rbatch)
  MVal_v_n <- MVal_v+ log(count2/count1,2)
  tmp <- cbind(PairGenes_v,expVals1_v,expVals2_v,MVal_v,MVal_v_n,DfGenes(rbatch))
  tmp_string <- paste("\"Signature(abs(log2(Fold_change) normalized) > ", threshold,")\"",sep="")
  dimnames(tmp) <- list(c(),c("\"GeneNames\"","\"value1\"","\"value2\"","\"log2(Fold_change)\"","\"log2(Fold_change) normalized\"",tmp_string))
  if(file != "none"){
     write.table(tmp,file=file,append=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  }
}


GetMeanAndSD <-function(rbatch){   #1%~100%
  if(is(rbatch, "PairNorm")){
    rbatch<-Norm2Result(rbatch)
  }else{
    cat("Something wrong.\n")
  }
  

  M <- as.vector(MVal(rbatch))
  A <- as.vector(AVal(rbatch))
  M <- M[order(A)]
  A <- A[order(A)]
  M <- M[!is.na(A)]
  A <- A[!is.na(A)]
  min_a <- min(A)    # new add
  max_a <- max(A)
  tmp <- min(A)+(max(A)-min(A))/10
  a_1_200 <- (max(A)-min(A))/200  
  M <- M[A > tmp]
  A <- A[A > tmp]  # new add 
  M <- M[A>0]
  A <- A[A>0]
  count_per_block <- ceiling(length(M)/100)
  block_num <-1
  block_mean <- rep(0,100)
  block_sd   <- rep(0,100)
  block_a_mean <- rep(0,100) 
  block_a_start <- rep(0,100)
  block_a_end <- rep(0,100)
  val_n <- 100;
  for(i in (1:100)){
    percent_start <- i
    #percent_start <- i - 5
    if(percent_start < 1){
      percent_start <-1
    }
    percent_end <- i
    #percent_end <- i + 4
    if(percent_start > 100){
      percent_start <- 100
    }
    block_start <- (percent_start-1)*count_per_block+1
    block_end <- percent_end*count_per_block
    if(block_start > length(M)){
      val_n <- i-1
      break
    }
    if(block_end > length(M)){
       block_end <- length(M)
    }
    window <- M[block_start:block_end]
    #window_a <- A[block_start:block_end]
    a_start <- (i-1)*count_per_block+1
    a_end   <- i*count_per_block
    
    if(a_start > length(M)){
      val_n <- i-1
      break
    }
    if(a_end > length(M)){
       a_end <- length(M)
    }
    #window_a <- A[a_start:a_end]
    window_a <- A[c(a_start,a_end)]
    
    block_mean[i] <- mean(window)
    #cat(file="./log",append = TRUE,"all:",length(M),"block_start:",block_start,"block_end:",block_end,"\n")
    #cat(file="./log",append = TRUE,i,":",window,"\n");
    block_sd[i]   <- sd(window)
    block_a_mean[i] <- mean(window_a)
    block_a_start[i] <- A[a_start]
    block_a_end[i] <- A[a_end]
  }
  
  block_mean <- block_mean[(1:val_n)]
  block_sd   <- block_sd[(1:val_n)]
  block_a_mean <- block_a_mean[(1:val_n)]
  block_a_start <- block_a_start[(1:val_n)]
  block_a_end <- block_a_end[(1:val_n)]
  
  
  block_mean_2 <- block_mean
  block_sd_2 <- block_sd
  block_a_mean_2 <- block_a_mean
  
  index2 <- 0
  for(index in (1:val_n)){
     a_seq <- seq(from=block_a_start[index],to=block_a_end[index],by=a_1_200)
     #cat("from:",block_a_start[index],"to",block_a_end[index],"\n")
     for(index_j in (1:length(a_seq))){
        block_a_mean_2[index2] <- a_seq[index_j]
        block_sd_2[index2] <- block_sd[index]
        block_mean_2[index2] <- block_mean[index]
        index2 <- index2+1
        #cat(index,"\t",index2,"\t",a_seq[index_j],"\n")
     }
  }
  
  list(block_mean_2,block_sd_2,block_a_mean_2)
} 


FindDiff_SD<-function(rbatch, cbatch, pvalue_t=0.0001, threshold=4, picfile=NULL,qvalue_t=0.0001,
                      thresholdKind=1, new_window=TRUE, filter=TRUE, line_col=4, control=1){   #1%~100%
  if(is(rbatch, "PairNorm")){
    rbatch<-Norm2Result(rbatch)
  }else{
    cat("Something wrong.\n")
  }
  if(threshold <= 0){
     threshold <- 4
  }

  if(!is.null(picfile))
    png_new(picfile)
  if(new_window==TRUE){
    plot(rbatch)
  }
  M <- as.vector(MVal(rbatch))
  A <- as.vector(AVal(rbatch))
  M <- M[order(A)]
  A <- A[order(A)]
  M <- M[!is.na(A)]
  A <- A[!is.na(A)]
  
  result_list <- GetMeanAndSD(cbatch)
  block_mean <- result_list[[1]]
  block_sd   <- result_list[[2]]
  block_a_mean <- result_list[[3]]


  
  #points(block_a_mean,block_mean+threshold*block_sd,col=5)
  #points(block_a_mean,block_mean-threshold*block_sd,col=5)
 
  #cat(file="./log",append = TRUE, "valn",val_n,"\n")
  #cat(file="./log",append = TRUE, "\npoints draw done!\n")
  loess.args=list(span=0.2, degree=1, family="symmetric",
                          control=loess.control(trace.hat="approximate",
                            iterations=5,surface="direct"))
  args1<-c(list((block_mean+4*block_sd)~block_a_mean),loess.args)   ## always use 4
  args2<-c(list((block_mean-4*block_sd)~block_a_mean),loess.args)   ## always use 4
  fit1<-do.call("loess",args1)
  fit2<-do.call("loess",args2)
 
  x_line<-seq(from=min(A),to=max(A),length=100)
  y_line1<-predict(fit1,x_line)
  y_line2<-predict(fit2,x_line)
  #lines(x_line,y_line1,lty=5,col=line_col)
  #lines(x_line,y_line2,lty=5,col=line_col)
  z_sore <- rep(NA,length(as.vector(AVal(rbatch))))
  
  for(i in (1:length(as.vector(AVal(rbatch))))){
    if(new_window!=TRUE){
       break
    }
    if(is.na(as.vector(AVal(rbatch))[i]))
       next
    yf1<-predict(fit1,as.vector(AVal(rbatch))[i])
    yf2<-predict(fit2,as.vector(AVal(rbatch))[i])
    sd_v <- (yf1-yf2)/(2*4)                                        ## always use 4
    mean_v <- (yf1+yf2)/2
    if(is.na(sd_v)){
       cat("yf1:",yf1,"yf2:",yf2,"threshold:",threshold,"\n")
    }
    if(sd_v!=0){
       z_sore[i] <- (as.vector(MVal(rbatch))[i]-mean_v)/sd_v
    }   
  }

  Pvalue_v <- 2*pnorm(-abs(z_sore))
  Qvalue_v1 <- getQvalue1(Pvalue_v)
  Qvalue_v2 <- getQvalue2(Pvalue_v)
  DiffGenes <- rep(FALSE, length(as.vector(AVal(rbatch))))
  
  # Pvalue_v[is.na(Pvalue_v)] <- 1
  # Qvalue_v1[is.na(Qvalue_v1)] <- 1
  # Qvalue_v2[is.na(Qvalue_v2)] <- 1
  
  for(i in (1:length(as.vector(AVal(rbatch))))){
    if(is.na(as.vector(AVal(rbatch))[i]))
       next
    
    if((thresholdKind == 1)||(thresholdKind == 2)){
        if(Pvalue_v[i] <  pvalue_t){
           points(AVal(rbatch)[i],MVal(rbatch)[i],pch=".",col="red")
        }
    }
    if(thresholdKind == 3){
        if(Qvalue_v1[i] <  qvalue_t){
           points(AVal(rbatch)[i],MVal(rbatch)[i],pch=".",col="red")
        }
    }
    if(thresholdKind == 4){
        if(Qvalue_v2[i] <  qvalue_t){
           points(AVal(rbatch)[i],MVal(rbatch)[i],pch=".",col="red")
        }
    }
    
  }
  
  if((thresholdKind == 1)||(thresholdKind == 2)){
     DiffGenes <- Pvalue_v < pvalue_t
  }
  if(thresholdKind == 3){
     DiffGenes <- Qvalue_v1 < qvalue_t
  }
  if(thresholdKind == 4){
     DiffGenes <- Qvalue_v2 < qvalue_t
  }
  DiffGenes[is.na(DiffGenes)] <- FALSE
  
  if((thresholdKind == 1)||(thresholdKind == 2)){
      lines(x_line,y_line1,lty=5,col=line_col)
      lines(x_line,y_line2,lty=5,col=line_col)			
  }
 
  
  #points(block_a_mean,block_mean+threshold*block_sd,col=3)
  #points(block_a_mean,block_mean-threshold*block_sd,col=3)
  if(control == 1){
     points(AVal(cbatch),MVal(cbatch),pch=".",col=4)
  }
  if(!is.null(picfile))
    DEV_OFF();

  Zscore(rbatch) <- cbind(z_sore)
  DfGenes(rbatch) <- cbind(DiffGenes)
  Pvalue(rbatch) <- cbind(Pvalue_v)
  Qvalue(rbatch) <- cbind(Qvalue_v1)
  Qvalue2(rbatch) <- cbind(Qvalue_v2)

  rbatch
} 

FindDiff_MARS <-function(rbatch, count1=0, count2=0, threshold=2, pvalue_t=0.0001,
                         qvalue_t=0.0001, thresholdKind=1, picfile=NULL, new_window=TRUE, dev_off=TRUE){
  if(is(rbatch, "PairNorm")){
    rbatch<-Norm2Result(rbatch)
  }else{
    cat("Something wrong.\n")
  }
  
  if(threshold <= 0){
     threshold <- 2
  }
  
  if(!is.null(picfile))
    png_new(picfile)
  if(new_window == TRUE){
     plot(rbatch)
  }
  A <- as.vector(AVal(rbatch))
  A <- A[!is.na(A)]
  
  seq_a <- seq(min(A), max(A)+5, length=5000)
  p_seq <- 2^seq_a/((count1*count2)^0.5)
  
  #p_max <- 2^max(A)/((count1*count2)^0.5)
  #count_max <- p_max*(count1+count2)+100 
  #p_v <- (1:count_max)/(count1+count2)
  #p_v <- seq(1/(count1+count2), count_max/(count1+count2), by=0.001)

  p_v <- p_seq
  UA <- (log(count1*p_v)+log(count2*p_v))/(2*log(2))
  EM <- (log(count1*p_v)-log(count2*p_v))/log(2)
  SD <- (4*(1-p_v)/((count1+count2)*p_v))^0.5/log(2)
  #lines(UA,EM+threshold*SD,col=2,lty=5)
  #lines(UA,EM-threshold*SD,col=2,lty=5)
  #lines(UA,EM,col=2,lty=3)
  z_sore <- rep(NA,length(as.vector(AVal(rbatch))))
  
  for(i in (1:length(as.vector(AVal(rbatch))))){
    if(is.na(as.vector(AVal(rbatch))[i]))
       next
    a <- as.vector(AVal(rbatch))[i]
    p_i <- 2^a/((count1*count2)^0.5)
    sd_i <- (4*(1-p_i)/((count1+count2)*p_i))^0.5/log(2)
    mean_i <- (log(count1*p_i)-log(count2*p_i))/log(2)
    #cat("a:",a,"p_i:",p_i,"sd_i:",sd_i,"\n")
    if(sd_i != 0){
       z_sore[i] <- (as.vector(MVal(rbatch))[i]-mean_i)/sd_i
    }
  }

  Pvalue_v <- 2*pnorm(-abs(z_sore))
  Qvalue_v1 <- getQvalue1(Pvalue_v)
  Qvalue_v2 <- getQvalue2(Pvalue_v)
  DiffGenes <- rep(FALSE, length(as.vector(AVal(rbatch))))
  
  # Pvalue_v[is.na(Pvalue_v)] <- 1
  # Qvalue_v1[is.na(Qvalue_v1)] <- 1
  # Qvalue_v2[is.na(Qvalue_v2)] <- 1
  
  for(i in (1:length(as.vector(AVal(rbatch))))){
    if(is.na(as.vector(AVal(rbatch))[i]))
       next
    if(dev_off == FALSE){
       break
    }
    if((thresholdKind == 1)||(thresholdKind == 2)){
        if(Pvalue_v[i] <  pvalue_t){
           points(AVal(rbatch)[i],MVal(rbatch)[i],pch=".",col="red")
        }
    }
    if(thresholdKind == 3){
        if(Qvalue_v1[i] <  qvalue_t){
           points(AVal(rbatch)[i],MVal(rbatch)[i],pch=".",col="red")
        }
    }
    if(thresholdKind == 4){
        if(Qvalue_v2[i] <  qvalue_t){
           points(AVal(rbatch)[i],MVal(rbatch)[i],pch=".",col="red")
        }
    }
    if(thresholdKind == 5){
        if((Qvalue_v2[i] <  qvalue_t)&&(abs(MVal(rbatch)[i]+log(count2/count1,2)) > threshold)&&(1)&&(1)){
           points(AVal(rbatch)[i],MVal(rbatch)[i],pch=".",col="red")
           points(AVal(rbatch)[i],MVal(rbatch)[i],pch="+",col="red")
        }else{
           if(Qvalue_v2[i] <  qvalue_t){
                points(AVal(rbatch)[i],MVal(rbatch)[i],pch=".",col=6)
           }
        }
    }
  }
  if((thresholdKind == 1)||(thresholdKind == 2)){
     DiffGenes <- Pvalue_v < pvalue_t
  }
  if(thresholdKind == 3){
     DiffGenes <- Qvalue_v1 < qvalue_t
  }
  if(thresholdKind == 4){
     DiffGenes <- Qvalue_v2 < qvalue_t
  }
  if(thresholdKind == 5){
     DiffGenes <- (Qvalue_v2 < qvalue_t)&(abs(MVal(rbatch)+log(count2/count1,2)) > threshold)
     abline(threshold-log(count2/count1,2),0,col=8,lwd=2)
     abline(-threshold-log(count2/count1,2),0,col=8,lwd=2)
  }
  DiffGenes[is.na(DiffGenes)] <- FALSE  
  if((thresholdKind == 1)||(thresholdKind == 2)){
      lines(UA,EM+threshold*SD,col=2,lty=5)
      lines(UA,EM-threshold*SD,col=2,lty=5)
  }
  
  
  if((!is.null(picfile))&&(dev_off==TRUE))
    DEV_OFF();
  
  Zscore(rbatch) <- cbind(z_sore)
  DfGenes(rbatch) <- cbind(DiffGenes)
  Pvalue(rbatch) <- cbind(Pvalue_v)
  Qvalue(rbatch) <- cbind(Qvalue_v1)
  Qvalue2(rbatch) <- cbind(Qvalue_v2)

  rbatch
}


FindDiff_LRT <- function(rbatch, count1=0, count2=0, picfile=NULL, pvalue_t=0.0001, 
                         qvalue_t=0.0001, thresholdKind=1, new_window=TRUE){
  
  if(!(is(rbatch, "lanePair"))){
     cat("Something wrong in funciton Ptest!\n")
  }
  
  
  if(!is.null(picfile))
    png_new(picfile)
  
  if(new_window==TRUE){
    plot(rbatch)
  } 
  value1 <- expVals1(rbatch)
  value2 <- expVals2(rbatch)
  valueA <- AVal(rbatch)
  valueM <- MVal(rbatch)
  
  p <- rep(NA,length(value1))

  count1 <- ceiling(count1)
  count2 <- ceiling(count2)
  for(i in (1:length(value1))){
     ob1 <- value1[i]
     ob2 <- value2[i]
     
     if(is.na(ob1)||is.na(ob2)){
       p[i] <- NA
       next
     }
     if((ob1 == 0)&&(ob2 == 0)){
       p[i] <- NA
       next
     }
     ob1 <- ceiling(ob1)
     ob2 <- ceiling(ob2)
     lamda1 <- ob1
     lamda2 <- ob2
     elamda1<-(ob1+ob2)*count1/(count1+count2)
     elamda2<-(ob1+ob2)*count2/(count1+count2)
     x2<- -2*(log(dpois(ob1,elamda1))+log(dpois(ob2,elamda2))-log(dpois(ob1,lamda1))-log(dpois(ob2,lamda2)))
     p[i]<- 1-pchisq(x2,1)
     if(is.na(p[i])){
        cat("ob1:",ob1,"ob2:",ob2,"\n")
     }
  }
  
  p[p>1] <- 1
  p[p<0] <- 0

  Pvalue_v <- rep(NA,length(AVal(rbatch)))
  Pvalue_v <- p
  
  Qvalue_v1 <- getQvalue1(Pvalue_v)
  Qvalue_v2 <- getQvalue2(Pvalue_v)
    
  # Pvalue_v[is.na(Pvalue_v)] <- 1
  # Qvalue_v1[is.na(Qvalue_v1)] <- 1
  # Qvalue_v2[is.na(Qvalue_v2)] <- 1
  
  rbatch<-Pair2Norm(rbatch)
  rbatch<-Norm2Result(rbatch)
  
  DiffGenes <- rep(FALSE, length(as.vector(AVal(rbatch))))
  for(i in (1:length(as.vector(AVal(rbatch))))){
      if(is.na(as.vector(AVal(rbatch))[i])){
         next
      }
      if((thresholdKind == 1)||(thresholdKind == 2)){
          if(Pvalue_v[i] <  pvalue_t){
             points(AVal(rbatch)[i],MVal(rbatch)[i],pch=".",col="red")
          }
      }
      if(thresholdKind == 3){
          if(Qvalue_v1[i] <  qvalue_t){
             points(AVal(rbatch)[i],MVal(rbatch)[i],pch=".",col="red")
          }
      }
      if(thresholdKind == 4){
          if(Qvalue_v2[i] <  qvalue_t){
             points(AVal(rbatch)[i],MVal(rbatch)[i],pch=".",col="red")
          }
      }
  }
  if((thresholdKind == 1)||(thresholdKind == 2)){
     DiffGenes <- Pvalue_v < pvalue_t
  }
  if(thresholdKind == 3){
     DiffGenes <- Qvalue_v1 < qvalue_t
  }
  if(thresholdKind == 4){
     DiffGenes <- Qvalue_v2 < qvalue_t
  }
  DiffGenes[is.na(DiffGenes)] <- FALSE
 
  
  Pvalue(rbatch) <- cbind(Pvalue_v)
  Qvalue(rbatch) <- cbind(Qvalue_v1)
  Qvalue2(rbatch) <- cbind(Qvalue_v2)
  DfGenes(rbatch) <- cbind(DiffGenes)  
  rbatch
}

FindDiff_FET <- function(rbatch, count1=0, count2=0, pvalue_t=0.001, 
                         qvalue_t=0.0001, thresholdKind=1, new_window=TRUE, picfile=NULL){
  
  if(!(is(rbatch, "lanePair"))){
     cat("Something wrong in funciton Ptest!\n")
  }
  
  
  if(!is.null(picfile))
    png_new(picfile)
  
  if(new_window == TRUE){
    plot(rbatch)
  } 
  value1 <- expVals1(rbatch)
  value2 <- expVals2(rbatch)
  valueA <- AVal(rbatch)
  valueM <- MVal(rbatch)

  Pvalue_v <- rep(NA,length(valueA))
  Qvalue_v <- rep(NA,length(valueA))
  Significant_v <- rep(NA,length(valueA))  

  count1 <- ceiling(count1)
  count2 <- ceiling(count2)
  for(i in (1:length(as.vector(value1)))){
      if(is.na(value1[i])){
         next
      }
      if(is.na(value2[i])){
         next
      }
      o1 <- ceiling(value1[i])
      o2 <- ceiling(value2[i])
      alleles <- matrix(c(o1,o2,count1-o1,count2-o2), nr=2)
      fisher_result <- fisher.test(alleles)
      Pvalue_v[i] <- fisher_result$p.value
  }

  Qvalue_v1 <- getQvalue1(Pvalue_v)
  Qvalue_v2 <- getQvalue2(Pvalue_v)
  
  # Pvalue_v[is.na(Pvalue_v)] <- 1
  # Qvalue_v1[is.na(Qvalue_v1)] <- 1
  # Qvalue_v2[is.na(Qvalue_v2)] <- 1
  
  rbatch<-Pair2Norm(rbatch)
  rbatch<-Norm2Result(rbatch)

  DiffGenes <- rep(FALSE, length(as.vector(AVal(rbatch))))
  for(i in (1:length(as.vector(AVal(rbatch))))){
      if(is.na(as.vector(AVal(rbatch))[i])){
         next
      }
      if((thresholdKind == 1)||(thresholdKind == 2)){
          if(Pvalue_v[i] <  pvalue_t){
             points(AVal(rbatch)[i],MVal(rbatch)[i],pch=".",col="red")
          }
      }
      if(thresholdKind == 3){
          if(Qvalue_v1[i] <  qvalue_t){
             points(AVal(rbatch)[i],MVal(rbatch)[i],pch=".",col="red")
          }
      }
      if(thresholdKind == 4){
          if(Qvalue_v2[i] <  qvalue_t){
             points(AVal(rbatch)[i],MVal(rbatch)[i],pch=".",col="red")
          }
      }
  }
  if((thresholdKind == 1)||(thresholdKind == 2)){
     DiffGenes <- Pvalue_v < pvalue_t
  }
  if(thresholdKind == 3){
     DiffGenes <- Qvalue_v1 < qvalue_t
  }
  if(thresholdKind == 4){
     DiffGenes <- Qvalue_v2 < qvalue_t
  }
  DiffGenes[is.na(DiffGenes)] <- FALSE
  Pvalue(rbatch) <- cbind(Pvalue_v)
  Qvalue(rbatch) <- cbind(Qvalue_v1)
  Qvalue2(rbatch) <- cbind(Qvalue_v2)
  DfGenes(rbatch) <- cbind(DiffGenes)
  rbatch
}


FindDiff_FC<-function(rbatch, threshold=3, picfile=NULL, count1, count2){   
  if(is(rbatch, "PairNorm")){
    rbatch<-Norm2Result(rbatch)
  }else{
    cat("Something wrong.\n")
  }
  
  M <- as.vector(MVal(rbatch))
  A <- as.vector(AVal(rbatch))
  z_sore <- rep(NA,length(A))
  if(!is.null(picfile))
    png_new(picfile)
  plot(rbatch)
  for(i in (1:length(A))){
    if(is.na(A[i]))
       next
    z_sore[i] <- 0
    if(abs(M[i]+log(count2/count1,2)) > threshold){
       points(A[i],M[i],pch=".",col="red")
    }
  }
  abline(h=threshold-log(count2/count1,2),col="red",lwd=1,lty=3)
  abline(h=-threshold+log(count2/count1,2),col="red",lwd=1,lty=3)
  DiffGenes <- (abs(M+log(count2/count1,2)) > threshold)
  DiffGenes[is.na(DiffGenes)] <- FALSE
  Zscore(rbatch) <- cbind(z_sore)
  DfGenes(rbatch) <- cbind(DiffGenes)
  if(!is.null(picfile))
    DEV_OFF();
  rbatch
} 

Check_TR <- function(PairS1S2_norm, count1=0, count2=0, picfile, pvalue_t){
     pValue_t <- 2*pnorm(-4)
     FindDiff_MARS(PairS1S2_norm, count1, count2, 4, pvalue_t=pvalue_t, picfile=picfile, new_window=TRUE, dev_off=FALSE)
     FindDiff_CTR(PairS1S2_norm, PairS1S2_norm, 4)
     DEV_OFF();
}


FindDiff_CTR<-function(rbatch, cbatch, threshold=4, picfile=NULL, new_window=TRUE,
                       filter=TRUE, line_col=2, control=1){   #1%~100%
  if(is(rbatch, "PairNorm")){
     rbatch<-Norm2Result(rbatch)
  }else{
     cat("Something wrong.\n")
  }
  
  M <- as.vector(MVal(rbatch))
  A <- as.vector(AVal(rbatch))
  M <- M[order(A)]
  A <- A[order(A)]
  M <- M[!is.na(A)]
  A <- A[!is.na(A)]
  result_list <- GetMeanAndSD(cbatch)
  block_mean <- result_list[[1]]
  block_sd   <- result_list[[2]]
  block_a_mean <- result_list[[3]]
  loess.args=list(span=0.2, degree=1, family="symmetric",
                  control=loess.control(trace.hat="approximate",
          iterations=5,surface="direct"))
  args1<-c(list((block_mean+threshold*block_sd)~block_a_mean),loess.args)
  args2<-c(list((block_mean-threshold*block_sd)~block_a_mean),loess.args)
  fit1<-do.call("loess",args1)
  fit2<-do.call("loess",args2)

  x_line<-seq(from=min(A),to=max(A),length=100)
  y_line1<-predict(fit1,x_line)
  y_line2<-predict(fit2,x_line)
  lines(x_line,y_line1,lty=6,col=4)
  lines(x_line,y_line2,lty=6,col=4)
}
