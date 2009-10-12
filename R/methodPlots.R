############################################################################
# methodPlots.R
# 
#
# S4 methods
# 
# The functions are modified from the package marray (Yang Y.H et al.).
#
# functions to generate graphs
###########################################################################


plot.lanePair <-
  function (x, xvar = "AVal", yvar = "MVal", zvar=NULL, lines.func=NULL,text.func,legend.func=NULL,...)
  {maPlot(m=x, x=xvar, y=yvar, z=zvar, lines.func=lines.func, text.func=text.func, legend.func=legend.func,...)}

plot.PairNorm <-
  function (x, xvar = "AVal", yvar = "MVal", zvar=NULL, lines.func=NULL,text.func,legend.func=NULL,...)
  {maPlot(m=x, x=xvar, y=yvar, z=zvar, lines.func=lines.func, text.func=text.func, legend.func=legend.func,...)}
##################################
plot.PairResult <- function(x, xvar = "AVal", yvar = "MVal", zvar=NULL, lines.func=NULL,text.func,legend.func=NULL,...)
{
   rs <- x
   M <- as.vector(MVal(rs))
   A <- as.vector(AVal(rs))
   z_sore <- as.vector(Zscore(rs))
   plot(A,M,pch=".",col=1,main=as.character(PairNotes(rs)[1]))
   abline(h=0,col="gray",lwd=0.5)
}
         

###boxplot
   
if( is.null(getGeneric("boxplot")))  setGeneric("boxplot")
setMethod("boxplot", signature(x="lanePair"), function (x, xvar = NULL, yvar = "MVal", ...){
            maBoxplot(m=x, x=xvar, y=yvar, ...)})
setMethod("boxplot", signature(x="PairNorm"), function (x, xvar = NULL, yvar = "MVal", ...){
            maBoxplot(m=x, x=xvar, y=yvar, ...)})
setMethod("boxplot", signature(x="laneRaw"), function (x, xvar = NULL, yvar = "LogVal", ...){
            maBoxplot(m=x, x=xvar, y=yvar, ...)})


## points, lines, text
if( is.null(getGeneric("points")))  setGeneric("points")
setMethod("points", signature(x="lanePair"), function (x, xvar = "AVal", yvar = "MVal", ...)
          {addPoints(object=x, xvar=xvar, yvar=yvar, ...)})
setMethod("points", signature(x="PairNorm"), function (x, xvar = "AVal", yvar = "MVal", ...)
          {addPoints(object=x, xvar=xvar, yvar=yvar, ...)})

if( is.null(getGeneric("text")))  setGeneric("text")
setMethod("text", signature(x="lanePair"), function (x, xvar = "AVal", yvar = "MVal",...)
          {addText(object=x, xvar=xvar, yvar=yvar, ...)})
setMethod("text", signature(x="PairNorm"), function (x, xvar = "AVal", yvar = "MVal", ...)
          {addText(object=x, xvar=xvar, yvar=yvar, ...)})

if( is.null(getGeneric("lines")))  setGeneric("lines")
setMethod("lines", signature(x="lanePair"), function (x, xvar = "AVal", yvar = "MVal", zvar=NULL,lines.func,legend.func,...)
          {addLines(object=x, xvar=xvar, yvar=yvar,lines.func=lines.func,legend.func=legend.func, ...)})
setMethod("lines", signature(x="PairNorm"), function (x, xvar = "AVal", yvar = "MVal", zvar=NULL,lines.func,legend.func,...)
          {addLines(object=x, xvar=xvar, yvar=yvar,lines.func=lines.func,legend.func=legend.func, ...)})


addText <- function(object, xvar="AVal", yvar="MVal", subset=NULL, labels = as.character(1:length(subset)), ...)
{
    text.func <- maText(subset=subset, labels=labels, ...)
    text.func(as.numeric(eval(call(xvar,object))), as.numeric(eval(call(yvar,object))))
}


addPoints <- function(object, xvar="AVal", yvar="MVal", subset=TRUE, ...)
{
    xx <- as.numeric(eval(call(xvar,object)))
    yy <- as.numeric(eval(call(yvar,object)))
    points(xx[subset], yy[subset], ...)
}

addLines <- function(object, xvar="AVal", yvar="MVal", zvar=NULL, subset=TRUE,lines.func,legend.func, ...)
{
    object <- object
    defs <- maDefaultPar(m=object,x=xvar,y=yvar,z=zvar)
    if(missing(lines.func)){
      lines.func1<-do.call("maLoessLines",
                        c(list(subset=subset, loess.args=list(span=0.2, degree=1, family="symmetric",
                                              control=loess.control(trace.hat="approximate",
                                                iterations=5,surface="direct"))),defs$def.lines))
      lines.func = list(lines.func1)
    }
    if(missing(legend.func)){
      legend.func1 <- do.call("maLegendLines",defs$def.legend)
      legend.func = list(legend.func1)
    }
    xx <- as.numeric(eval(call(xvar,object)))
    yy <- as.numeric(eval(call(yvar,object)))
    if(!is.null(zvar))
      zz <- as.numeric(eval(call(zvar,object)))
    else
      zz<-rep(1,length(xx))
      
    m <- object
    if(length(lines.func)>0)
    { 
      for(func in lines.func)
        func(xx, yy, zz)
    }
    if(length(legend.func)>0)
    { 
      for(func in legend.func)
        func(xx, yy)
    }
  }
 

###########################################################################
# Default plotting parameters for microarray objects
###########################################################################
# Compare default and ... plotting parameters and let ... overwrite defaults

maDotsDefaults<-function(dots, defaults)
{
  args<-c(dots,defaults[setdiff(names(defaults),names(dots))])
  return(args)
}

#########################
# Default parameters for microarray objects of class lanePair and PairNorm
maDefaultPar<-function(m,x,y,z)
{   
    if (is(m, "laneRaw"))
       main<-as.character(laNotes(m)[1])
    if ((is(m, "lanePair"))||(is(m, "PairNorm")))
       main<-as.character(PairNotes(m)[1])
    xlab<-ylab<-zlab<-""
    col<-2
    lty<-1
    lwd<-2.5
    las<-1
    names<-""
    def.legend<-list()
    ylab<-strsplit(y,"Val")[[1]]
    if(ylab[1]=="")
     ylab<-ylab[2]

    if(!is.null(x))
    {
      xlab<-strsplit(x,"Val")[[1]]
      if(xlab[1]=="")
        xlab<-xlab[2]
    }
    if(!is.null(z))
    { 
      zz<-eval(call(z,m))
      zlab<-strsplit(z,"Val")[[1]]
      if(zlab[1]=="")
        zlab<-zlab[2]

      if(z!="maPrintTip")
      {
        ## names<-paste(zlab,unique(zz),sep=" ")  ## BUG [modified Feb 27, 2004]
        ifelse(is.factor(zz), tmp <- levels(zz), tmp<- unique(zz))
        names<-paste(zlab, tmp, sep=" ")  
        col<-(1:length(names))
        lty<-rep(1,length(names))
        las<-1
        ncol<-1
        ord<-order(lty,col)
        def.legend<-list(legend=names[ord],col=col[ord], lty=lty[ord], lwd=lwd, ncol=ncol)
      }
   }

    def.box<-list(xlab=xlab,ylab=ylab,names=names,col=col,las=las,main=main)
    def.plot<-list(xlab=xlab,ylab=ylab,pch=".",col=1,main=main)
    def.lines<-list(col=col,lty=lty,lwd=lwd)
    def.text<-list(pch=16,col="purple")
    return(list(def.box=def.box,def.plot=def.plot,def.lines=def.lines,
                def.legend=def.legend,def.text=def.text))
}

###########################################################################
# maBoxplot: Boxplot methods
###########################################################################
# Boxplots for a single and multiple arrays

maBoxplot<- function(m, x=NULL, y="MVal", ...) {
    opt<-list(...)
    N <- 1
    if (is(m, "laneRaw"))
        N <- laNlanes(m)
    if (is(m, "lanePair")||is(m, "PairNorm"))
        N <- Npairs(m)
        
    if(N == 1)
    {
      yy<-as.numeric(eval(call(y,m)))
      if(is.null(x))
        xx<-rep(1,length(yy))
      if(!is.null(x))
        xx<-eval(call(x,m))
      # Optional graphical parameter defaults
      def<-maDefaultPar(m,x,y,x)$def.box
      if(!is.null(opt))
        def<-maDotsDefaults(opt,def)
      args<-c(list(yy~xx),def)
      do.call("boxplot",args)
    }
    if(N > 1)
    {
     yy<-as.data.frame(eval(call(y,m)))

     # Optional graphical parameter defaults
     if (is(m, "laneRaw")){
       if(length(laNotes(m)) != 0)
          def <- list(names=laNotes(m),ylab="log2(Number of reads mapped to a gene)",col=2)
       else
          def <- list(names=dimnames(yy)[[2]], ylab="log2(Number of reads mapped to a gene)",col=2)
     }
     if ((is(m, "lanePair"))||(is(m, "PairNorm"))){
       #cat("here2\n")
       if(length(PairNotes(m)) != 0)
          def <- list(names=PairNotes(m),ylab=strsplit(y,"Val")[[1]][1],col=2)
       else
          def <- list(names=dimnames(yy)[[2]], ylab=strsplit(y,"Val")[[1]][1],col=2)
     }
     #if(!is.null(opt))
     #   def<-maDotsDefaults(opt,def)
     
      args<-c(list(yy),def)
      do.call("boxplot",args)
    } 
    ##if(y=="MVal") abline(h=0,col="gray",lwd=2.5)
}

###########################################################################
##
# maPlot: Scatter-plot methods with fitted lines and points highlighted
##
###########################################################################
# General function for scatter-plot of y vs. x with fitted lines within
# values of z and subset of points highlighted

maPlot.func<-function(x, y, z,
    lines.func=NULL,
    text.func=maText(),
    legend.func=NULL,
    ...)
{
  plot(x,y,...)
  abline(h=0,col="gray",lwd=0.5)
  # Plot fitted curves
  if(length(lines.func)>0)
  { 
    for(func in lines.func)
        func(x, y, z)
  }

  # Legend
  if(length(legend.func)>0)
  { 
    for(func in legend.func)
        func(x, y)
  }

  # Label a subset of points
  if(!is.null(text.func))
    text.func(x,y)
}

#########################
# Label a subset of points
# Jean: Modify (Oct 21, 2002)  line "tmp <- length(c(1:length(subset))[subset])"

maText <-function (subset = NULL, labels = as.character(1:length(subset)),
    ...)
{
    function(x, y) {
      tmp <- length(c(1:length(subset))[subset])
      if (tmp > 0) {
            if (length(subset) < length(labels))
                text(x[subset], y[subset], labels[subset], ...)
            if (length(subset) > length(labels))
                text(x[subset], y[subset], labels, ...)
            if ((length(subset) == length(labels)) & is.logical(subset))
                text(x[subset], y[subset], labels[subset], ...)
            if ((length(subset) == length(labels)) & is.numeric(subset))
                text(x[subset], y[subset], labels, ...)
        }
    }
}


#########################
# Plot fitted lines

#  Med line

maMedLines<-function(subset=TRUE,col=2)
{ 
  function(x,y,z)
  {
    subset<-maNum2Logic(length(y), subset)
    yfit<-rep(NA,length(y))
    yfit[subset] <- median(y[subset], na.rm=TRUE)
    abline(a=median(y[subset], na.rm=TRUE),b=0,col=col)
  }
}

# Lowess
maLowessLines<-function(subset=TRUE, f=0.2, col=2, lty=1, lwd=2.5,...)
{
  function(x,y,z)
  {
    subset<-maNum2Logic(length(x), subset)
    g<-unique(z[subset])
    if(length(col)<length(g))
      col<-rep(col[1],length(g))
    if(length(lty)<length(g))
      lty<-rep(lty[1],length(g))

    for(i in (1:length(g)))
    {
      which<-z[subset]==g[i]
      xx<-x[subset][which]
      yy<-y[subset][which]
      ind <- is.na(xx) | is.na(yy) | is.infinite(xx) | is.infinite(yy)
      fit<- lowess(xx[!ind], yy[!ind], f=f)
      lines(fit,col=col[i],lty=lty[i],lwd=lwd,...)
    }
  }
}

# Loess

maLoessLines<-function(subset=TRUE, weights=NULL,
                       loess.args=list(span=0.2, degree=1, family="symmetric",
                       control=loess.control(trace.hat="approximate",
                       iterations=5,surface="direct")),col=2, lty=1, lwd=2.5, ...)
{
  function(x,y,z)
  {
    subset<-maNum2Logic(length(x), subset)
    g<-unique(z[subset])
    if(length(col)<length(g))
      col<-rep(col[1],length(g))
    if(length(lty)<length(g))
      lty<-rep(lty[1],length(g))

    for(i in (1:length(g)))
    {
      which<-z[subset]==g[i]
      if(sum(which) >= 50000){
##        print(sum(which));      print(length(which))
        tmp <- sample(1:sum(which), 50000)
        which[which][-tmp] <- FALSE 
##        print(sum(which));        print(length(which))
      }
      xx<-x[subset][which]
      yy<-y[subset][which]
      ww<-weights[subset][which]
      args<-c(list(yy ~ xx, weights=ww),loess.args)
      fit<-do.call("loess",args)
      #xf<-seq(quantile(xx,0.005,na.rm=TRUE),quantile(xx,0.995,na.rm=TRUE),length=100)
      xx <- x
      xx <- xx[order(xx)]
      ind <- is.na(xx) | is.infinite(xx) ## add
      yf<-predict(fit,xx[!ind])    ##add
      lines(xx[!ind],yf,col=col[i],lty=lty[i],lwd=lwd,...)
      #points(xx[!ind],yf,col=7)
      #cat(length(xx[!ind]),":number2\n")
    }
  }
}


#########################
# Add legend to existing plot

maLegendLines<-function(legend="", col=2, lty=1, lwd=2.5, ncol=1, ...)
{
  function(x,y)
  {
    a<-min(x[!(is.na(x)|is.infinite(x))])
    b<-max(y[!(is.na(y)|is.infinite(y))])
    legend(a,b,legend=as.character(legend),col=col,lty=lty,lwd=lwd,ncol=ncol,...)
  }
}

###########################################################################
# Methods for microarray objects: wrapper around maPlot.func
##
## Jean April 9,2003 modified default to maLoessLines
##

maPlot <- function(m, x="AVal", y="MVal", z=NULL,lines.func=NULL,text.func,legend.func=NULL, ...)
{

  # Default plotting arguments
  defs<-maDefaultPar(m,x,y,z)

  if(missing(text.func))
    text.func<-maText()


  xx<-as.numeric(eval(call(x,m)))
  yy<-as.numeric(eval(call(y,m)))
  if(is.null(z))
     zz<-rep(1,length(xx))
  if(!is.null(z))
    zz<-eval(call(z,m))

  opt<-list(...)
  if(!is.null(opt))
    def.plot<-maDotsDefaults(opt,defs$def.plot)

  do.call("maPlot.func", c(list(x=xx,y=yy,z=zz,lines.func=lines.func,text.func=text.func,legend.func=legend.func),def.plot))
  #if(y=="MVal") abline(h=0,col="gray",lwd=0.5)
}

###########################################


##############################
LineOnPlot<-function(pair,
    lineKind=list("median", "loess","Lowess"),showlegend = TRUE,...)
{
  #plot(pair)
  col <- c();
  legend <- c();
  for(oneKind in lineKind){
     switch(oneKind,
         "median" =   { lines(pair,lines.func=list(maMedLines(col=2)),legend.fun=NULL);
                   legend <- cbind(legend,c("Median"));
                   col <- cbind(col,c(2))},
         "loess" = { lines(pair,lines.func=list(maLoessLines(col=3)),legend.func=NULL);
                   legend <- cbind(legend,c("Loess"));
                   col <- cbind(col,c(3))},
         "Lowess" = { lines(pair,lines.func=list(maLowessLines(col=5)),legend.func=NULL)
                    legend <- cbind(legend,c("Lowess"));
                    col <- cbind(col,c(5))}
         )
  }
  if(showlegend)
     lines(pair,lines.func=NULL,legend.func=list(maLegendLines(col=col,legend=legend)))
}


##############################

png_new <- function(file,file_dir=""){
 # png(file=file)
 if((file_dir != "none")&&(file != "none")){ 
   if(file_dir != ""){
      file <- paste(file_dir,"/",file,sep="")
   }
   options(warn=-1)
   tmp <- try(png(file=file,units="px",width=480, height=480),silent = TRUE)
   if(class(tmp) =="try-error") {
      tmp2 <- try(bitmap(file=file, type = "png16",res=144, width=4, height=4),silent = TRUE)
      if(class(tmp2) =="try-error"){
         cat("Your system must support the function png() or bitmap()\n!")
         stop("\n")
      }else{
         #cat("Use bitmap()\n")
      }
   }else{
         #cat("Use png()\n")
   }
   par(cex.axis = 1.5,cex.main=2,cex.lab=1.5,font.axis=2,font.lab=2,lwd=3)
   options(warn=0)
 }else{
   options(warn=-1)
   ##tmp <- try(X11(),silent = TRUE)
   options(warn=0)
 }
}
