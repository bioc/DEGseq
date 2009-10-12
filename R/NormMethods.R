###########################################################################
# NormMethods.R
#
# 
# The functions are modified from the package marray (Yang Y.H et al.).
#
# functions to do normalization
###########################################################################
###########################################################################
# Convert integer vector of indices to logical vector

maNum2Logic<-function(n=length(subset), subset=TRUE)
{
  if(is.logical(subset))
    return(subset)
  if(is.numeric(subset))
  {
    which<-rep(FALSE,n)
    which[subset]<-TRUE
    return(which)
  }
}

maDotsMatch <- function(dots, defaults)
  {
    ind <- intersect(names(dots), setdiff(names(defaults), "..."))
    for(i in ind)
      defaults[[i]] <- dots[[i]]
    return(defaults[!names(defaults)=="..."])
  }



###########################################################################
############################################################################
# maNormMain
# Main within-slide location and scale normalization function

maNormMain<-function(pbatch,
                     f.loc=list(maNormLoess()),
                     f.scale=NULL,
                     a.loc=maCompNormEq(),
                     a.scale=maCompNormEq(),
                     Mloc=TRUE, Mscale=TRUE, echo=FALSE)
{
    if (is(pbatch, "lanePair"))
    {
        mnorm<-Pair2Norm(pbatch)
        M<-Ml<-Ms<-NULL
    }
    if (is(pbatch, "PairNorm"))
        mnorm<-pbatch

  slot(mnorm, "PairNormCall") <- match.call()

  if(length(f.loc)>0)
    M<-Ml<-NULL
  if(length(f.scale)>0)
    M<-Ms<-NULL

  for(i in 1:ncol(MVal(pbatch)))
  {

    if(echo)
      cat(paste("Normalizing array ", i, ".\n", sep=""))
    m<-pbatch
    
    M1<-M2<-NULL
    Mnorm<-MVal(m)

    # Location
    if(length(f.loc)>0)
    { 
      for(func in f.loc)
        M1<-cbind(M1, func(m))
      #Mtmp <- M1  ##
      if(length(f.loc)>1){
          M1 <- rowSums(M1*a.loc(AVal(m),length(f.loc)))
      }
      Ml<-cbind(Ml,M1)
      Mnorm<-(Mnorm - M1)
      #Ml<-cbind(Ml,Mtmp)  ##
    }


    # Scale
    if(length(f.scale)>0)
    {
      # Scale computed for location normalized data
      m<-mnorm
      slot(m, "MVal")<-Mnorm
      for(func in f.scale)
        M2<-cbind(M2, func(m))

      if(length(f.scale)>1)
          M2 <- rowSums(M2*a.scale(AVal(m),length(f.scale)))
      Ms<-cbind(Ms,M2)
      Mnorm<-Mnorm/M2
    }

    M<-cbind(M, Mnorm)
  }

  slot(mnorm, "MVal")<-M
  if(length(f.loc)>0 & Mloc)
    slot(mnorm, "pMloc")<-Ml
  if(length(f.scale)>0 & Mscale)
    slot(mnorm, "pMscale")<-Ms

  return(mnorm)
}

############################################################################
# maNorm
# Wrapper function around maNormMain for simple normalization procedures
#
# From Gordon's Mar 17 E-mail
# loess(y~x,span=0.2,degree=1,family="symmetric",control=loess.control(trace.hat="ap
# proximate",iterations=5,surface="direct"))
#

maNorm<-function(pbatch,
    norm=c("none", "median", "loess"),
    subset=TRUE, span=0.2, Mloc=TRUE, Mscale=TRUE, echo=FALSE, ...)
{
  ## Set normalization defaults
  opt <- list(...)

  ## norm <- unlist(strsplit(norm, split=""))[1]
  norm.method <- match.arg(norm)
  if(echo)
    cat(paste("Normalization method: ", norm.method, ".\n", sep=""))

  switch(norm.method,
         "none" = maNormMain(pbatch, f.loc=NULL,
           Mloc=Mloc, Mscale=Mscale, echo=echo),

         "median" = maNormMain(pbatch,
           f.loc=list(maNormMed(x=NULL, y="MVal", subset=subset)),
           Mloc=Mloc, Mscale=Mscale, echo=echo),

         "loess" = maNormMain(pbatch,
           f.loc=list(maNormLoess(x="AVal", y="MVal", z=NULL, w=NULL,
             subset=subset, span=span, ...)),
           Mloc=Mloc, Mscale=Mscale, echo=echo)
    )
}


###########################################################################
# Median location normalization
###########################################################################
# maMed
# General function: median of y within values of x

maMed<-function(x, y, subset=TRUE)
{
  subset<-maNum2Logic(length(y), subset)
  yfit<-rep(NA,length(y))
  for(i in unique(x))
    yfit[x==i]<-median(y[(x==i)&subset], na.rm=TRUE)
  return(yfit)
}

#########################
# maNormMed
# Function for objects of class marrayRaw

maNormMed<-function(x=NULL, y="MVal", subset=TRUE)
{
  function(m)
  {
    if(is.character(x))
      maMed(x=eval(call(x,m)), y=eval(call(y,m)), subset=subset)
    else
      maMed(x=TRUE, y=eval(call(y,m)), subset=subset)
  }
}

###########################################################################
# Loess location normalization
###########################################################################
# maLoess
# General function: regress y on x within values of z
# In our case, x=NA iff y=NA
## Jean modify Nov 5th, April 9, 2003
## Very ugly fix
## Jean modify Sep 14, 2004
## Sample Loess if the number is too large

maLoess <-
function (x, y, z, w = NULL, subset = TRUE, span = 0.2, ...)
{
  opt <- list(...)
  subset <- maNum2Logic(length(y), subset)
  yfit <- rep(NA, length(y))
  good <- !(is.infinite(x) | is.infinite(y) | is.na(x) | is.na(y))
  if(length(z) == 1) info <- z
  else
    info <- z[subset]
  for (i in unique(info))
    {
      which <- z == i
      ##      cat(i, " :: ", sum(which), "::", unique(z), "\n")
      if((sum(which) == 1) | (length(unique(z)) == 1)) {
        samplesub <- rep(TRUE, length(x))
        if(length(x[subset&good]) > 50000){
           samplesub[-sample(1:length(x), min(length(x), 50000))] <- FALSE
           #cat("in if ", sum(samplesub), " :: ", length(samplesub), "\n")
        }
      }
      else
        samplesub <- TRUE

      defs <- list(degree=1,
                   family="symmetric",
                   control=loess.control(trace.hat="approximate",iterations=5, surface="direct"),
                   subset = samplesub & which & subset & good,
                   span=span,
                   na.action = na.omit,
                   weights=w)
      args <- maDotsMatch(c(defs, opt), formals(args("loess")))
      #cat(paste("span: ", args$span, ".\n", sep=""))
      #cat(paste("length(y): ", length(y), ".\n", sep=""))
      #cat(paste("length(x): ", length(x), ".\n", sep=""))
      fit <- loess(y ~ x,  weights = args$w, subset = args$subset,
                   span = args$span, na.action = args$na.action,
                   degree = args$degree, family=args$family,
                   control=args$control)
      ##      gc(TRUE)
      ##      fit <- loess(y ~ x, weights = w, subset = which & subset &
      ##                   good, span = span, na.action = na.omit)
      yfit[which & good] <- predict(fit, x[which & good])
      #cat("here2\n");
    }
    #cat(length(yfit[which & good]),":number2\n")
    return(yfit)
}


#########################
# maNormLoess


maNormLoess<-function(x="AVal", y="MVal", z=NULL, w=NULL, subset=TRUE, span=0.2, ...)
{
  function(m)
  {
    if(is.character(z))
      maLoess(x=eval(call(x,m)), y=eval(call(y,m)), z=eval(call(z,m)), w=w, subset=subset, span=span, ...)
    else
      maLoess(x=eval(call(x,m)), y=eval(call(y,m)), z=TRUE, w=w, subset=subset, span=span, ...)
  }
}
########################

###########################################################################
# Weights for composite normalization
###########################################################################
# maCompNormEq
# Equal weights for composite normalization
maCompNormEq<-function()
{
  function(x,n)
  {
    return(matrix(1/n,length(x),n))
  }
}

#########################
# maCompNormA
# Weights for composite intensity dependent normalization

maCompNormA<-function()
{
  function(x,n)
  {
    if(n!=2)
      stop("Can only do composite location normalization for 2 functions in f.loc")
    f<-ecdf(x)
    a<-f(x)
    return(as.matrix(cbind(a,1-a)))
  }
}

###########################################################################
