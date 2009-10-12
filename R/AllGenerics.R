####################################
#########  AllGenerics.R
#########
#########      all generics in DEGseq
####################################
###########################################################################


###########################################################################
# laneRaw

# Accessor methods for laneRaw class

if(!isGeneric("expVals"))
   setGeneric("expVals", function(object) standardGeneric("expVals"))
setMethod("expVals", "laneRaw", function(object) slot(object, "expVals"))

if(!isGeneric("laGenes"))
   setGeneric("laGenes", function(object) standardGeneric("laGenes"))
setMethod("laGenes", "laneRaw", function(object) slot(object, "laGenes"))

if(!isGeneric("laNotes"))
   setGeneric("laNotes", function(object) standardGeneric("laNotes"))
setMethod("laNotes", "laneRaw", function(object) slot(object, "laNotes"))

if(!isGeneric("laNlanes"))
   setGeneric("laNlanes", function(object) standardGeneric("laNlanes"))
setMethod("laNlanes", "laneRaw", function(object) ncol(expVals(object)))

if(!isGeneric("LogVal"))
   setGeneric("LogVal", function(object) standardGeneric("LogVal"))
setMethod("LogVal", "laneRaw",
   function(object){
     r<-expVals(object)
     if(length(expVals(object))!=0)
       {
        r<-log(ifelse(r>0, r, NA),2)
       }
     r
   }
)

######################################
# Assignment methods for laneRaw class

if( !isGeneric("expVals<-") )
      setGeneric("expVals<-", function(object, value)
               standardGeneric("expVals<-"))

setReplaceMethod("expVals", signature("laneRaw", "matrix"),
  function(object, value) {
     slot(object,"expVals")<- value
     object
  })

if( !isGeneric("laGenes<-") )
      setGeneric("laGenes<-", function(object, value)
               standardGeneric("laGenes<-"))

setReplaceMethod("laGenes", signature("laneRaw", "matrix"),
  function(object, value) {
     slot(object,"laGenes")<- value
     object
  })

if( !isGeneric("laNotes<-") )
      setGeneric("laNotes<-", function(object, value)
               standardGeneric("laNotes<-"))

setReplaceMethod("laNotes", signature("laneRaw", "character"),
  function(object, value) {
     slot(object,"laNotes")<- value
     object
  })
  
####################################### 
###########################################################################
## lanePair
###############

# Accessor methods for lanePair class

if(!isGeneric("expVals1"))
   setGeneric("expVals1", function(object) standardGeneric("expVals1"))
setMethod("expVals1", "lanePair", function(object) slot(object, "expVals1"))

if(!isGeneric("expVals2"))
   setGeneric("expVals2", function(object) standardGeneric("expVals2"))
setMethod("expVals2", "lanePair", function(object) slot(object, "expVals2"))

if(!isGeneric("PairGenes"))
   setGeneric("PairGenes", function(object) standardGeneric("PairGenes"))
setMethod("PairGenes", "lanePair", function(object) slot(object, "PairGenes"))

if(!isGeneric("PairNotes"))
   setGeneric("PairNotes", function(object) standardGeneric("PairNotes"))
setMethod("PairNotes", "lanePair", function(object) slot(object, "PairNotes"))

############
if(!isGeneric("Npairs"))
   setGeneric("Npairs", function(object) standardGeneric("Npairs"))
setMethod("Npairs", "lanePair", function(object) ncol(expVals1(object)))


if(!isGeneric("LogVal1"))
   setGeneric("LogVal1", function(object) standardGeneric("LogVal1"))
setMethod("LogVal1", "lanePair",
   function(object){
     r<-expVals1(object)
     if(length(expVals1(object))!=0){
       r<-log(ifelse(r>0, r, NA),2)
     }
     r
   }
)

if(!isGeneric("LogVal2"))
   setGeneric("LogVal2", function(object) standardGeneric("LogVal2"))
setMethod("LogVal2", "lanePair",
   function(object){
     r<-expVals2(object)
     if(length(expVals2(object))!=0){
       r<-log(ifelse(r>0, r, NA),2)
     }
     r
   }
)

if(!isGeneric("MVal"))
   setGeneric("MVal", function(object) standardGeneric("MVal"))
   setMethod("MVal", "lanePair",
   function(object){
     M<-matrix(nr=0,nc=0)
     if((length(expVals1(object))!=0) & (length(expVals2(object))!=0))
         M<-LogVal1(object)-LogVal2(object)
     M
   }
)

if(!isGeneric("AVal"))
   setGeneric("AVal", function(object) standardGeneric("AVal"))
   setMethod("AVal", "lanePair",
   function(object){
     A<-matrix(nr=0,nc=0)
     if((length(expVals1(object))!=0) & (length(expVals2(object))!=0))
        A<-(LogVal1(object)+LogVal2(object))/2
     A
   }
)
#########################
# Assignment methods for lanePair class

if( !isGeneric("expVals1<-") )
     setGeneric("expVals1<-", function(object, value)
               standardGeneric("expVals1<-"))

setReplaceMethod("expVals1", signature("lanePair", "matrix"),
  function(object, value) {
     slot(object,"expVals1")<- value
     object
  })

if( !isGeneric("expVals2<-") )
      setGeneric("expVals2<-", function(object, value)
               standardGeneric("expVals2<-"))

setReplaceMethod("expVals2", signature("lanePair", "matrix"),
  function(object, value) {
     slot(object,"expVals2")<- value
     object
  })

if( !isGeneric("PairGenes<-") )
      setGeneric("PairGenes<-", function(object, value)
               standardGeneric("PairGenes<-"))

setReplaceMethod("PairGenes", signature("lanePair", "matrix"),
  function(object, value) {
     slot(object,"PairGenes")<- value
     object
  })
  
if( !isGeneric("PairNotes<-") )
      setGeneric("PairNotes<-", function(object, value)
               standardGeneric("PairNotes<-"))

setReplaceMethod("PairNotes", signature("lanePair", "character"),
  function(object, value) {
     slot(object,"PairNotes")<- value
     object
  })
  
#################################################

###########################################################################
## PairNorm
##################
# Accessor methods for PairNorm class

if(!isGeneric("AVal"))
   setGeneric("AVal", function(object) standardGeneric("AVal"))
setMethod("AVal", "PairNorm", function(object) slot(object, "AVal"))

if(!isGeneric("MVal"))
   setGeneric("MVal", function(object) standardGeneric("MVal"))
setMethod("MVal", "PairNorm", function(object) slot(object, "MVal"))

if(!isGeneric("PairGenes"))
   setGeneric("PairGenes", function(object) standardGeneric("PairGenes"))
setMethod("PairGenes", "PairNorm", function(object) slot(object, "PairGenes"))

if(!isGeneric("PairNotes"))
   setGeneric("PairNotes", function(object) standardGeneric("PairNotes"))
setMethod("PairNotes", "PairNorm", function(object) slot(object, "PairNotes"))

if(!isGeneric("pMloc"))
   setGeneric("pMloc", function(object) standardGeneric("pMloc"))
setMethod("pMloc", "PairNorm", function(object) slot(object, "pMloc"))

if(!isGeneric("pMscale"))
   setGeneric("pMscale", function(object) standardGeneric("pMscale"))
setMethod("pMscale", "PairNorm", function(object) slot(object, "pMscale"))


if(!isGeneric("Npairs"))
   setGeneric("Npairs", function(object) standardGeneric("Npairs"))
setMethod("Npairs", "PairNorm", function(object) ncol(MVal(object)))


###############
# Methods for quantities that are not slots of PairNorm

if(!isGeneric("expVals1"))
   setGeneric("expVals1", function(object) standardGeneric("expVals1"))
setMethod("expVals1", "PairNorm",
   function(object){
     M <- MVal(object)
     A <- AVal(object)
     r <- M/2 + A
     2^r
   }
)

if(!isGeneric("expVals2"))
  setGeneric("expVals2", function(object) standardGeneric("expVals2"))
setMethod("expVals2", "PairNorm",
  function(object){
    M <- MVal(object)
    A <- AVal(object)
    g <- A - M/2
    2^g
  }
)

#########################

# Assignment methods for PairNorm class

if( !isGeneric("AVal<-") )
      setGeneric("AVal<-", function(object, value)
               standardGeneric("AVal<-"))

setReplaceMethod("AVal", signature("PairNorm", "matrix"),
  function(object, value) {
     slot(object,"AVal")<- value
     object
  })

if( !isGeneric("MVal<-") )
      setGeneric("MVal<-", function(object, value)
               standardGeneric("MVal<-"))

setReplaceMethod("MVal", signature("PairNorm", "matrix"),
  function(object, value) {
     slot(object,"MVal")<- value
     object
  })

if( !isGeneric("PairGenes<-") )
      setGeneric("PairGenes<-", function(object, value)
               standardGeneric("PairGenes<-"))

setReplaceMethod("PairGenes", signature("PairNorm", "matrix"),
  function(object, value) {
     slot(object,"PairGenes")<- value
     object
  })
  
if( !isGeneric("PairNotes<-") )
      setGeneric("PairNotes<-", function(object, value)
               standardGeneric("PairNotes<-"))

setReplaceMethod("PairNotes", signature("PairNorm", "character"),
  function(object, value) {
     slot(object,"PairNotes")<- value
     object
  })

if( !isGeneric("pMloc<-") )
      setGeneric("pMloc<-", function(object, value)
               standardGeneric("pMloc<-"))

setReplaceMethod("pMloc", signature("PairNorm", "matrix"),
  function(object, value) {
     slot(object,"pMloc")<- value
     object
  })
  
if( !isGeneric("pMscale<-") )
      setGeneric("pMscale<-", function(object, value)
               standardGeneric("pMscale<-"))

setReplaceMethod("pMscale", signature("PairNorm", "matrix"),
  function(object, value) {
     slot(object,"pMscale")<- value
     object
  })

############################################################

#######################################################
            
###########################################################################
## PairResult
######################

##########################################
###### Access methods for class PairResult
#################################
if(!isGeneric("Zscore"))
   setGeneric("Zscore", function(object) standardGeneric("Zscore"))
setMethod("Zscore", "PairResult", function(object) slot(object, "Zscore"))

if(!isGeneric("Pvalue"))
   setGeneric("Pvalue", function(object) standardGeneric("Pvalue"))
setMethod("Pvalue", "PairResult", function(object) slot(object, "Pvalue"))

if(!isGeneric("Qvalue"))
   setGeneric("Qvalue", function(object) standardGeneric("Qvalue"))
setMethod("Qvalue", "PairResult", function(object) slot(object, "Qvalue"))

if(!isGeneric("Qvalue2"))
   setGeneric("Qvalue2", function(object) standardGeneric("Qvalue2"))
setMethod("Qvalue2", "PairResult", function(object) slot(object, "Qvalue2"))

if(!isGeneric("DfGenes"))
   setGeneric("DfGenes", function(object) standardGeneric("DfGenes"))
setMethod("DfGenes", "PairResult", function(object) slot(object, "DfGenes"))

if(!isGeneric("DfGeneNames"))
   setGeneric("DfGeneNames", function(object) standardGeneric("DfGeneNames"))
setMethod("DfGeneNames", "PairResult",
   function(object){
     PairGenes(object)[DfGenes(object)]
   }
)

############################################################
###########Assignment methods for PairResult class
################################
if( !isGeneric("Zscore<-") )
      setGeneric("Zscore<-", function(object, value)
                 standardGeneric("Zscore<-"))

setReplaceMethod("Zscore", signature("PairResult", "matrix"),
  function(object, value) {
     slot(object,"Zscore")<- value
     object
  })

if( !isGeneric("Pvalue<-") )
      setGeneric("Pvalue<-", function(object, value)
                 standardGeneric("Pvalue<-"))

setReplaceMethod("Pvalue", signature("PairResult", "matrix"),
  function(object, value) {
     slot(object,"Pvalue")<- value
     object
  })
  
if( !isGeneric("Qvalue<-") )
      setGeneric("Qvalue<-", function(object, value)
                 standardGeneric("Qvalue<-"))

setReplaceMethod("Qvalue", signature("PairResult", "matrix"),
  function(object, value) {
     slot(object,"Qvalue")<- value
     object
  })

if( !isGeneric("Qvalue2<-") )
      setGeneric("Qvalue2<-", function(object, value)
                 standardGeneric("Qvalue2<-"))

setReplaceMethod("Qvalue2", signature("PairResult", "matrix"),
  function(object, value) {
     slot(object,"Qvalue2")<- value
     object
  })

if( !isGeneric("DfGenes<-") )
      setGeneric("DfGenes<-", function(object, value)
                 standardGeneric("DfGenes<-"))

setReplaceMethod("DfGenes", signature("PairResult", "matrix"),
  function(object, value) {
     slot(object,"DfGenes")<- value
     object
  })

######################################################################
######################
# Accessor methods for PairResult class

if(!isGeneric("AVal"))
   setGeneric("AVal", function(object) standardGeneric("AVal"))
setMethod("AVal", "PairResult", function(object) slot(object, "AVal"))

if(!isGeneric("MVal"))
   setGeneric("MVal", function(object) standardGeneric("MVal"))
setMethod("MVal", "PairResult", function(object) slot(object, "MVal"))

if(!isGeneric("PairGenes"))
   setGeneric("PairGenes", function(object) standardGeneric("PairGenes"))
setMethod("PairGenes", "PairResult", function(object) slot(object, "PairGenes"))

if(!isGeneric("PairNotes"))
   setGeneric("PairNotes", function(object) standardGeneric("PairNotes"))
setMethod("PairNotes", "PairResult", function(object) slot(object, "PairNotes"))

if(!isGeneric("pMloc"))
   setGeneric("pMloc", function(object) standardGeneric("pMloc"))
setMethod("pMloc", "PairResult", function(object) slot(object, "pMloc"))

if(!isGeneric("pMscale"))
   setGeneric("pMscale", function(object) standardGeneric("pMscale"))
setMethod("pMscale", "PairResult", function(object) slot(object, "pMscale"))


if(!isGeneric("Npairs"))
   setGeneric("Npairs", function(object) standardGeneric("Npairs"))
setMethod("Npairs", "PairResult", function(object) ncol(MVal(object)))

###############
# Methods for quantities that are not slots of PairResult

if(!isGeneric("expVals1"))
   setGeneric("expVals1", function(object) standardGeneric("expVals1"))
setMethod("expVals1", "PairResult",
   function(object){
     M <- MVal(object)
     A <- AVal(object)
     r <- M/2 + A
     r
   }
)

if(!isGeneric("expVals2"))
  setGeneric("expVals2", function(object) standardGeneric("expVals2"))
setMethod("expVals2", "PairResult",
  function(object){
    M <- MVal(object)
    A <- AVal(object)
    g <- A - M/2
    g
  }
)

#####################################

# Assignment methods for PairResult class

if( !isGeneric("AVal<-") )
      setGeneric("AVal<-", function(object, value)
                 standardGeneric("AVal<-"))

setReplaceMethod("AVal", signature("PairResult", "matrix"),
  function(object, value) {
     slot(object,"AVal")<- value
     object
  })

if( !isGeneric("MVal<-") )
      setGeneric("MVal<-", function(object, value)
                 standardGeneric("MVal<-"))

setReplaceMethod("MVal", signature("PairResult", "matrix"),
  function(object, value) {
     slot(object,"MVal")<- value
     object
  })

if( !isGeneric("PairGenes<-") )
      setGeneric("PairGenes<-", function(object, value)
                 standardGeneric("PairGenes<-"))

setReplaceMethod("PairGenes", signature("PairResult", "matrix"),
  function(object, value) {
     slot(object,"PairGenes")<- value
     object
  })
  
if( !isGeneric("PairNotes<-") )
      setGeneric("PairNotes<-", function(object, value)
                 standardGeneric("PairNotes<-"))

setReplaceMethod("PairNotes", signature("PairResult", "character"),
  function(object, value) {
     slot(object,"PairNotes")<- value
     object
  })

if( !isGeneric("pMloc<-") )
      setGeneric("pMloc<-", function(object, value)
                 standardGeneric("pMloc<-"))

setReplaceMethod("pMloc", signature("PairResult", "matrix"),
  function(object, value) {
     slot(object,"pMloc")<- value
     object
  })
  
if( !isGeneric("pMscale<-") )
      setGeneric("pMscale<-", function(object, value)
                 standardGeneric("pMscale<-"))

setReplaceMethod("pMscale", signature("PairResult", "matrix"),
  function(object, value) {
     slot(object,"pMscale")<- value
     object
  })
  
########################
#######################################################
