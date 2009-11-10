####################################
#########  AllClasses.R
#########
#########      all classes in DEGseq
####################################
###########################################################################
## laneRaw

setClass("laneRaw",
   representation(
     expVals="matrix",
     laGenes="matrix",
     laNotes="character"
   )
)


###########################################################################
## lanePair

setClass("lanePair",
   representation(
      expVals1="matrix",
      expVals2="matrix",
      PairGenes="matrix",
      PairNotes="character"
   )
)


##########################################################################
## PairNorm

setClass("PairNorm",
   representation(
      AVal="matrix",
      MVal="matrix",
      PairGenes="matrix",
      PairNotes="character",
      pMloc="matrix",
      pMscale="matrix",
      PairNormCall="call"
   )
)

######################
## PairResult

setClass("PairResult",
   representation(
      AVal="matrix",
      MVal="matrix",
      PairGenes="matrix",
      PairNotes="character",
      pMloc="matrix",
      pMscale="matrix",
      PairNormCall="call",
      Zscore ="matrix",
      DfGenes="matrix",
      Pvalue="matrix",
      Qvalue="matrix",
      Qvalue2="matrix"
   )
)
#######################################################################
