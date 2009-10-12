\name{samWrapper}
\alias{samWrapper}
\title{samWrapper: A Wrapper (with some modification) of the functions in the package samr 
       to identify differentially expressed genes for the RNA-seq data from two
       groups of paired or unpaired samples.}
\description{
   This function is a wrapper of the functions in \emph{samr}. It is used to identify differentially expressed genes 
   for two sets of samples with multiple replicates or two groups of samples from different individuals (e.g. disease samples vs. control samples).
   For the advanced users, please see \emph{samr} \url{http://cran.r-project.org/web/packages/samr/index.html}
   for detail.
}
\usage{
samWrapper(geneExpFile1, geneCol1=1, expCol1=2, measure1=rep(1, length(expCol1)),
           geneExpFile2, geneCol2=1, expCol2=2, measure2=rep(2, length(expCol2)),
           header=TRUE, sep="", paired=FALSE, s0=NULL, s0.perc=NULL, nperms=100, 
           testStatistic=c("standard","wilcoxon"), max.qValue=1e-3, min.foldchange=0,
           logged2=FALSE, output)
}
\arguments{
  \item{geneExpFile1}{file containing gene expression values for group1.}
  \item{geneCol1}{gene id column in geneExpFile1.}
  \item{expCol1}{expression value \emph{columns} in geneExpFile1. See the example.}
  \item{measure1}{numeric vector of outcome measurements for group1.
                  \cr like c(1,1,1...) when \code{paired=FALSE},
                  \cr or like c(-1,-2,-3,...) when \code{paired=TRUE}.
                 }
  \item{geneExpFile2}{file containing gene expression values for group2.}
  \item{geneCol2}{gene id column in geneExpFile2.}
  \item{expCol2}{expression value \emph{columns} in geneExpFile2. See the example.}
  \item{measure2}{numeric vector of outcome measurements for group2.
                  \cr like c(2,2,2...) when \code{paired=FALSE},
                  \cr or like c(1,2,3,...) when \code{paired=TRUE}.
                 }
  \item{header}{a logical value indicating whether geneExpFile1 and geneExpFile2
                contain the names of the variables as its first line. See \code{?read.table}.}
  \item{sep}{the field separator character. If sep = "" (the default for read.table)
             the separator is \emph{white space}, that is one or more spaces, tabs, newlines or carriage returns.
             See \code{?read.table}.}
  \item{paired}{a logical value indicating whether the samples are paired.}
  \item{s0}{exchangeability factor for denominator of test statistic; Default is automatic choice.}
  \item{s0.perc}{percentile of standard deviation values to use for s0; default is automatic choice; 
                 \cr -1 means s0=0 (different from s0.perc=0, meaning s0=zeroeth percentile of standard 
                 \cr deviation values= min of sd values.}
  \item{nperms}{number of permutations used to estimate false discovery rates.}
  \item{testStatistic}{test statistic to use in two class unpaired case. Either \code{"standard"} (t-statistic) 
                       \cr or \code{"wilcoxon"} (Two-sample wilcoxon or Mann-Whitney test).
                       \cr recommend \code{"standard"}.}
  \item{max.qValue}{the max qValue desired; shoube be <1.}
  \item{min.foldchange}{the minimum fold change desired; should be >1.
                       \cr default is zero, meaning no fold change criterion is applied.}
  \item{logged2}{a logical value indicating whether the expression values are logged2.}
  \item{output}{the output file.}
}

\references{
  Tusher,V., and et al. (2001): Significance analysis of microarrays applied to the ionizing radiation response \emph{PNAS} \bold{98}, 5116-5121.
  
  Tibshirani,R, and et al.: samr \url{http://cran.r-project.org/web/packages/samr/index.html}.
  
  A more complete description is given in the SAM manual at \url{http://www-stat.stanford.edu/~tibs/SAM}.
}
\seealso{
  \code{\link{DEGexp}},
  \code{\link{DEGseq}},
  \code{\link{GeneExpExample1000}},
  \code{\link{GeneExpExample5000}}.
}
     
\examples{
  ## If the data files are collected in a zip archive, the following
  ## commands will first extract them to the temporary directory.
  
  geneExpFile <- system.file("data", "GeneExpExample1000.txt", package="DEGseq")
  if(geneExpFile == ""){
     zipFile <- system.file("data", "Rdata.zip", package="DEGseq")
     if(zipFile != ""){
        unzip(zipFile, "GeneExpExample1000.txt", exdir=tempdir())
        geneExpFile <- file.path(tempdir(), "GeneExpExample1000.txt")
     }
  }
  set.seed(100)
  geneExpFile1 <- geneExpFile 
  geneExpFile2 <- geneExpFile
  output <- file.path(tempdir(), "samWrapperOut.txt")
  exp <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(7,9,12,15,18))
  exp[30:35,]
  exp <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(8,10,11,13,16))
  exp[30:35,]
  samWrapper(geneExpFile1=geneExpFile1, geneCol1=1, expCol1=c(7,9,12,15,18), measure1=c(-1,-2,-3,-4,-5),
             geneExpFile2=geneExpFile2, geneCol2=1, expCol2=c(8,10,11,13,16), measure2=c(1,2,3,4,5),
             nperms=100, min.foldchange=2, max.qValue=1e-4, output=output, paired=TRUE)
  cat("output:", output, "\n")
}
\keyword{methods}
