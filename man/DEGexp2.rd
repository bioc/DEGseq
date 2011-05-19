\name{DEGexp2}
\alias{DEGexp2}
\title{DEGexp2: Identifying Differentially Expressed Genes from gene expression data}
\description{
   This function is another (old) version of DEGexp. It takes the gene expression files as input instead of gene expression matrixs.  
}
\usage{
DEGexp2(geneExpFile1, geneCol1=1, expCol1=2, depth1=rep(0, length(expCol1)), groupLabel1="group1",
        geneExpFile2, geneCol2=1, expCol2=2, depth2=rep(0, length(expCol2)), groupLabel2="group2",
        header=TRUE, sep="", method=c("LRT", "CTR", "FET", "MARS", "MATR", "FC"), 
        pValue=1e-3, zScore=4, qValue=1e-3, foldChange=4, 
        thresholdKind=1, outputDir="none", normalMethod=c("none", "loess", "median"),
        replicate1="none", geneColR1=1, expColR1=2, depthR1=rep(0, length(expColR1)), replicateLabel1="replicate1",
        replicate2="none", geneColR2=1, expColR2=2, depthR2=rep(0, length(expColR2)), replicateLabel2="replicate2", rawCount=TRUE)
}
\arguments{
  \item{geneExpFile1}{file containing gene expression values for replicates of sample1 (or replicate1 when \code{method="CTR"}).}
  \item{geneCol1}{gene id column in geneExpFile1.}
  \item{expCol1}{expression value \emph{columns} in geneExpFile1 for replicates of sample1 (numeric vector).
                 \cr \emph{Note}: Each column corresponds to a replicate of sample1.
                }
  \item{depth1}{the total number of reads uniquely mapped to genome for each replicate of sample1 (numeric vector),
                \cr default: take the total number of reads mapped to all annotated genes as the depth for each replicate.}
  \item{groupLabel1}{label of group1 on the plots.}
  \item{geneExpFile2}{file containing gene expression values for replicates of sample2 (or replicate2 when \code{method="CTR"}).}
  \item{geneCol2}{gene id column in geneExpFile2.}
  \item{expCol2}{expression value \emph{columns} in geneExpFile2 for replicates of sample2 (numeric vector).
                 \cr \emph{Note}: Each column corresponds to a replicate of sample2.
                }
  \item{depth2}{the total number of reads uniquely mapped to genome for each replicate of sample2 (numeric vector),
                \cr default: take the total number of reads mapped to all annotated genes as the depth for each replicate.}
  \item{groupLabel2}{label of group2 on the plots.}
  \item{header}{a logical value indicating whether geneExpFile1 and geneExpFile2
                contain the names of the variables as its first line. See \code{?read.table}.}
  \item{sep}{the field separator character. If sep = "" (the default for read.table)
             the separator is \emph{white space}, that is one or more spaces, tabs, newlines or carriage returns.
             See \code{?read.table}.}
  \item{method}{method to identify differentially expressed genes. Possible methods are:
                 \itemize{
                  \item \code{ "LRT"}:  Likelihood Ratio Test (Marioni et al. 2008),
                  \item \code{ "CTR"}:  Check whether the variation between Technical Replicates
                                        can be explained by the random sampling model (Wang et al. 2009),
                  \item \code{ "FET"}:  Fisher's Exact Test (Joshua et al. 2009),
                  \item \code{"MARS"}:  MA-plot-based method with Random Sampling model (Wang et al. 2009),
                  \item \code{"MATR"}:  MA-plot-based method with Technical Replicates (Wang et al. 2009),
                  \item \code{ "FC" }:  Fold-Change threshold on MA-plot.
                 }
               }
  \item{pValue}{pValue threshold (for the methods: \code{LRT, FET, MARS, MATR}). 
                \cr only used when \code{thresholdKind=1}.}
  \item{zScore}{zScore threshold (for the methods: \code{MARS, MATR}). 
                \cr only used when \code{thresholdKind=2}.}
  \item{qValue}{qValue threshold (for the methods: \code{LRT, FET, MARS, MATR}).
                \cr only used when \code{thresholdKind=3} or \code{thresholdKind=4}.}
  \item{thresholdKind}{the kind of threshold. Possible kinds are:
                 \itemize{
                  \item \code{1}:  pValue threshold, 
                  \item \code{2}:  zScore threshold,
                  \item \code{3}:  qValue threshold (Benjamini et al. 1995),
                  \item \code{4}:  qValue threshold (Storey et al. 2003),
                  \item \code{5}:  qValue threshold (Storey et al. 2003) and Fold-Change threshold on MA-plot are both required (can be used only when \code{method="MARS"}).
                 }
               }
  \item{foldChange}{fold change threshold on MA-plot (for the method: \code{FC}).}
  \item{outputDir}{the output directory.}
  \item{normalMethod}{the normalization method: \code{"none", "loess", "median"} (Yang et al. 2002). 
                      \cr recommend: \code{"none"}. }
  \item{replicate1}{file containing gene expression values for replicate batch1 (only used when \code{method="MATR"}).
                    \cr \emph{Note}: replicate1 and replicate2 are two (groups of) technical replicates of a sample.}
  \item{geneColR1}{gene id column in the expression file for replicate batch1 (only used when \code{method="MATR"}).}
  \item{expColR1}{expression value \emph{columns} in the expression file for replicate batch1 (numeric vector) (only used when \code{method="MATR"}).}
  \item{depthR1}{the total number of reads uniquely mapped to genome for each replicate in replicate batch1 (numeric vector),
                 \cr default: take the total number of reads mapped to all annotated genes as the depth for each replicate (only used when \code{method="MATR"}).}
  \item{replicateLabel1}{label of replicate batch1 on the plots (only used when \code{method="MATR"}).}
  \item{replicate2}{file containing gene expression values for replicate batch2 (only used when \code{method="MATR"}).
                    \cr \emph{Note}: replicate1 and replicate2 are two (groups of) technical replicates of a sample.}
  \item{geneColR2}{gene id column in the expression file for replicate batch2 (only used when \code{method="MATR"}).}
  \item{expColR2}{expression value \emph{columns} in the expression file for replicate batch2 (numeric vector) (only used when \code{method="MATR"}).}
  \item{depthR2}{the total number of reads uniquely mapped to genome for each replicate in replicate batch2 (numeric vector),
                 \cr default: take the total number of reads mapped to all annotated genes as the depth for each replicate (only used when \code{method="MATR"}).}
  \item{replicateLabel2}{label of replicate batch2 on the plots (only used when \code{method="MATR"}).}
  \item{rawCount}{a logical value indicating the gene expression values are based on raw read counts or normalized values.}
}

\references{
  Benjamini,Y. and Hochberg,Y (1995). Controlling the false discovery rate: a practical and
  powerful approach to multiple testing. \emph{J. R. Stat. Soc. Ser. B} \bold{57}, 289-300. 
  
  Jiang,H. and Wong,W.H. (2008) Statistical inferences for isoform expression in RNA-seq.
  \emph{Bioinformatics}, \bold{25}, 1026-1032.

  Bloom,J.S. et al. (2009) Measuring differential gene expression by short read sequencing: quantitative comparison to
  2-channel gene expression microarrays. \emph{BMC Genomics},  \bold{10}, 221.

  Marioni,J.C. et al. (2008) RNA-seq: an assessment of technical reproducibility and comparison with gene expression arrays.
  \emph{Genome Res.}, \bold{18}, 1509-1517.
  
  Storey,J.D. and Tibshirani,R. (2003) Statistical significance for genomewide studies. \emph{Proc. Natl. Acad. Sci.} \bold{100}, 9440-9445.

  Wang,L.K. and et al. (2010) DEGseq: an R package for identifying differentially expressed genes from RNA-seq data, \emph{Bioinformatics} \bold{26}, 136 - 138.
    
  Yang,Y.H. et al. (2002) Normalization for cDNA microarray data: a robust composite method addressing single and multiple
  slide systematic variation. \emph{Nucleic Acids Research}, \bold{30}, e15.
}
\seealso{
 \code{\link{DEGexp}},
 \code{\link{DEGseq}},
 \code{\link{getGeneExp}},
 \code{\link{readGeneExp}},
 \code{\link{GeneExpExample1000}},
 \code{\link{GeneExpExample5000}}.
}
     
\examples{
  
  ## kidney: R1L1Kidney, R1L3Kidney, R1L7Kidney, R2L2Kidney, R2L6Kidney 
  ## liver: R1L2Liver, R1L4Liver, R1L6Liver, R1L8Liver, R2L3Liver
  
  geneExpFile <- system.file("extdata", "GeneExpExample5000.txt", package="DEGseq")
  outputDir <- file.path(tempdir(), "DEGexpExample")
  exp <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(7,9,12,15,18))
  exp[30:35,]
  exp <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(8,10,11,13,16))
  exp[30:35,]
  DEGexp2(geneExpFile1=geneExpFile, geneCol1=1, expCol1=c(7,9,12,15,18), groupLabel1="kidney",
          geneExpFile2=geneExpFile, geneCol2=1, expCol2=c(8,10,11,13,16), groupLabel2="liver",
          method="MARS", outputDir=outputDir)
  cat("outputDir:", outputDir, "\n")
}
\keyword{methods}
