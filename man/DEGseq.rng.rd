\name{DEGseq.rng}
\alias{DEGseq.rng}
\title{DEGseq.rng: Identify Differentially Expressed Genes from RNA-seq data (deal with the objects of class RangedData)}
\description{
  This function is used to identify differentially expressed genes from RNA-seq data. It takes uniquely mapped
  reads from RNA-seq data for the two samples with a gene annotation as input. 
  So users should map the reads (obtained from sequencing libraries of the samples) to the corresponding genome in advance.
}
\usage{
DEGseq.rng(rngBatch1, rngBatch2, 
           strandInfo=FALSE, refFlat, groupLabel1="group1", groupLabel2="group2",
           method=c("LRT", "CTR", "FET", "MARS", "MATR", "FC"), 
           pValue=1e-3, zScore=4, qValue=1e-3, foldChange=4, thresholdKind=1,
           outputDir="none", normalMethod=c("none", "loess", "median"),
           depthKind=1, replicateRngBatch1=NULL, replicateRngBatch2=NULL,
           replicateLabel1="replicate1", replicateLabel2="replicate2")
}
\arguments{
  \item{rngBatch1}{list containing uniquely mapping reads (objects of class RangedData) for technical replicates of sample1 (or replicate1 when \code{method="CTR"}).}
  \item{rngBatch2}{list containing uniquely mapping reads (objects of class RangedData) for technical replicates of sample2 (or replicate2 when \code{method="CTR"}).}
  \item{strandInfo}{whether the strand information was retained during the cloning of the cDNAs.
                    \itemize{
                    \item \code{"TRUE" }: retained,
                    \item \code{"FALSE"}: not retained.
                    }
                   }
  \item{refFlat}{gene annotation file in UCSC refFlat format. 
                 \cr See \url{http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat}. 
                } 
  \item{groupLabel1}{label of group1 on the plots.}
  \item{groupLabel2}{label of group2 on the plots.}
  \item{method}{method to identify differentially expressed genes. Possible methods are:
                 \itemize{
                  \item \code{ "LRT"}:  Likelihood Ratio Test (Marioni et al. 2008),
                  \item \code{ "CTR"}:  Check whether the variation between two Technical Replicates
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
                  \item \code{4}:  qValue threshold (Storey et al. 2003).
                 }
               }
  \item{foldChange}{fold change threshold on MA-plot (for the method: \code{FC}).}
  \item{outputDir}{the output directory.}
  \item{normalMethod}{the normalization method: \code{"none", "loess", "median"} (Yang,Y.H. et al. 2002). \cr
                      recommend: \code{"none"}. }
  \item{depthKind}{\code{1}: take the total number of reads uniquely mapped to genome as the depth for each replicate, \cr
                   \code{0}: take the total number of reads uniquely mapped to all annotated genes as the depth for each replicate. \cr
                   We recommend taking \code{depthKind=1}, 
                   especially when the genes in annotation file are part of all genes.} 
  \item{replicateRngBatch1}{list containing uniquely mapping reads (objects of class RangedData) obtained from replicate batch1 (only used when \code{method="MATR"}).}
  \item{replicateRngBatch2}{list containing uniquely mapping reads (objects of class RangedData) obtained from replicate batch2 (only used when \code{method="MATR"}).}
  \item{replicateLabel1}{label of replicate batch1 on the plots (only used when \code{method="MATR"}).}
  \item{replicateLabel2}{label of replicate batch2 on the plots (only used when \code{method="MATR"}).}
}
\note{
  Users should use \code{\link{DEGseq}} instead of this function when the short reads are too many to be loaded into memory.
}
\references{
  Benjamini,Y. and Hochberg,Y. (1995) Controlling the false discovery rate: a practical and
  powerful approach to multiple testing. \emph{J. R. Stat. Soc. Ser. B} \bold{57}, 289-300.

  Jiang,H. and Wong,W.H. (2009) Statistical inferences for isoform expression in RNA-seq.
  \emph{Bioinformatics}, \bold{25}, 1026-1032.

  Bloom,J.S. et al. (2009) Measuring differential gene expression by short read sequencing: quantitative comparison to
  2-channel gene expression microarrays. \emph{BMC Genomics},  \bold{10}, 221.

  Marioni,J.C. et al. (2008) RNA-seq: an assessment of technical reproducibility and comparison with gene expression arrays.
  \emph{Genome Res.}, \bold{18}, 1509-1517.

  Storey,J.D. and Tibshirani,R. (2003) Statistical significance for genomewide studies. \emph{Proc. Natl. Acad. Sci.} \bold{100}, 9440-9445.
  
  Wang,L.K. and et al. (2009) DEGseq: an R package to identify differentially expressed genes from RNA-seq data. Submitted.
  
  Yang,Y.H. et al. (2002) Normalization for cDNA microarray data: a robust composite method addressing single and multiple
  slide systematic variation. \emph{Nucleic Acids Research}, \bold{30}, e15.
}
\seealso{
 \code{\link{DEGexp}},
 \code{\link{DEGseq}},
 \code{\link{DEGseq.aln}},
 \code{\link{getGeneExp.rng}},
 \code{\link{readGeneExp}},
 \code{\link{kidneyChr21.bed}},
 \code{\link{liverChr21.bed}},
 \code{\link{refFlatChr21}}.
}
\examples{
  kidneyR1L1 <- system.file("extdata", "kidneyChr21.bed.txt", package="DEGseq")
  liverR1L2  <- system.file("extdata", "liverChr21.bed.txt", package="DEGseq")
  refFlat    <- system.file("extdata", "refFlatChr21.txt", package="DEGseq")
  kidneyR1L1_track <- rtracklayer::import(kidneyR1L1, format="bed")
  liverR1L2_track <- rtracklayer::import(liverR1L2, format="bed")
  rngBatch1 <- list(kidneyR1L1_track)  ## only use the data from kidneyR1L1 and liverR1L2
  rngBatch2 <- list(liverR1L2_track)
  outputDir <- file.path(tempdir(), "DEGseqRngExample")
  DEGseq.rng(rngBatch1, rngBatch2, refFlat=refFlat, outputDir=outputDir, method="LRT")
  cat("outputDir:", outputDir, "\n")
}
\keyword{methods}
