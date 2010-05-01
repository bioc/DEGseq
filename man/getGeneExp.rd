\name{getGeneExp}
\alias{getGeneExp}
\title{getGeneExp: Count the number of reads and calculate the RPKM for each gene}
\description{
  This function is used to count the number of reads and calculate the RPKM for each gene.
  It takes uniquely mapped reads from RNA-seq data for a sample with a gene annotation file as input. 
  So users should map the reads (obtained from sequencing library of the sample) to the corresponding genome in advance.
}
\usage{
getGeneExp(mapResultBatch, fileFormat="bed", readLength=32, strandInfo=FALSE,
           refFlat, output=paste(mapResultBatch[1],".exp",sep=""), min.overlapPercent=1)
}
\arguments{
  \item{mapResultBatch}{vector containing uniquely mapping result files for a sample.
                        \cr \emph{Note}: The sample can have multiple technical replicates.}
  \item{fileFormat}{file format: \code{"bed"} or \code{"eland"}.
                    \cr example of \code{"bed"} format: \code{chr12    7    38    readID    2    +}
                    \cr example of \code{"eland"} format: \code{readID    chr12.fa    7    U2    F}
                    \cr \emph{Note}: The field separator character is \code{TAB}. And the files must
                        follow the format as one of the examples.
                   }
  \item{readLength}{the length of the reads (only used if \code{fileFormat="eland"}).}
  \item{strandInfo}{whether the strand information was retained during the cloning of the cDNAs.
                     \itemize{
                       \item \code{"TRUE" }: retained,
                       \item \code{"FALSE"}: not retained.
                     }
                   }
  \item{refFlat}{gene annotation file in UCSC refFlat format.
                 \cr See \url{http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat}.
                }
  \item{output}{the output file.}
  \item{min.overlapPercent}{the minimum percentage of the overlapping length for a read and an exon over 
                            the length of the read itself, for counting this read from the exon. should be <=1. 
                            \cr \code{0}: at least \code{1} bp overlap between a read and an exon.
                           }
}
\note{
  This function sums up the numbers of reads coming from all exons of a specific gene 
  (according to the known gene annotation) as the gene expression value. 
  The exons may include the 5'-UTR, protein coding region, and 3'-UTR of a gene. 
  All introns are ignored for a gene for the sequenced reads are from the spliced transcript library. 
  If a read falls in an exon (usually, a read is shorter than an exon), the read count for this exon plus 1. 
  If a read is crossing the boundary of an exon, users can tune the parameter \code{min.overlapPercent}, 
  which is the minimum percentage of the overlapping length for a read and an exon over the length of the read 
  itself, for counting this read from the exon. 
  The method use the union of all possible exons for calculating the length for each gene.
}
\references{
  Mortazavi,A. et al. (2008) Mapping and quantifying mammalian transcriptomes by RNA-seq. \emph{Nat. Methods}, \bold{5}, 621-628.
}
\seealso{
 \code{\link{DEGexp}},
 \code{\link{DEGseq}},
 \code{\link{getGeneExp.aln}},
 \code{\link{readGeneExp}},
 \code{\link{kidneyChr21.bed}},
 \code{\link{liverChr21.bed}},
 \code{\link{refFlatChr21}}.
}
     
\examples{
  kidneyR1L1 <- system.file("extdata", "kidneyChr21.bed.txt", package="DEGseq")
  refFlat    <- system.file("extdata", "refFlatChr21.txt", package="DEGseq")
  mapResultBatch <- list(kidneyR1L1)
  output <- file.path(tempdir(), "kidneyChr21.bed.exp")
  exp <- getGeneExp(mapResultBatch, refFlat=refFlat, output=output)
  write.table(exp[30:35,], row.names=FALSE)
  cat("output: ", output, "\n")
}
\keyword{methods}
