\name{getGeneExp}
\alias{getGeneExp}
\title{getGeneExp: count the number of reads and calculate the RPKM for each gene}
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
  \item{mapResultBatch}{a vector containing uniquely mapping result files for a sample.
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
  The method use \emph{the average length} of all the isoforms as a gene length when 
  calculating the RPKM values. This is different from Mortazavi,A. et al. (2008).
}
\references{
  Mortazavi,A. et al. (2008) Mapping and quantifying mammalian transcriptomes by RNA-seq. \emph{Nat. Methods}, \bold{5}, 621-628.
}
\seealso{
 \code{\link{DEGexp}},
 \code{\link{DEGseq}},
 \code{\link{readGeneExp}},
 \code{\link{kidneyChr21.bed}},
 \code{\link{liverChr21.bed}},
 \code{\link{refFlatChr21}}.
}
     
\examples{
  ## If the data files are collected in a zip archive, the following
  ## commands will first extract them to the temporary directory.
  
  kidneyR1L1 <- system.file("data", "kidneyChr21.bed.txt", package="DEGseq")
  refFlat    <- system.file("data", "refFlatChr21.txt", package="DEGseq")
  if((kidneyR1L1 == "")||(refFlat == "")){
     zipFile <- system.file("data", "Rdata.zip", package="DEGseq")
     if(zipFile != ""){
        unzip(zipFile, c("kidneyChr21.bed.txt", "refFlatChr21.txt"), exdir=tempdir())
        kidneyR1L1 <- file.path(tempdir(), "kidneyChr21.bed.txt")
        refFlat    <- file.path(tempdir(), "refFlatChr21.txt")
     }
  }
  mapResultBatch <- c(kidneyR1L1)
  output <- file.path(tempdir(), "kidneyChr21.bed.exp")
  getGeneExp(mapResultBatch, refFlat=refFlat, output=output)
  exp <- readGeneExp(file=output, geneCol=1, valCol=c(2,3), 
                     label=c("raw count", "RPKM"))
  exp[30:35,]
  cat("output: ", output, "\n")
}
\keyword{methods}
