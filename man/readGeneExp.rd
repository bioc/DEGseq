\name{readGeneExp}
\alias{readGeneExp}
\title{readGeneExp: read gene expression values to a matrix}
\description{
   This method is used to read gene expression values from a file to a matrix in R workspace. 
   So that the matrix can be used as input of other packages, such as \emph{edgeR}.
   The input of the method is a file that contains gene expression values. 
}
\usage{
   readGeneExp(file, geneCol=1, valCol=2, label = NULL, header=TRUE, sep="")
}
\arguments{
  \item{file}{file containing gene expression values.}
  \item{geneCol}{gene id column in file.}
  \item{valCol}{expression value \emph{columns} to be read in the file.}
  \item{label}{label for the columns.}
  \item{header}{a logical value indicating whether the file contains 
                the names of the variables as its first line. See \code{?read.table}.}
  \item{sep}{the field separator character. If sep = "" (the default for read.table)
             the separator is \emph{white space}, that is one or more spaces, tabs, newlines or carriage returns.
             See \code{?read.table}.}
}
\seealso{
 \code{\link{getGeneExp}},
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
  exp <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(7,9,12,15,18,8,10,11,13,16))
  exp[30:35,]
}
\keyword{methods}
