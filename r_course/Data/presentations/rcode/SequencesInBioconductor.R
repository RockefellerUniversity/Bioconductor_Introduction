## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE

## ----load, echo=TRUE,eval=FALSE------------------------------------------
## library(BSgenome.Mmusculus.UCSC.mm10)
## class(BSgenome.Mmusculus.UCSC.mm10)

## ----load1, echo=FALSE,eval=TRUE-----------------------------------------
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
class(BSgenome.Mmusculus.UCSC.mm10)

## ----theObject, echo=TRUE,eval=TRUE,collapse=FALSE-----------------------
BSgenome.Mmusculus.UCSC.mm10

## ----contignames, echo=TRUE,eval=TRUE,collapse=FALSE---------------------
contigNames <- seqnames(BSgenome.Mmusculus.UCSC.mm10)
contigNames[1:22]

## ----contiglengths, echo=TRUE,eval=TRUE,collapse=FALSE-------------------
contigLengths <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
contigLengths[1:4]

## ----access, echo=TRUE,eval=TRUE,collapse=FALSE--------------------------
chr19_Seq <- BSgenome.Mmusculus.UCSC.mm10$chr19
chr19_Seq

## ----access2, echo=TRUE,eval=TRUE,collapse=FALSE-------------------------
chr19_Seq <- BSgenome.Mmusculus.UCSC.mm10[["chr19"]]
chr19_Seq

## ----dnastring, echo=TRUE,eval=TRUE,collapse=FALSE-----------------------
class(chr19_Seq)

## ----dnastringSub, echo=TRUE,eval=TRUE,collapse=FALSE--------------------
chr19_Seq[1:10000000]

## ----dnastringSub2, echo=TRUE,eval=TRUE,collapse=FALSE-------------------
subseq(chr19_Seq,start=1,end=10000000)

## ----dnastringSub3, echo=TRUE,eval=TRUE,collapse=FALSE-------------------
alphabetFrequency(chr19_Seq)

## ----dnastringSub4, echo=TRUE,eval=TRUE,collapse=FALSE-------------------
chr19_SeqComp <- complement(chr19_Seq)
chr19_SeqRev <- reverse(chr19_Seq)
chr19_SeqRevComp <- reverseComplement(chr19_Seq[10000000:10000010])
chr19_Seq[10000000:10000010]
chr19_SeqRevComp

## ----dnastringSub5, echo=TRUE,eval=TRUE,collapse=FALSE-------------------
chr19_Count <- countPattern(pattern="ATCTGCAATG",chr19_Seq)
chr19_Count

## ----dnastringSub611, echo=TRUE,eval=FALSE,collapse=FALSE----------------
## chr19_SeqSet <- DNAStringSet(chr19_Seq[10000000:10000010])
## names(chr19_SeqSet) <- "chr19"
## writeXStringSet(chr19_SeqSet,filepath = "Data/chr19_Seq.fa")
## 

## ----dnastringSub61s1, echo=FALSE,eval=TRUE,collapse=FALSE---------------
chr19_SeqSet <- DNAStringSet(chr19_Seq[10000000:10000010])
names(chr19_SeqSet) <- "chr19"
writeXStringSet(chr19_SeqSet,filepath = "../../Data/chr19_Seq.fa")


## ----dnastringSub62, echo=TRUE,eval=FALSE,tidy=TRUE----------------------
## chr19_FromFASTA <- readDNAStringSet(filepath = "Data/chr19_Seq.fa" )
## chr19_FromFASTA$chr19

## ----dnastringSubxzz62, echo=FALSE,eval=TRUE,tidy=TRUE-------------------
chr19_FromFASTA <- readDNAStringSet(filepath = "../../Data/chr19_Seq.fa" )
chr19_FromFASTA$chr19

