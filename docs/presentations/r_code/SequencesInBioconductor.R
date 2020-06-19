## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE


## ----load, echo=TRUE,eval=FALSE-----------------------------------------------
## library(BSgenome.Mmusculus.UCSC.mm10)
## class(BSgenome.Mmusculus.UCSC.mm10)


## ----load1, echo=FALSE,eval=TRUE----------------------------------------------
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
class(BSgenome.Mmusculus.UCSC.mm10)


## ----theObject, echo=TRUE,eval=TRUE,collapse=FALSE----------------------------
BSgenome.Mmusculus.UCSC.mm10


## ----contignames, echo=TRUE,eval=TRUE,collapse=FALSE--------------------------
contigNames <- seqnames(BSgenome.Mmusculus.UCSC.mm10)
contigNames[1:22]


## ----contiglengths, echo=TRUE,eval=TRUE,collapse=FALSE------------------------
contigLengths <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
contigLengths[1:4]


## ----access, echo=TRUE,eval=TRUE,collapse=FALSE-------------------------------
chr19_Seq <- BSgenome.Mmusculus.UCSC.mm10$chr19
chr19_Seq


## ----access2, echo=TRUE,eval=TRUE,collapse=FALSE------------------------------
chr19_Seq <- BSgenome.Mmusculus.UCSC.mm10[["chr19"]]
chr19_Seq


## ----DNAString, echo=TRUE,eval=TRUE,collapse=FALSE----------------------------
class(chr19_Seq)


## ----DNAStringSub, echo=TRUE,eval=TRUE,collapse=FALSE-------------------------
chr19_Seq[1:10000000]


## ----DNAStringSub2, echo=TRUE,eval=TRUE,collapse=FALSE------------------------
subseq(chr19_Seq,start=1,end=10000000)


## ----DNAStringSub3, echo=TRUE,eval=TRUE,collapse=FALSE------------------------
alphabetFrequency(chr19_Seq)


## ----DNAStringSub4, echo=TRUE,eval=TRUE,collapse=FALSE------------------------
chr19_SeqComp <- complement(chr19_Seq)
chr19_SeqRev <- reverse(chr19_Seq)
chr19_SeqRevComp <- reverseComplement(chr19_Seq[10000000:10000010])
chr19_Seq[10000000:10000010]
chr19_SeqRevComp


## ----DNAStringSub4a, echo=TRUE,eval=TRUE,collapse=FALSE-----------------------
length(chr19_Seq[10000000:10000008])
chr19_SeqTranslation <- translate(chr19_Seq[10000000:10000008])
chr19_SeqTranslation


## ----DNAStringSub5k, echo=TRUE,eval=TRUE,collapse=FALSE,tidy=FALSE------------
chr19_Count <- countPattern(pattern="ATCTGCAATG",
                            subject=chr19_Seq)
chr19_Count


## ----DNAStringSubdd5k, echo=TRUE,eval=TRUE,collapse=FALSE,tidy=FALSE----------
chr19_Count <- countPattern(pattern="ATCTGCAATG",
                            subject=chr19_Seq,
                            max.mismatch = 2,
                            min.mismatch = 0)
chr19_Count


## ----DNAStringSub5ssk, echo=TRUE,eval=TRUE,collapse=FALSE,tidy=FALSE----------
IUPAC_CODE_MAP


## ----DNAStringSub5ssssk, echo=TRUE,eval=TRUE,collapse=FALSE,tidy=FALSE--------
chr19_Count <- countPattern(pattern="RYKHBNKYSRR",
                            subject=chr19_Seq,
                            fixed=FALSE)
chr19_Count


## ----DNAStringSub5l, echo=TRUE,eval=TRUE,collapse=FALSE,tidy=FALSE------------
chr19_Count_Watson <- countPattern(pattern="ATCTGCAATG",
                                    subject=chr19_Seq)
chr19_Count_Crick <- countPattern(pattern="ATCTGCAATG",
                                    subject=reverseComplement(chr19_Seq)
                                   )
Total_chr19_Count <- chr19_Count_Watson+chr19_Count_Crick



## ----DNAStringSub611, echo=TRUE,eval=FALSE,collapse=FALSE---------------------
## chr19_SeqSet <- DNAStringSet(chr19_Seq[10000000:10000010])
## names(chr19_SeqSet) <- "chr19"
## writeXStringSet(chr19_SeqSet,filepath = "data/chr19_Seq.fa")
## 


## ----DNAStringSub61s1, echo=FALSE,eval=TRUE,collapse=FALSE--------------------
chr19_SeqSet <- DNAStringSet(chr19_Seq[10000000:10000010])
names(chr19_SeqSet) <- "chr19"
writeXStringSet(chr19_SeqSet,filepath = "data/chr19_Seq.fa")



## ----DNAStringSub62, echo=TRUE,eval=FALSE,tidy=TRUE---------------------------
## chr19_FromFASTA <- readDNAStringSet(filepath = "data/chr19_Seq.fa" )
## chr19_FromFASTA$chr19


## ----DNAStringSubxzz62, echo=FALSE,eval=TRUE,tidy=TRUE------------------------
chr19_FromFASTA <- readDNAStringSet(filepath = "data/chr19_Seq.fa" )
chr19_FromFASTA$chr19

