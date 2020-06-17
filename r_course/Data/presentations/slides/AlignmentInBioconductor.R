## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE

## ----setup2, include=FALSE,eval=FALSE,echo=FALSE-------------------------
## library(ShortRead)
## temp <- readFastq("~/Projects/Results/RNAseqPipeTest/FirstTest/FQs/ENCFF000CXH.fastq.gz")
## 
## ~/Downloads/out.fq
## temp <- readFastq("~/Downloads/out.fq")
## tAlin <- temp[sample(1:length(temp),10000)]
## writeFastq(tAlin,"~/Downloads/sampledActin.fq.gz")
## BiocInstaller::biocLite("QuasR")
## myFile <- data.frame(FileName="../../Data/sampled_ENCFF000CXH.fastq.gz",SampleName="ENCFF000CXH")
## library("BSgenome.Hsapiens.UCSC.hg19")
## tis <- BSgenome.Hsapiens.UCSC.hg19[["chr5"]]
## writeXStringSet(DNAStringSet(list(chr5=tis)),"chr5.fa")
## write.table(myFile,"samples.txt",sep="\t",row.names=FALSE,quote=FALSE)
## qAlign("samples.txt","chr5.fa")
## library(Rsamtools)
## Rsamtools::sortBam("../../Data/sampled_ENCFF000CXH_29a7bd074f7.bam","Sorted_sampled_ENCFF000CXH")
## Rsamtools::indexBam("Sorted_sampled_ENCFF000CXH.bam")
## myCoverage <- coverage("Sorted_sampled_ENCFF000CXH.bam")
## export.bw(myCoverage,con = "Sorted_sampled_ENCFF000CXH.bw")
## 
## Rsamtools::indexBam("~/Downloads/ENCFF846QSN.bam")
## 
## myFile <- data.frame(FileName="~/Downloads/sampledActin.fq.gz",SampleName="ENCFF000CXH")
## library("BSgenome.Hsapiens.UCSC.hg19")
## tis <- BSgenome.Hsapiens.UCSC.hg19[["chr7"]]
## writeXStringSet(DNAStringSet(list(chr7=tis)),"chr7.fa")
## write.table(myFile,"samples.txt",sep="\t",row.names=FALSE,quote=FALSE)
## qAlign("samples.txt","chr7.fa",splicedAlignment = TRUE)
## library(Rsamtools)
## Rsamtools::sortBam("~/Downloads/sampledActin_29a70b5f1d3.bam","sampledActinSpliced")
## Rsamtools::indexBam("sampledActinSpliced.bam")
## myCoverage <- coverage("Sorted_sampled_ENCFF000CXH.bam")
## export.bw(myCoverage,con = "Sorted_sampled_ENCFF000CXH.bw")
## 
## Rsamtools::indexBam("~/Downloads/ENCFF846QSN.bam")
## 

## ----load, echo=TRUE,eval=FALSE------------------------------------------
## source("https://bioconductor.org/biocLite.R")
## biocLite("QuasR")
## library(QuasR)

## ----load1q, echo=FALSE,eval=TRUE----------------------------------------
suppressPackageStartupMessages(library(QuasR))

## ----load11, echo=TRUE,eval=FALSE----------------------------------------
## source("https://bioconductor.org/biocLite.R")
## biocLite("BSgenome.Hsapiens.UCSC.hg38")
## library(BSgenome.Hsapiens.UCSC.hg38)

## ----load1, echo=FALSE,eval=TRUE-----------------------------------------
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))

## ----sampleTable1, echo=TRUE,eval=FALSE----------------------------------
## FileName <- "../../Data/sampledActin.fq.gz"
## SampleName <- "sampledActin"
## sampleTable <- data.frame(FileName,SampleName)
## write.table(sampleTable,file="sampleTable.txt",sep="\t",quote=FALSE,row.names = FALSE)

## ----sampleTable1s, echo=FALSE,eval=TRUE---------------------------------
FileName <- "../../Data/sampledActin.fq.gz"
SampleName <- "sampledActin" 
data.frame(FileName,SampleName)

## ----bsgenome, echo=TRUE,eval=FALSE--------------------------------------
## library(QuasR)
## qAlign("sampleTable.txt","BSgenome.Hsapiens.UCSC.hg38")

## ----bsgenomeffft, echo=TRUE,eval=FALSE----------------------------------
## library(BSgenome.Hsapiens.UCSC.hg38)
## chr7hg38 <- BSgenome.Hsapiens.UCSC.hg38[["chr7"]]
## chr7hg38Set <- DNAStringSet(list(chr7=chr7hg38))
## writeXStringSet(chr7hg38Set,file="chr7.fa")

## ----bsgenomek, echo=TRUE,eval=FALSE-------------------------------------
## qAlign("sampleTable.txt","chr7.fa")

## ----loadq, echo=TRUE,eval=FALSE-----------------------------------------
## source("https://bioconductor.org/biocLite.R")
## biocLite("Rsamtools")
## library(Rsamtools)

## ----loadq1, echo=FALSE,eval=TRUE----------------------------------------
suppressPackageStartupMessages(library(Rsamtools))

## ----sort, echo=TRUE,eval=TRUE-------------------------------------------
sortBam("../../Data/sampledActin_29a3a870bf.bam","SortedActB")

## ----index, echo=TRUE,eval=TRUE------------------------------------------
indexBam("SortedActB.bam")

## ----sum, echo=TRUE,eval=TRUE,collapse=TRUE------------------------------
quickBamFlagSummary("SortedActB.bam")

## ----bsgenomqek, echo=TRUE,eval=FALSE------------------------------------
## qAlign("sampleTable.txt","chr7.fa",splicedAlignment = TRUE)

## ----sortS, echo=TRUE,eval=TRUE------------------------------------------
sortBam("../../Data/sampledActin_29a42342e34.bam","SortedActBSpliced")
indexBam("SortedActBSpliced.bam")

## ----sum2, echo=TRUE,eval=TRUE,collapse=TRUE-----------------------------
quickBamFlagSummary("SortedActBSpliced.bam")

## ----sum2wcww, echo=TRUE,eval=FALSE,include=TRUE-------------------------
## library(Rsubread)

## ----sum2wcwssw, echo=TRUE,eval=FALSE,include=TRUE-----------------------
## buildindex("chr7","chr7.fa")

## ----sum2waacww, echo=TRUE,eval=FALSE,include=TRUE-----------------------
## subjunc("chr7","../../Data/sampledActin.fq.gz",
##         output_format = "BAM",
##         output_file = "../../Data/RsubreadsampledActin.bam")

## ----sum2wcaww, echo=TRUE,eval=FALSE,include=TRUE------------------------
## sortBam("../../Data/RsubreadsampledActin.bam",
##         "Sorted_RsubreadsampledActin")
## indexBam("Sorted_RsubreadsampledActin.bam")

## ----sum2wcwsssw, echo=TRUE,eval=TRUE,include=TRUE-----------------------
quickBamFlagSummary("Sorted_RsubreadsampledActin.bam")


## ----sum2ww, echo=FALSE,eval=FALSE,collapse=TRUE,include=FALSE-----------
## temp <- scanBam("SortedActBSpliced.bam")
## reads <- temp[[1]]$qname[is.na(temp[[1]]$mapq)]
## temp2 <- readGAlignments("~/Downloads/out.bam",param=ScanBamParam(what=c("qname", "strand", "pos", "qwidth")))
## temp3 <- temp2[temp2@elementMetadata$qname %in% reads]
## export(temp3,"missing.bam")
## 
## 
## library(Rsubread)
## Rsubread::buildindex("chr7","chr7.fa",gappedIndex = FALSE,indexSplit = FALSE)
## Rsubread::subjunc("/Users/tcarroll/Projects/Software/Github/RUBioconductor_Introduction/r_course/presentations/slides/chr7","/Users/tcarroll/Projects/Software/Github/RUBioconductor_Introduction/r_course/Data/sampledActin.fq.gz",output_format = "BAM",output_file = "/Users/tcarroll/Projects/Software/Github/RUBioconductor_Introduction/r_course/Data/RsubreadsampledActin.bam")
## 
## temp <- scanBam("../../Data/RsubreadsampledActin.bam")
## reads <- temp[[1]]$qname[is.na(temp[[1]]$mapq)]
## temp2 <- readGAlignments("~/Downloads/out.bam",param=ScanBamParam(what=c("qname", "strand", "pos", "qwidth")))
## temp3 <- temp2[temp2@elementMetadata$qname %in% reads]
## export(temp3,"missingRsubread.bam")
## 
## qAlign("sampleTable.txt","chr7.fa",splicedAlignment = TRUE,maxHits = 100)
## 
## 
## 

