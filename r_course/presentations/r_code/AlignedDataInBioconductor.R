params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE


## ----setup2, include=FALSE,eval=FALSE,echo=FALSE------------------------------
## library(ShortRead)
## temp <- readFastq("~/Projects/Results/RNAseqPipeTest/FirstTest/FQs/ENCFF000CXH.fastq.gz")
## 
## ~/Downloads/out.fq
## temp <- readFastq("~/Downloads/out.fq")
## tAlin <- temp[sample(1:length(temp),10000)]
## writeFastq(tAlin,"~/Downloads/sampledActin.fq.gz")
## BiocInstaller::biocLite("QuasR")
## myFile <- data.frame(FileName="data/sampled_ENCFF000CXH.fastq.gz",SampleName="ENCFF000CXH")
## library("BSgenome.Hsapiens.UCSC.hg19")
## tis <- BSgenome.Hsapiens.UCSC.hg19[["chr5"]]
## writeXStringSet(DNAStringSet(list(chr5=tis)),"chr5.fa")
## write.table(myFile,"samples.txt",sep="\t",row.names=FALSE,quote=FALSE)
## qAlign("samples.txt","chr5.fa")
## library(Rsamtools)
## Rsamtools::sortBam("data/sampled_ENCFF000CXH_29a7bd074f7.bam","Sorted_sampled_ENCFF000CXH")
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


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("
class: inverse, center, middle

# Aligned Sequences

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("

# Aligned Sequences

---
"    
  )
  
}


## ----a1,echo=TRUE,eval=FALSE--------------------------------------------------
## library(Rsamtools)


## ----a2,echo=FALSE,eval=TRUE--------------------------------------------------
suppressPackageStartupMessages(library(Rsamtools))


## ----b1,echo=TRUE,eval=TRUE---------------------------------------------------
coordSorted <- sortBam("data/liver.bodyMap.bam",
                       "Sorted_liver")
coordSorted


## ----c1,echo=TRUE,eval=TRUE---------------------------------------------------
readnameSorted <- sortBam("data/liver.bodyMap.bam",
                          "SortedByName_liver",
                          byQname=TRUE)
readnameSorted


## ----d1,echo=TRUE,eval=TRUE---------------------------------------------------
coordSorted <- sortBam("data/liver.bodyMap.bam",
                          "Sorted_liver",
                          maxMemory=1)
coordSorted


## ----e1,echo=TRUE,eval=TRUE---------------------------------------------------
indexBam("Sorted_liver.bam")


## ----f1,echo=TRUE,eval=TRUE---------------------------------------------------
quickBamFlagSummary("Sorted_liver.bam")


## ----g1,echo=TRUE,eval=TRUE---------------------------------------------------
idxstatsBam("Sorted_liver.bam")


## ----aa1,echo=TRUE,eval=FALSE-------------------------------------------------
## BiocManager::install('GenomicAlignments')
## library(GenomicAlignments)


## ----aa2,echo=FALSE,eval=TRUE-------------------------------------------------
suppressPackageStartupMessages(library(GenomicAlignments))


## ----ba2,echo=TRUE,eval=TRUE--------------------------------------------------
myHeader <- scanBamHeader("Sorted_liver.bam")
str(myHeader)


## ----ca2,echo=TRUE,eval=TRUE--------------------------------------------------
names(myHeader)
names(myHeader$Sorted_liver.bam)


## ----da2,echo=TRUE,eval=TRUE--------------------------------------------------
myHeader$Sorted_liver.bam$targets


## ----ea2,echo=TRUE,eval=TRUE--------------------------------------------------
myHeader$Sorted_liver.bam$text


## ----fa2,echo=TRUE,eval=TRUE--------------------------------------------------
myHeader$Sorted_liver.bam$text["@HD"]


## ----ga2,echo=TRUE,eval=TRUE--------------------------------------------------
myHeader <- scanBamHeader("SortedByName_liver.bam")
myHeader$SortedByName_liver.bam$text["@HD"]


## ----ha2,echo=TRUE,eval=TRUE--------------------------------------------------
myHeader <- scanBamHeader("Sorted_liver.bam")
myHeader$Sorted_liver.bam$text["@PG"]


## ----ia2,echo=TRUE,eval=TRUE--------------------------------------------------
myReads <- readGAlignments("Sorted_liver.bam")
class(myReads)


## ----ja2,echo=TRUE,eval=TRUE--------------------------------------------------
myReads[1:2,]


## ----la2,echo=TRUE,eval=TRUE--------------------------------------------------
myReads[1:2,]


## ----ma2,echo=TRUE,eval=TRUE--------------------------------------------------
myReads[1:2,]


## ----na2,echo=TRUE,eval=TRUE--------------------------------------------------
seqnames(myReads)
start(myReads)[1:2]


## ----na222,echo=TRUE,eval=TRUE------------------------------------------------
cigar(myReads)[1:2]
njunc(myReads)[1:2]


## ----oa2,echo=TRUE,eval=TRUE--------------------------------------------------
myReads[strand(myReads) == "+"]


## ----pa2,echo=TRUE,eval=TRUE--------------------------------------------------
my5primeReads <- narrow(myReads, start=1, width = 1)
my5primeReads[1:2]


## ----qa2,echo=TRUE,eval=TRUE--------------------------------------------------
myReadsPos <- narrow(myReads[strand(myReads) == "+"],
                     start=1, width = 1)
myReadsNeg <- narrow(myReads[strand(myReads) == "-"],
                     end=-1, width = 1)

my5primeReads <- c(myReadsPos,myReadsNeg)
my5primeReads[1:2]


## ----ra2,echo=TRUE,eval=TRUE--------------------------------------------------
myReadAsGRanges <- granges(myReads,use.mcols = TRUE)
myReadAsGRanges


## ----sa2,echo=TRUE,eval=TRUE--------------------------------------------------
myReadAsGRangesList <- grglist(myReads,use.mcols = TRUE)
myReadAsGRangesList[njunc(myReads) == 1]


## ----ta2,echo=TRUE,eval=TRUE--------------------------------------------------
myReadAsGRanges <- granges(myReads, use.mcols = TRUE)
myReadsAgain <- as(myReadAsGRanges, "GAlignments")
myReadsAgain[1:2]


## ----ua2,echo=TRUE,eval=TRUE--------------------------------------------------
myReadAsGRanges <- granges(myReads, use.mcols = TRUE)
my5Prime <- resize(myReadAsGRanges, fix = "start", width = 1)
my5PrimeAsReads <- as(my5Prime, "GAlignments")
my5PrimeAsReads


## ----va2,echo=TRUE,eval=FALSE-------------------------------------------------
## library(rtracklayer)
## export(my5PrimeAsReads, con="myModifiedReads.bam")


## ----wa2,echo=TRUE,eval=TRUE--------------------------------------------------
myRanges <- GRanges("chr12", IRanges(98591400,98608400))
myParam <- ScanBamParam(which=myRanges)
myParam


## ----xa2,echo=TRUE,eval=TRUE--------------------------------------------------
filteredReads <- readGAlignments("Sorted_liver.bam", param = myParam)
filteredReads


## ----ya2,echo=TRUE,eval=TRUE--------------------------------------------------
myParam <- ScanBamParam(what=c("qname", "seq", "qual"))
infoInReads <- readGAlignments("Sorted_liver.bam", param = myParam)
infoInReads[1]


## ----za2,echo=TRUE,eval=TRUE--------------------------------------------------
mcols(infoInReads)


## ----aaa1,echo=TRUE,eval=TRUE-------------------------------------------------
bamHeader <- scanBamHeader("Sorted_liver.bam")
myChromosomes <- bamHeader$Sorted_liver.bam$targets
for(i in 1:length(myChromosomes)){
  grangesForImport <- GRanges(names(myChromosomes)[i],
                              IRanges(1,myChromosomes)[i])
  myParam <- ScanBamParam(which = grangesForImport)
  myReads <- readGAlignments("Sorted_liver.bam", 
                             param=myParam)
  print(length(myReads))
}

