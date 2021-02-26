params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE


## ----setup2, include=FALSE,eval=FALSE,echo=FALSE------------------------------
## library(ShortRead)
## temp <- readFastq("~/Projects/Results/RNAseqPipeTest/FirstTest/FQs/ENCFF000CXH.fastq.gz")
## fastqSample <- temp[1:100000]
## writeFastq(fastqSample,file = "~/Projects/Software/Github/RUBioconductor_Introduction/r_course/data/sampled_ENCFF000CXH.fastq.gz",mode = "w")


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("
class: inverse, center, middle

# Summarizing Scores

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("

# Summarizing Scores

---
"    
  )
  
}


## ----load, echo=TRUE,eval=FALSE-----------------------------------------------
## library(GenomicAlignments)


## ----load1, echo=FALSE,eval=TRUE,warning=FALSE--------------------------------
suppressPackageStartupMessages(library(GenomicAlignments))


## ----a1-----------------------------------------------------------------------
sortedHeart <- sortBam("data/heart.bodyMap.bam","Heart")
indexBam(sortedHeart)



## ----a1232,include=FALSE------------------------------------------------------
sortedHeart <- sortBam("data/liver.bodyMap.bam","Liver")
indexBam(sortedHeart)



## ----a1.2---------------------------------------------------------------------
heartCoverage <- coverage("Heart.bam")
class(heartCoverage)


## ----a2-----------------------------------------------------------------------
heartCoverage


## ----a3,fig.width=8,fig.height=2----------------------------------------------
chr12Cov <- heartCoverage[["chr12"]]
signalDepth <- chr12Cov[98591400:98608400]
signalDepthScaled <- data.frame(Position=98591400:98608400,
                                Signal=signalDepth*1000)
library(ggplot2)
ggplot(signalDepthScaled,aes(x=Position,y=Signal))+
  geom_line()+theme_minimal()



## ----a4-----------------------------------------------------------------------
heartAln <- readGAlignments("Heart.bam")
heartCov1 <- coverage(heartAln)


## ----frfga3,fig.width=8,fig.height=2,echo=FALSE-------------------------------
chr12Cov <- heartCov1[["chr12"]]
signalDepth <- chr12Cov[98591400:98608400]
signalDepthScaled <- data.frame(Position=98591400:98608400,
                                Signal=signalDepth*1000)
library(ggplot2)
ggplot(signalDepthScaled,aes(x=Position,y=Signal))+
  geom_line()+theme_minimal()


## ----a5-----------------------------------------------------------------------
heartGR <- granges(heartAln)
heartCov2 <- coverage(heartGR)
heartGRL <- grglist(heartAln)
heartCov3 <- coverage(heartGRL)


## ----frfa3,fig.width=8,fig.height=2,echo=FALSE--------------------------------
chr12Cov <- heartCov3[["chr12"]]
signalDepth <- chr12Cov[98591400:98608400]
signalDepthScaled <- data.frame(Position=98591400:98608400,
                                Signal=signalDepth*1000)
library(ggplot2)
ggplot(signalDepthScaled,aes(x=Position,y=Signal))+
  geom_line()+theme_minimal()


## ----a6s,echo=TRUE,eval=FALSE-------------------------------------------------
## heartAlnPos <- heartAln[strand(heartAln) == "+"]
## heartAlnPos <- coverage(heartAlnPos)
## heartAlnPos["chr12"]
## export.bw(heartAlnPos,con="heartPos.bw")


## ----a6,echo=FALSE,eval=TRUE--------------------------------------------------
heartAlnPos <- heartAln[strand(heartAln) == "+"]
heartAlnPos <- coverage(heartAlnPos)
heartAlnPos["chr12"]


## ----a7,collapse=TRUE---------------------------------------------------------
heartCoverageX10 <- coverage("Heart.bam",
                          weight = 10)
heartCoverageX10["chr12"]
heartCoverage["chr12"]


## ----a8-----------------------------------------------------------------------
allChromosomeStats <- idxstatsBam("Heart.bam")
allChromosomeStats


## ----a8k----------------------------------------------------------------------
mapped <- sum(allChromosomeStats[,"mapped"])
heartCoverageNorm <- coverage("Heart.bam",
                          weight = (10^6)/mapped)
heartCoverageNorm["chr12"]


## ----a91,echo=FALSE-----------------------------------------------------------
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))


## ----a9-----------------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
exonsOfGenes <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene,
                        by="gene")
slc25A3 <- exonsOfGenes[["5250"]]
slc25A3


## ----a9lk,echo=TRUE,eval=FALSE------------------------------------------------
## heartCoverageNorm <- coverage("Heart.bam")
## myMeanCovOverExons <- mean(heartCoverageNorm[slc25A3])
## myMeanCovOverExons


## ----a9lkD,echo=FALSE,eval=TRUE-----------------------------------------------
heartCoverageNorm <- coverage("Heart.bam")
myMeanCovOverExons <- BiocGenerics::mean(heartCoverageNorm[slc25A3])
myMeanCovOverExons


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("---
class: inverse, center, middle

# Summarizing counts in regions from alignments

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("---

# Summarizing counts in regions from alignments

---
"    
  )
  
}


## ----a11----------------------------------------------------------------------
geneBody <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
TSS <- promoters(geneBody,500,500)
myTSScounts <- summarizeOverlaps(TSS,"Heart.bam")
class(myTSScounts)


## ----a11a---------------------------------------------------------------------
myTSScounts


## ----a12,echo=TRUE,eval=FALSE-------------------------------------------------
## countMatrix <- assay(myTSScounts)
## countMatrix["5250",]


## ----a12s,echo=FALSE,eval=TRUE------------------------------------------------
countMatrix <- assay(myTSScounts)
countMatrix["5250",,drop=FALSE]


## ----a13,echo=TRUE------------------------------------------------------------
Granges <- rowRanges(myTSScounts)
Granges


## ----a131---------------------------------------------------------------------

geneExons <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene,by="gene")
geneExons["5250"]


## ----a14a---------------------------------------------------------------------
myGeneCounts <- summarizeOverlaps(geneExons,"Heart.bam")
myGeneCounts


## ----a122a,echo=TRUE,eval=FALSE-----------------------------------------------
## countMatrix <- assay(myGeneCounts)
## countMatrix["5250",]


## ----a122b,echo=FALSE,eval=TRUE-----------------------------------------------
countMatrix <- assay(myGeneCounts)
countMatrix["5250",,drop=FALSE]


## ----a1211a-------------------------------------------------------------------
grgList <- rowRanges(myGeneCounts)
grgList


## ----a1312a-------------------------------------------------------------------
allGeneCounts <- summarizeOverlaps(geneExons,
                                   c("Heart.bam","Liver.bam"))


## ----a122ssa,echo=TRUE,eval=FALSE---------------------------------------------
## countMatrix <- assay(allGeneCounts)
## countMatrix["5250",]


## ----a1dd22b,echo=FALSE,eval=TRUE---------------------------------------------
countMatrix <- assay(allGeneCounts)
countMatrix["5250",,drop=FALSE]


## ----a20----------------------------------------------------------------------
myBam <- BamFile("Heart.bam")
class(myBam)


## ----a21q---------------------------------------------------------------------
myBam <- BamFile("Heart.bam", yieldSize = 1000)
heartGeneCounts <- summarizeOverlaps(geneExons,myBam)
heartGeneCounts


## ----a22w---------------------------------------------------------------------
myBam <- BamFileList(c("Heart.bam","Liver.bam"),
                     yieldSize = 1000)
allGeneCounts <- summarizeOverlaps(geneExons,myBam)
allGeneCounts

