## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE


## ----setssssup, eval=FALSE----------------------------------------------------
## BiocManager::install(version = "3.10")
## BiocManager::install("GenomicRanges")


## ----sedfs,eval=TRUE,echo=FALSE-----------------------------------------------
suppressPackageStartupMessages(library(GenomicRanges))


## ----ss,eval=FALSE------------------------------------------------------------
## library(GenomicRanges)


## ----sttsss-------------------------------------------------------------------
myIntervals <- IRanges(start=c(1,10,20),end=c(2,11,30))
class(myIntervals)


## ----sttxxxsx-----------------------------------------------------------------
myIntervals


## ----sttxxsjjx----------------------------------------------------------------
myGenomicIntervals <- GRanges(seqnames=c("chr1","chr1","chr2"),
                              myIntervals)
class(myGenomicIntervals)


## ----stxxtsdx-----------------------------------------------------------------
myGenomicIntervals


## ----sttsxxx------------------------------------------------------------------
myGenomicIntervals[c(2,3),]


## ----stxxtsx------------------------------------------------------------------
myGenomicIntervals[c(2,3)]


## ----stddxxtsxss--------------------------------------------------------------
start(myGenomicIntervals)


## ----stddxxtssss--------------------------------------------------------------
end(myGenomicIntervals)


## ----stddxxtsxssss------------------------------------------------------------
start(myGenomicIntervals) <- c(1,5,15)


## ----stddxxtsxssds------------------------------------------------------------
start(myGenomicIntervals)


## ----stddxxtsxssssddss,error=TRUE---------------------------------------------
start(myGenomicIntervals) <- c(10,50,150)


## ----stddxxtsxdddsssddss,error=TRUE-------------------------------------------
contigNames <- seqnames(myGenomicIntervals)
contigNames


## ----stddaaxxtsxssddss,error=TRUE---------------------------------------------
as.character(contigNames)


## ----stddssxxtsxssddss,error=TRUE---------------------------------------------
seqnames(myGenomicIntervals) <- c("chr2","chr2","chr1")
seqnames(myGenomicIntervals)


## ----stddxxaadtsxassddddss,error=TRUE-----------------------------------------
seqnames(myGenomicIntervals) <- c("chr3","chr2","chr1")


## ----stddxxddstsxsddsddss,error=TRUE------------------------------------------
seqlevels(myGenomicIntervals)


## ----stddxxssdtsxsddsddss,error=TRUE------------------------------------------
seqlevels(myGenomicIntervals) <- c("chr1","chr2","chr3")
seqlevels(myGenomicIntervals)


## ----stddxxddwwstsxsddsddss,error=TRUE----------------------------------------
seqnames(myGenomicIntervals) <- c("chr1","chr2","chr3")
myGenomicIntervals


## ----aaeroe,error=TRUE--------------------------------------------------------

myGenomicIntervals <- GRanges(seqnames=c("chr1","chr1","chr2"), 
                              myIntervals,
                              strand=c("+","-","*"))
myGenomicIntervals


## ----atra,error=TRUE----------------------------------------------------------

strand(myGenomicIntervals) <- c("+","+","-")
strand(myGenomicIntervals)


## ----aera,error=TRUE----------------------------------------------------------

myGenomicIntervals <- GRanges(seqnames=c("chr1","chr1","chr2"),
                              myIntervals,
                              strand=c("+","+","+"))
strand(myGenomicIntervals)


## ----aaaq,error=TRUE----------------------------------------------------------
names(myGenomicIntervals) <- c("peak1","peak2","peak3")
myGenomicIntervals


## ----afffa,error=TRUE---------------------------------------------------------
myGenomicIntervals["peak2"]


## ----asa,error=TRUE-----------------------------------------------------------
myIntervals <- IRanges(start=c(1,2,2),end=c(2,11,30))
myGenomicIntervals <- GRanges(seqnames=c("chr1","chr1","chr2"),
                              myIntervals,strand=c("+","+","+"),
                              Score=c(10,20,40),
                              Comment=c("GoodQC","GoodQC","BadQC"))



## ----aasaaaa,error=TRUE-------------------------------------------------------
myGenomicIntervals$Score
myGenomicIntervals$Comment


## ----aasaaa,error=TRUE--------------------------------------------------------
mcols(myGenomicIntervals)


## ----aaszzzzzaaaaa,error=TRUE-------------------------------------------------
as.data.frame(myGenomicIntervals)


## ----aaszzzzzaaadsdaa,error=TRUE----------------------------------------------
GRanges("chr1:110-120")


## ----aaszzzzzaxaaaaaca,error=TRUE---------------------------------------------
library(stringi)
myRange <- "chr1:110:120"
newRange <- stri_replace_last_fixed(myRange, ':', '-')
newRange


## ----aaszzzzzaaassaca,error=TRUE----------------------------------------------
library(stringi)
myRange <- "MyID:chr1:110:120"
newRange <- stri_replace_last_fixed(myRange, ':', '-')
newerRange <- stri_replace_first_regex(newRange, '\\w*:', '')
newerRange
GRanges(newerRange)


## ----aaaskjbvasa,error=TRUE---------------------------------------------------
myGenomicIntervals[1]
shift(myGenomicIntervals[1],shift = 10)


## ----aaasatcsa,error=TRUE-----------------------------------------------------
myGenomicIntervals[3]
resize(myGenomicIntervals[3], width=20)


## ----aaasdasalpi,error=TRUE---------------------------------------------------
myGenomicIntervals[3]
resize(myGenomicIntervals[3], width=20, fix="end")


## ----aaasakjsfvrfa,error=TRUE-------------------------------------------------
resize(myGenomicIntervals[3],width=20, fix="start")
strand(myGenomicIntervals)[3] <- "-"
resize(myGenomicIntervals[3],width=20, fix="start")


## ----ajfddkkmsa,error=TRUE----------------------------------------------------
myGenomicIntervals <- GRanges(seqnames=c("chr1","chr1","chr2"),
                              ranges=IRanges(start=c(1,2,2),end=c(2,11,30)),
                              strand=c("+","+","+"))
myGenomicIntervals


## ----ajfddslka,error=TRUE-----------------------------------------------------
mergedGenomicIntervals <- reduce(myGenomicIntervals)
mergedGenomicIntervals


## ----ajfddujlsa,error=TRUE----------------------------------------------------
strand(myGenomicIntervals) <- c("+","-","+")
mergedGenomicIntervals <- reduce(myGenomicIntervals)
mergedGenomicIntervals


## ----ajfddusxjihlsa,error=TRUE------------------------------------------------
strand(myGenomicIntervals) <- c("+","-","+")
mergedGenomicIntervals <- reduce(myGenomicIntervals,
                                 ignore.strand=TRUE)
mergedGenomicIntervals


## ----ajfddfucfdjlsa,error=TRUE------------------------------------------------
myGenomicIntervals1 <- GRanges(seqnames=c("chr1","chr1"),
                              ranges=IRanges(start=c(1,25),
                                             end=c(20,30)),
                              strand=c("+","+"))

myGenomicIntervals2 <- GRanges(seqnames=c("chr1","chr1"),
                              ranges=IRanges(start=c(22,100),
                                             end=c(27,130)),
                              strand=c("+","+"))


## ----ajfddfuhjlsa,error=TRUE--------------------------------------------------
myGenomicIntervals1
myGenomicIntervals2


## ----ajfddfuhjcfvlsa,error=TRUE-----------------------------------------------
myGenomicIntervals1 %over% myGenomicIntervals2


## ----ajfdcfdfuhjlsa,error=TRUE------------------------------------------------
myGenomicIntervals1[myGenomicIntervals1 %over% myGenomicIntervals2]


## ----dccdc,error=TRUE---------------------------------------------------------
myOverlaps <- findOverlaps(myGenomicIntervals1,myGenomicIntervals2)
class(myOverlaps)
myOverlaps 


## ----cdkkcdcdc,error=TRUE-----------------------------------------------------
queryHits(myOverlaps)
subjectHits(myOverlaps)


## ----rekkpkvre,error=TRUE-----------------------------------------------------
myGenomicIntervals1[queryHits(myOverlaps)]
myGenomicIntervals2[subjectHits(myOverlaps)]


## ----ajfddfusscdcscfdjlsa,error=TRUE------------------------------------------
myGenomicIntervals1 <- GRanges(seqnames=c("chr1","chr1"),
                              ranges=IRanges(start=c(10,20),
                                             end=c(25,30)),
                              strand=c("+","+"))

myGenomicIntervals2 <- GRanges(seqnames=c("chr1","chr1"),
                              ranges=IRanges(start=c(1,10000),
                                             end=c(2,10002)),
                              strand=c("+","+"))



## ----ajfddfussscfdxcscsjlsa,error=TRUE----------------------------------------
myGenomicIntervals1
myGenomicIntervals2


## ----ajfddfussseecrfdjlsa,error=TRUE------------------------------------------
indexOfNearest <- nearest(myGenomicIntervals1,myGenomicIntervals2)
indexOfNearest


## ----ajfddfussscrssfdjlsa,error=TRUE------------------------------------------
myGenomicIntervals2[indexOfNearest]



## ----ajfddfussscrfdjlecdsa,error=TRUE-----------------------------------------
precedeIndex <- precede(myGenomicIntervals1,myGenomicIntervals2)
followIndex <- follow(myGenomicIntervals1,myGenomicIntervals2)
myGenomicIntervals2[precedeIndex]


## ----ajfddfussscscscrfdjlsa,error=TRUE----------------------------------------
distances <- distanceToNearest(myGenomicIntervals1,myGenomicIntervals2)
distances
mcols(distances)


## ----revruujje,error=TRUE-----------------------------------------------------
library(rtracklayer)
mySicerPeaks <- import.bed(con="data/SicerPeaks.bed")
mySicerPeaks


## ----revkkjkre,error=TRUE-----------------------------------------------------
export.bed(mySicerPeaks, con="moreSicerPeaks.bed")


## ----revkdedkjkaare,echo=FALSE------------------------------------------------
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))


## ----revkdedkjkre,error=TRUE--------------------------------------------------
library(BSgenome.Mmusculus.UCSC.mm10)
subseq(BSgenome.Mmusculus.UCSC.mm10$chr10,1,100)


## ----revkvfydedkjkre,error=TRUE-----------------------------------------------

sicerSeq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, names=mySicerPeaks)
sicerSeq


## ----revkdedkjkkkre,error=TRUE------------------------------------------------
writeXStringSet(sicerSeq,"sicerSeq.fa")

