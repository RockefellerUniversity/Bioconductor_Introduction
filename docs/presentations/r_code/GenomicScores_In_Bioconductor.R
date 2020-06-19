## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE


## ----makingData, include=FALSE,cache=FALSE,eval=FALSE-------------------------
## # library(rtracklayer)
## # myRanges <- GRanges("chr1",ranges = IRanges(72811055,72856974))
## # myRanges2 <- c(myRanges,GRanges("chr1",ranges = IRanges(1,100)))
## # mycPeaks <- import.bed("data/Myc_Ch12_1_withInput_Input_Ch12_summits.bed")
## # mycPeaks[mycPeaks %over% myRanges]
## # toRun <- import("~/Downloads/ENCFF940MBK.bigWig",selection=BigWigSelection(myRanges),as="RleList")
## # export(toRun,con="data/TSS_ENCFF940MBK.bedGraph")
## # export(toRun,con="data/TSS_ENCFF940MBK.bw")
## # export.ucsc(toRun,con="data/TSS_ENCFF940MBK.wig",subformat ="wig")


## ----sedfs,eval=TRUE,echo=FALSE-----------------------------------------------
suppressPackageStartupMessages(library(rtracklayer))


## ----ss,eval=FALSE------------------------------------------------------------
## library(rtracklayer)


## ----c,eval=TRUE--------------------------------------------------------------
myBedG <- import.bedGraph("data/TSS_ENCFF940MBK.bedGraph")


## ----b,eval=TRUE--------------------------------------------------------------
myBedG[1:3]


## ----d,eval=TRUE--------------------------------------------------------------
strand(myBedG)


## ----e,eval=TRUE--------------------------------------------------------------
myBigWig <- import.bw("data/TSS_ENCFF940MBK.bw")


## ----f,eval=TRUE--------------------------------------------------------------
myBigWig[1:3]


## ----g,eval=TRUE--------------------------------------------------------------
myGRanges <- GRanges("chr1",IRanges(72823698,72824485))
filteredBigWig <- myBigWig[myBigWig %over% myGRanges]
filteredBigWig[1:3]


## ----h,eval=TRUE--------------------------------------------------------------
myBigWig <- import.bw("data/TSS_ENCFF940MBK.bw",
                      as = "RleList")
class(myBigWig)


## ----hqqk,eval=TRUE-----------------------------------------------------------
mycPeaks <- import.bed("data/Myc_Ch12_1_withInput_Input_Ch12_summits.bed")
seqnames(mycPeaks)


## ----hqq,eval=TRUE------------------------------------------------------------
myNumbers <- c(0,0,0,0,0,1,1,1,0,0,0,0,0)
Rle(myNumbers)


## ----hss,eval=TRUE------------------------------------------------------------
myNumbers2 <- c(0,0,0,0,0,1,1,1,2,2,2,2,2)
chr1Scores <- Rle(myNumbers)
chr2Scores <- Rle(myNumbers2)
myRleList <- RleList(chr1=chr1Scores,chr2=chr2Scores)
myRleList


## ----i,eval=TRUE--------------------------------------------------------------
myBigWig[1:2]


## ----j,eval=TRUE--------------------------------------------------------------
chr1_rle <- myBigWig$chr1
# Or
chr1_rle <- myBigWig[["chr1"]]
chr1_rle


## ----k,eval=TRUE--------------------------------------------------------------
chr1_rle[1:10]


## ----l,eval=TRUE--------------------------------------------------------------
chr1_rle[1:10] <- 100
chr1_rle


## ----m,eval=TRUE--------------------------------------------------------------
rleAsVector <- as.vector(chr1_rle[1:10])
rleAsVector


## ----n,eval=TRUE--------------------------------------------------------------
rleAsDF <- as.data.frame(chr1_rle[1:10])
rleAsDF


## ----o,eval=TRUE--------------------------------------------------------------
chr1_rle+1000
(chr1_rle+1000)*10


## ----o1,eval=TRUE-------------------------------------------------------------
chr1_rle < 10


## ----o3,eval=TRUE-------------------------------------------------------------
chr1_rle[chr1_rle < 10] <- 0
chr1_rle


## ----o2,eval=TRUE-------------------------------------------------------------
mean(chr1_rle)
max(chr1_rle)
sum(chr1_rle)


## ----o4,eval=TRUE-------------------------------------------------------------
myBigWig <- import.bw("data/TSS_ENCFF940MBK.bw",as="RleList")
myBigWig+10


## ----o4s,eval=TRUE------------------------------------------------------------
chromosomeMax <- max(myBigWig)
chromosomeMax[1:4]


## ----oa4,eval=TRUE------------------------------------------------------------
myRanges <- GRanges("chr1",ranges = IRanges(72811055,72856974))
mycPeaks <- import.bed("data/Myc_Ch12_1_withInput_Input_Ch12_summits.bed")
mycPeaks <- resize(mycPeaks,50,fix="center")
newMycPeaks <- mycPeaks[mycPeaks %over% myRanges]
newMycPeaks


## ----v,eval=TRUE--------------------------------------------------------------
rleOverGranges <- myBigWig[newMycPeaks]
rleOverGranges


## ----w,eval=TRUE--------------------------------------------------------------
sum(rleOverGranges)


## ----p1,eval=TRUE-------------------------------------------------------------
myRleList <- RleList(chr1=chr1_rle)
myRleList


## ----q,eval=TRUE--------------------------------------------------------------
export.bw(myRleList,con="chr1_Myc.bw")


## ----x,eval=TRUE--------------------------------------------------------------
newMycPeaks


## ----z,eval=TRUE--------------------------------------------------------------

mySelection <- BigWigSelection(newMycPeaks)
import.bw("data/TSS_ENCFF940MBK.bw", 
          selection=mySelection, 
          as="RleList")

