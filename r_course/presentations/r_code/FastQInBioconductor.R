## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE



## ----setup2, include=FALSE,eval=FALSE,echo=FALSE------------------------------
## library(ShortRead)
## temp <- readFastq("~/Projects/Results/RNAseqPipeTest/FirstTest/FQs/ENCFF000CXH.fastq.gz")
## fastqSample <- temp[1:100000]
## writeFastq(fastqSample,file = "~/Projects/Software/Github/RUBioconductor_Introduction/r_course/data/sampled_ENCFF000CXH.fastq.gz",mode = "w")


## ----load, echo=TRUE,eval=FALSE-----------------------------------------------
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## 
## BiocManager::install("ShortRead")


## ----load1, echo=FALSE,eval=TRUE----------------------------------------------
suppressPackageStartupMessages(library(ShortRead))


## ----setup2112,eval=TRUE,echo=TRUE--------------------------------------------
library(ShortRead)
fastQ <- readFastq("data/sampled_ENCFF000CXH.fastq.gz")
class(fastQ)


## ----setup22------------------------------------------------------------------
fastQ


## ----setup23------------------------------------------------------------------
length(fastQ)


## ----setup2q3-----------------------------------------------------------------
readLengths <- width(fastQ)
readLengths[1:10]


## ----setup2q3a----------------------------------------------------------------
fastQ[1:10]


## ----setup231zzw--------------------------------------------------------------
sequenceOfReads <- sread(fastQ)
class(sequenceOfReads)


## ----setup231ccw3drr----------------------------------------------------------
sequenceOfReads 


## ----setup231w3dddrr----------------------------------------------------------
alpFreq <- alphabetFrequency(sequenceOfReads)
alpFreq[1:2,]


## ----setup23ew2E--------------------------------------------------------------
idsOfReads <- id(fastQ)
class(idsOfReads)
idsOfReads[1:2]


## ----setup23dd3re-------------------------------------------------------------
Ids <- as.character(idsOfReads)
Ids[1:4]


## ----setup23s3wwdre-----------------------------------------------------------
quals <- quality(fastQ)
class(quals)


## ----setup23s3reed------------------------------------------------------------
quals


## ----setup23we4---------------------------------------------------------------
qualityEncoding <- encoding(quals)
qualityEncoding


## ----setup23sswejjdccd4-------------------------------------------------------
quals[1]
quals[[1]]


## ----setup23sswedccdfds4------------------------------------------------------
toTranslateList <- strsplit(as.character(quals[[1]]),"")
toTranslate <- unlist(toTranslateList)
toTranslate


## ----setup23sswedcqcjjd4------------------------------------------------------
qualityEncoding[toTranslate]


## ----setup23sswedccjjwwd4-----------------------------------------------------
readScores<- alphabetScore(quals)
readScores[1]
sum(qualityEncoding[toTranslate])


## ----setup23sswesqdccjjwwd4---------------------------------------------------
matrixOfQualities <- as(quals,"matrix")
rowSums(matrixOfQualities)[1]


## ----setup23sswessdccrrjjd4---------------------------------------------------
alpByCyc <- alphabetByCycle(sequenceOfReads)
alpByCyc[1:4,1:15]


## ----setup23sswedccraarjjd4---------------------------------------------------
qualsByCyc <- alphabetByCycle(quals)
qualsByCyc[1:4,1:15]


## ----setup23sswedccrssrjjd4,eval=FALSE,echo=TRUE------------------------------
## readOccurence <- table(sequenceOfReads)
## sort(readOccurence,decreasing = TRUE)[1:2]


## ----setup23sswedccrssrdjjd44,eval=TRUE,echo=FALSE----------------------------
readOccurence <- BiocGenerics::table(sequenceOfReads)

## ----setup23sswedccrssrdjjd44_2,eval=TRUE,echo=FALSE--------------------------
as.data.frame(sort(readOccurence,decreasing = TRUE)[1:2])


## ----setup23sswedsqccrrjjd4---------------------------------------------------
duplicates <- srduplicated(fastQ)
duplicates[1:3]


## ----setup23sswedcwdcrrjjd4---------------------------------------------------
table(duplicates)


## ----setup2355a---------------------------------------------------------------
my_QA <- qa("data/sampled_ENCFF000CXH.fastq.gz")
my_QA


## ----setup236a----------------------------------------------------------------
myReport <- report(my_QA)
myReport


## -----------------------------------------------------------------------------
report(my_QA, dest="QC_report")


## ----setup23q7,eval=FALSE-----------------------------------------------------
## browseURL(myReport)


## ----setup2ee37---------------------------------------------------------------
TrimmedFastq <- trimTails(fastQ,20,"5")
TrimmedFastq


## ----setup237ede,eval=FALSE---------------------------------------------------
## writeFastq(TrimmedFastq,"myTrimmed_Fastq.fastq.gz")


## ----setup23kk7ede,echo=FALSE-------------------------------------------------
unlink("myTrimmed_Fastq.fastq.gz")
writeFastq(TrimmedFastq,"myTrimmed_Fastq.fastq.gz",mode="w")


## ---- warnings=F, message=FALSE-----------------------------------------------
library(Rfastp)

rfastp_report <- rfastp(read1 = "data/sampled_ENCFF000CXH.fastq.gz", outputFastq ="data/Rfastp_ENCFF000CXH")



## -----------------------------------------------------------------------------
dfsummary <- qcSummary(rfastp_report)
dfsummary


## -----------------------------------------------------------------------------
curvePlot(rfastp_report, curve="quality_curves")



## -----------------------------------------------------------------------------
curvePlot(rfastp_report, curve="content_curves")



## ----setup237fef--------------------------------------------------------------
sampleToRead <- FastqSampler("data/sampled_ENCFF000CXH.fastq.gz",
                             n=100)
yield(sampleToRead)


## ----setup23kk7ssa------------------------------------------------------------
sampleToRead <- FastqStreamer("data/sampled_ENCFF000CXH.fastq.gz",
                              n=100)
first100Reads <- yield(sampleToRead)
second100Reads <- yield(sampleToRead)


## ----setup2kk37---------------------------------------------------------------
fq <- FastqStreamer("data/sampled_ENCFF000CXH.fastq.gz",
                   n=25000)
while (length(fq_stream <- yield(fq)) > 0) {
    print(length(fq_stream ))
}

## ----closeconn, echo=F--------------------------------------------------------
close(fq)

