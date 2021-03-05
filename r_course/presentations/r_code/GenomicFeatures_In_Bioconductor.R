params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE


## ----results='asis',include=TRUE,echo=FALSE-----------------------------------
if(params$isSlides == "yes"){
  cat("
class: inverse, center, middle

# Genomic Features

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("

## Genomic Features

---
"    
  )
  
}


## ----setssssup, eval=FALSE,echo=FALSE-----------------------------------------
## library(TxDb.Hsapiens.UCSC.hg19.knownGene)
## export.gff(TxDb.Hsapiens.UCSC.hg19.knownGene,con = "data/test.gff")


## ----setssssuxxp, eval=FALSE--------------------------------------------------
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## 
## BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")


## ----sedfs,eval=TRUE,echo=FALSE-----------------------------------------------
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))


## ----ssaa,eval=FALSE----------------------------------------------------------
## library(TxDb.Hsapiens.UCSC.hg19.knownGene)


## ----sss,eval=TRUE------------------------------------------------------------
class(TxDb.Hsapiens.UCSC.hg19.knownGene)


## ----sssk,eval=TRUE-----------------------------------------------------------
TxDb.Hsapiens.UCSC.hg19.knownGene


## ----ssslkk,eval=TRUE---------------------------------------------------------
myGenes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
myGenes


## ----ssslkzassk1,eval=TRUE----------------------------------------------------
myExons <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
myExons


## ----ssslkzaak1,eval=TRUE-----------------------------------------------------
myTranscripts <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
myTranscripts


## ----ssslkzak1z,eval=TRUE-----------------------------------------------------
myCDS <- cds(TxDb.Hsapiens.UCSC.hg19.knownGene)
myCDS


## ----ssslkzzk,eval=TRUE-------------------------------------------------------
myTranscripts <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                             columns=c("gene_id","tx_id"))
myTranscripts[1:2]


## ----ssslfdkzzkjj,eval=TRUE---------------------------------------------------
myPromoters <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                             upstream=2000,downstream=50)
myPromoters[1:2]


## ----ssslkzzzjzkjj,eval=TRUE--------------------------------------------------
transcriptByGenes <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                   by="gene")
transcriptByGenes[1:2]


## ----ssslkzzszzkjj,eval=TRUE--------------------------------------------------
# transcriptByGenes$1 or
transcriptByGenes[[1]]


## ----ssslkzzzzfkjj,eval=TRUE--------------------------------------------------
exonsByTranscript <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                   by="tx")
exonsByTranscript[1:2]


## ----ssslkzzzzfkkjkjj,eval=TRUE-----------------------------------------------
myCustomTxDb <- makeTxDbFromGFF("data/Xkr4.gtf")
class(myCustomTxDb)


## ----ssslkzzzzfkaakjkjaj,eval=TRUE--------------------------------------------
myCustomTxDb


## ----sssleekzzzzfkjj,eval=TRUE------------------------------------------------
genes(myCustomTxDb)


## ----ssaaslkzzzzfkjj,eval=TRUE------------------------------------------------
exonsBy(myCustomTxDb,by="gene")


## ----skuysslkzzzzfkjj,eval=TRUE-----------------------------------------------
library(rtracklayer)
availableGenomes <- ucscGenomes()
availableGenomes[1:4,]


## ----ssslkzzzzfkssjj,eval=F,cache=FALSE,message=FALSE-------------------------
## hg18TxDb <- makeTxDbFromUCSC(genome="hg18")
## hg18TxDb


## ----ssslkzzzzfkssjj2, eval=TRUE, echo=F--------------------------------------

#saveDb(hg18TxDb, file="data/hg18txdb.sqlite")
hg18TxDb <- loadDb("data/hg18txdb.sqlite")
#load(file="data/hg18txdb.RData")
hg18TxDb


## ----ssslkzzzzddfkjj,eval=TRUE,cache=TRUE,message=FALSE,warning=FALSE---------
hg18Promoters <- promoters(hg18TxDb,2000,50)


## ----ssslkzzzzdddsdsfkjj,eval=TRUE,cache=TRUE,message=FALSE,warning=FALSE-----
export.gff(myCustomTxDb,con="customTxDbb.gff",format="gff")
export.gff(myCustomTxDb,con="customTxDbb.gtf",format="gtf")


## ----ssslkzzzzsdsddfkjj,eval=TRUE,cache=TRUE,message=FALSE,warning=FALSE------
transcriptByGenes <- exonsBy(hg18TxDb,by="gene")
length(transcriptByGenes)


## ----ssslkzzzzddcccfkjj,eval=TRUE,cache=TRUE,message=FALSE,warning=FALSE------
transcriptNumberPerGene <- lengths(transcriptByGenes)
transcriptNumberPerGene[1:5]


## ----ssslkzzxxzssszddfkjj,eval=TRUE,cache=TRUE,message=FALSE,warning=FALSE----
transcript_Lens <- transcriptLengths(hg18TxDb)
transcript_Lens[1:5,]


## ----ssslkzzzzddsddfkjj,eval=TRUE,cache=TRUE,message=FALSE,warning=FALSE------
library(BSgenome.Hsapiens.UCSC.hg19)
hg19TransSeq <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg19, 
                                      transcripts=TxDb.Hsapiens.UCSC.hg19.knownGene)
hg19TransSeq


## ----ssslkdszzzzddfkjj,eval=TRUE,cache=TRUE,message=FALSE,warning=FALSE-------
writeXStringSet(hg19TransSeq,"myTranscriptSequences.fa")


## -----------------------------------------------------------------------------
library(GenomeInfoDb)


## -----------------------------------------------------------------------------
allMappings <- genomeStyles()
names(allMappings)


## -----------------------------------------------------------------------------
#allMappings$Homo_sapiens or
allMappings[["Homo_sapiens"]]


## -----------------------------------------------------------------------------
myGenes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
seqlevelsStyle(myGenes)
myGenes[1:2,]



## -----------------------------------------------------------------------------
seqlevelsStyle(myGenes) <- "Ensembl"
myGenes[1:2,]



## ----ssslkdszzzdsdzddfkjjsds,eval=FALSE,echo=TRUE,cache=TRUE,message=FALSE,warning=FALSE----
## BiocManager::install("org.Hs.eg.db")
## library(org.Hs.eg.db)
## class(org.Hs.eg.db)


## ----ssslkdszzsdzzddfkjjsds,eval=TRUE,echo=FALSE,cache=FALSE,message=FALSE,warning=FALSE----
suppressPackageStartupMessages(library(org.Hs.eg.db))
class(org.Hs.eg.db)


## ---- echo=TRUE, eval=TRUE----------------------------------------------------
columns(org.Hs.eg.db)


## ---- echo=TRUE, eval=TRUE----------------------------------------------------
help(GENENAME)


## ---- echo=TRUE, eval=TRUE----------------------------------------------------
keytypes(org.Hs.eg.db)


## ---- echo=TRUE, eval=TRUE----------------------------------------------------
keys(org.Hs.eg.db, keytype="SYMBOL")[1:10]


## ---- echo=TRUE, eval=TRUE,message=FALSE,warning=FALSE------------------------
select(org.Hs.eg.db, keys = "A1BG", keytype = "SYMBOL", 
       columns = c("SYMBOL", "GENENAME", "CHR") )


## ----sqsq, echo=TRUE, eval=TRUE, dependson="ssslkzzzzfkssjj"------------------
geneLocations <- genes(hg18TxDb)
geneLocations


## ----sqsaq, echo=TRUE, eval=TRUE, dependson="ssslkzzzzfkssjj",message=FALSE,warning=FALSE----
IDs <- geneLocations$gene_id
myTable <- select(org.Hs.eg.db, keys = IDs, keytype = "ENTREZID",
                  columns = c("SYMBOL", "GENENAME", "ENTREZID") )
myTable[1:2,]


