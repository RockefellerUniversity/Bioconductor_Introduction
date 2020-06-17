## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE

## ----install, echo=TRUE, eval=FALSE--------------------------------------
## source("https://bioconductor.org/biocLite.R")
## biocLite("basecallQC")

## ----version1, include=FALSE,echo=FALSE, eval=TRUE-----------------------
source("https://bioconductor.org/biocLite.R")

## ----version2, echo=TRUE, eval=TRUE--------------------------------------
biocVersion()

