## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE


## ----installBioconductor, echo=TRUE, eval=FALSE-------------------------------
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## BiocManager::install(version = "3.10")


## ----installPackage, echo=TRUE, eval=FALSE------------------------------------
## BiocManager::install("basecallQC")


## ----version2, echo=TRUE, eval=TRUE-------------------------------------------
BiocManager::version()


## ----updateBioconductor, echo=TRUE, eval=FALSE--------------------------------
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## BiocManager::install(version = "3.10")

