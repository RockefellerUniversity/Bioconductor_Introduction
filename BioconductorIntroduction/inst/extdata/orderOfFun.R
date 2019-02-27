orderOfCourse <- function(path=getwd()){
  rmdFilesOrder <- c(8,7,1,
    6,5,3,
    4,2,9)
  rmdFiles <- c("AlignedDataInBioconductor.Rmd","AlignmentInBioconductor.Rmd","bioconductor_introduction.Rmd",
  "FastQInBioconductor.Rmd","GenomicFeatures_In_Bioconductor.Rmd","GenomicIntervals_In_Bioconductor.Rmd",
  "GenomicScores_In_Bioconductor.Rmd","SequencesInBioconductor.Rmd","Summarising_Scores_In_Bioconductor.Rmd")
  file.path(path,rmdFiles[order(rmdFilesOrder)])
}
