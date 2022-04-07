#source file of DBEC preprocessing
#background subtration of BD Rhapsody WTA data based on bimodal distribution of lognormal raw read counts 
#of each gene by using mclust package
#Written by Shigeyuki Shichino, ver0.2 20190603
library(reticulate)
reticulate::use_python(
  python = "C:/Users/AokiMPM/Anaconda3/python.exe",
  required = TRUE
)
suppressMessages(library(BiocParallel))
suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(Matrix))
suppressMessages(library(mclust))
suppressMessages(library(Seurat))
suppressMessages(library(recommenderlab))
suppressMessages(library(MASS))
suppressMessages(library(future.apply))
suppressMessages(library(rDBEC))
suppressMessages(library(ReductionWrappers))
set.seed(42)
options(future.globals.maxSize= 100*1024^3)

