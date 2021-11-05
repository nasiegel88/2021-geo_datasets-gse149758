if(!require(renv)){
  install.packages("renv")
  library('renv')
}
# bioconductor
if(!require(BiocManager)){
  install.packages("BiocManager")
  library('BiocManager')
}
if(!require(rmarkdown)){
  install.packages("rmarkdown", dep = TRUE)
  library('rmarkdown')
}
if(!require(tinytex)){
  install.packages("tinytex")
  tinytex::install_tinytex()
}
if(!require(devtools)){
  install.packages("devtools")
  library('devtools')
}
if(!require(remotes)){
  install.packages('remotes')
  library('remotes')
}
if(!require(sctransform)){
  devtools::install_github(repo = 'ChristophH/sctransform', force = TRUE)
  library(sctransform)
}
if(!require(Seurat)){
  install.packages('Seurat')
  library(Seurat)
}
if(!require(scCATCH)){
  devtools::install_github('ZJUFanLab/scCATCH', force = TRUE)
  library('scCATCH')
}
if(!require(patchwork)){
  install.packages("patchwork")
  library('patchwork')
}
if(!require(dplyr)){
  install.packages("dplyr")
  library('dplyr')
}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library('ggplot2')
}
if(!require(kableExtra)){
  install.packages("kableExtra")
  library('kableExtra')
}
if(!require(purrr)){
  install.packages("purrr")
  library('purrr')
}
if(!require(gdtools)){
  install.packages("gdtools")
  library('gdtools')
}
if(!require(flextable)){
  install.packages("flextable")
  library('flextable')
}
if(!require(pracma)){
  install.packages("pracma")
  library('pracma')
}
if(!require(stringr)){
  install.packages("stringr")
  library('stringr')
}
if(!require(data.table)){
  install.packages("data.table")
  library('data.table')
}
# Bioconductor
if(!require(multtest)){
  BiocManager::install("multtest")
  library('multtest')
}
if(!require(metap)){
  BiocManager::install("metap")
  library('metap')
}
if(!require(escape)){
  BiocManager::install("escape")
  library('escape')
}
if(!require(dittoSeq)){
  BiocManager::install("dittoSeq")
  library('dittoSeq')
}
if(!require(SingleCellExperiment)){
  BiocManager::install("SingleCellExperiment")
  library('SingleCellExperiment')
}
if(!require(pcaMethods)){
  BiocManager::install("pcaMethods")
  library(pcaMethods)
}
if(!require(rhdf5)){
  BiocManager::install("rhdf5")
  library('rhdf5')
}
if(!require(hdf5r)){
  install.packages('hdf5r')
  library('hdf5r')
}
