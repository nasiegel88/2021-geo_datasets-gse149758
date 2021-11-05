# install sccatch from github
if(!require(Seurat)){
  install.packages("Seurat", repos = "http://cran.us.r-project.org")
  library(Seurat)
}
if(!require(scCATCH)){
  devtools::install_github("ZJUFanLab/scCATCH", force = TRUE)
  library(scCATCH)
}
if(!require(ulimit)){
  devtools::install_github("krlmlr/ulimit", force = TRUE)
  library(ulimit)
}

int <- readRDS("./output/rdata/2021_November_05_12:32_gse_control_integrated.RData")

require('scCATCH')
clu_markers_combined <- findmarkergenes(object = int,
                                        species = 'Human',
                                        cluster = 'All',
                                        match_CellMatch = FALSE,
                                        cancer = NULL,
                                        tissue = NULL,
                                        cell_min_pct = 0.25,
                                        logfc = 0.25,
                                        pvalue = 0.05)
set.seed(1)
out_file <- file.path(getwd(), "output", "rdata",
          paste(format(Sys.time(),
                       '%Y_%B_%d_%H:%M'), "sccatch-output.RData", sep = "_"))
save.image(file = out_file)