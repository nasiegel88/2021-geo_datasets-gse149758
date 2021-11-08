# Initialize the Seurat object with the raw (non-normalized data)

source("R/download-data.r", echo = TRUE)

control_1 <- CreateSeuratObject(counts = control_LD09_h5,
                                project = "control_1", min.cells = 3,
                                min.features = 300)
control_2 <- CreateSeuratObject(counts = control_31810_h5,
                                project = "control_2", min.cells = 3,
                                min.features = 300)
control_3 <- CreateSeuratObject(counts = control_LE99_h5,
                                project = "control_3", min.cells = 3,
                                min.features = 300)

merglung <- merge(x = control_1, y = c(control_2, control_3),
                  add.cell.ids = c("control-1", "control-2", "control-3"),
                  project = "control_animals")

file <- paste(format(Sys.time(),'%Y_%B_%d_%H:%M'),
              "azimuth_input.rds", sep = "_")
full <- file.path(getwd(),
                      "output", "rdata", file)

saveRDS(object = merglung, file = full)