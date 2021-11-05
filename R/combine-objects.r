# Initialize the Seurat object with the raw (non-normalized data)
control_1 <- CreateSeuratObject(counts = control_31810_h5,
                       project = "control_1", min.cells = 3, min.features = 300)
control_2 <- CreateSeuratObject(counts = control_31810_h5,
                       project = "control_2", min.cells = 3, min.features = 300)
control_3 <- CreateSeuratObject(counts = control_LE99_h5,
                       project = "control_3", min.cells = 3, min.features = 300)

#DefaultAssay(object = control_3)
#samples <- lapply(ls(pattern = "control_"), get)
#
#res <- lapply(samples, function(x) {
#  
#  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
#  
#  x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
#  
#  x <- RunPCA(x, features = VariableFeatures(object = x))
#  
#  x <- FindNeighbors(x, dims = 1:10)
#  
#  x <- FindClusters(x, resolution = 0.5)
#  
#})

#for (i in seq_along(samples)){ 
#  samples[[i]][["percent.mt"]] <- PercentageFeatureSet(samples[[i]], pattern = "^mt-")
#  
#  samples[[i]] <- SCTransform(samples[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
#  
#  samples[[i]] <- RunPCA(samples[[i]], features = VariableFeatures(object = samples[[i]]))
#  
#  samples[[i]] <- FindNeighbors(samples[[i]], dims = 1:10)
#  
#  samples[[i]] <- FindClusters(samples[[i]], resolution = 0.5)
#  
#}

## control 3
control_3[["percent.mt"]] <- PercentageFeatureSet(control_3, pattern = "^mt-")

control_3 <- SCTransform(control_3, vars.to.regress = "percent.mt", verbose = FALSE)

control_3 <- RunPCA(control_3, features = VariableFeatures(object = control_3))

control_3 <- FindNeighbors(control_3, dims = 1:10)

control_3 <- FindClusters(control_3, resolution = 0.5)

## control 2

control_2[["percent.mt"]] <- PercentageFeatureSet(control_2, pattern = "^mt-")

control_2 <- SCTransform(control_2, vars.to.regress = "percent.mt", verbose = FALSE)

control_2 <- RunPCA(control_2, features = VariableFeatures(object = control_2))

control_2 <- FindNeighbors(control_2, dims = 1:10)

control_2 <- FindClusters(control_2, resolution = 0.5)

## control 1
control_1[["percent.mt"]] <- PercentageFeatureSet(control_1, pattern = "^mt-")

control_1 <- SCTransform(control_1, vars.to.regress = "percent.mt", verbose = FALSE)

control_1 <- RunPCA(control_1, features = VariableFeatures(object = control_1))

control_1 <- FindNeighbors(control_1, dims = 1:10)

control_1 <- FindClusters(control_1, resolution = 0.5)

control.list <- list(control_1, control_2, control_3)


# Integrate data           
# FindIntegrationAnchors() requires list argument for object, so created control.list instead
options(future.globals.maxSize = 4000 * 1024^2)
control.features <- SelectIntegrationFeatures(object.list = control.list, nfeatures = 3000)
control.list <- PrepSCTIntegration(object.list = control.list, anchor.features = control.features, 
                                   verbose = FALSE)

control.anchors <- FindIntegrationAnchors(object.list = control.list, normalization.method = "SCT", 
                                          anchor.features = control.features, verbose = FALSE)
control.integrated <- IntegrateData(anchorset = control.anchors, normalization.method = "SCT", 
                                    verbose = FALSE)
future.seed=TRUE  

# Perform integrated analysis
DefaultAssay(control.integrated) <- "integrated"
# Standard visualization and clustering workflow
control.integrated <- ScaleData(control.integrated, verbose = FALSE)
control.integrated <- RunPCA(control.integrated, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
control.integrated <- RunUMAP(control.integrated, reduction = "pca", dims = 1:20)
control.integrated <- FindNeighbors(control.integrated, reduction = "pca", dims = 1:20)
control.integrated <- FindClusters(control.integrated, resolution = 0.5)

file <- paste(format(Sys.time(),'%Y_%B_%d_%H:%M'),
              "gse_control_integrated.rds", sep = "_")
combined <- file.path(getwd(),
                      "output", "rdata", file)

saveRDS(object = control.integrated, file = combined)
