
csv_to_seurat <- function(ExpressionMarker, MetaData, AssayType='Akoya', FOVName='Image') {
  
  require(Seurat);require(data.table)
  options(Seurat.object.assay.version = "v5")
  
  
  ExpressionMarker <- fread(ExpressionMarker, data.table = F)
  MetaData <- fread(MetaData)
  
  rownames(ExpressionMarker) <- ExpressionMarker[,1]; ExpressionMarker <- ExpressionMarker[,-1]
  
  centroids <- data.frame(x= MetaData[,3],
                          y = MetaData[,2],
                          cell = MetaData[,1])
  
  obj <- CreateSeuratObject(
    counts = as.sparse(t(ExpressionMarker)),
    assay = AssayType,
    meta.data = MetaData
  )
  
  coords <- suppressWarnings(expr = CreateFOV(
    coords = centroids,
    type = 'centroids',
    key = FOVName,
    assay = AssayType
  ))
  
 obj[[FOVName]] <- coords
  
  return(obj)
}

### Maycon

## Seurat object
library(Seurat)
library(BPCells)
library(ggplot2)
library(data.table)

options(future.globals.maxSize = 1e8)
options(Seurat.object.assay.version = "v5")


# Define a function to read in Akoya data and metadata files
inputDir = '/home/lead/Akoya/BioHackathon_Samples'
readSpatialData <- function(inputDir) {
  
  # Load necessary packages
  library(data.table)
  
  # Read in marker expression data
  markerExprFile <- paste0(inputDir, "/Sample1_Group2_MarkerExpression.csv")
  ak <- fread(markerExprFile, data.table = F)
  rownames(ak) <- ak$Row
  ak <- ak[, -1]
  
  # Read in metadata
  metaFile <- paste0(inputDir, "/Sample1_Group2_MetaData.csv")
  ak_meta <- fread(metaFile)
  
  # Return list of data and metadata
  return(list(ak_data = ak, ak_meta = ak_meta))
  
}

SpatialData = readSpatialData(inputDir)







create_seurat_obj <- function(ak, ak_meta) {
  
  library(data.table)
  library(Seurat)
  library(Matrix)
  
  r <- 1/2
  R <- (2 / sqrt(3)) * r
  
  # read in data files
  ak <- fread(ak_file, data.table = F)
  rownames(ak) <- ak$Row
  ak <- ak[, -1]
  ak_meta <- fread(meta_file)
  
  # process spatial data
  ak_meta$y <- ak_meta$Y - min(ak_meta$Y) + 1
  ak_meta$x <- ak_meta$X - min(ak_meta$X) + 1
  ak_meta$y <- ak_meta$y * R * (3/2)
  ak_meta$x <- (ak_meta$x + 1) / 2
  ak_meta$y <- -ak_meta$y
  
  centroids <- data.frame(x = ak_meta$y, y = ak_meta$x, cell = ak_meta$Object_ID)
  
  # create Seurat object
  mtx <- t(ak)
  mtx <- as.sparse(mtx)
  obj <- CreateSeuratObject(
    counts = mtx,
    assay = 'Akoya',
    meta.data = ak_meta
  )
  
  # add image to Seurat object
  coords <- suppressWarnings(expr = CreateFOV(
    coords = centroids,
    type = 'centroids',
    key = 'Test',
    assay = 'Akoya'
  ))
  suppressWarnings(expr = obj[['Test']] <- coords)
  
  # remove DAPI channel and normalize data
  obj <- obj[-1, ]
  obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2)
  
  # sketch data
  obj <- SketchData(
    object = obj,
    ncells = 50000,
    method = "LeverageScore",
    sketched.assay = "sketch", verbose = T
  )
  
  return(obj)
  
}

obj = create_seurat_obj(ak = SpatialData$ak_data, ak_meta = SpatialData$ak_meta)






# Define module for Seurat object processing
processSeurat <- function(inputObj) {
  
  # Set default assay
  DefaultAssay(inputObj) <- "sketch"
  
  # Run standart seurat pipeline 
  inputObj <- ScaleData(inputObj)
  VariableFeatures(inputObj) <- rownames(inputObj)
  
  inputObj <- RunPCA(object = inputObj, npcs = 30)
  inputObj <- RunUMAP(object = inputObj, dims = 1:30)
  inputObj <- FindNeighbors(object = inputObj, dims = 1:30, verbose = FALSE)
  inputObj <- FindClusters(object = inputObj, resolution = 1, algorithm = 2)
  
  # Return processed Seurat object
  return(inputObj)
  
}

obj_processed = processSeurat(inputObj = obj)




# Plot 1
dimPlotModule <- function(seu_object, label = TRUE) {
  
  plot <- DimPlot(seu_object, label = label)
  
  return(plot)
  
}

UMAP_plot = dimPlotModule(obj_processed)






# Plot 2
ImagePlotModule <- function(seu_object, label = TRUE) {
  
  plot <- ImageDimPlot(seu_object, label = label)
  
  return(plot)
  
}

Image_plot = ImagePlotModule(obj_processed)



























































