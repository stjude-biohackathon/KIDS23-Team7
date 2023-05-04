## Seurat object
# remotes::install_github("bnprks/BPCells")
# library(BPCells)
library(Seurat)
library(ggplot2)
library(data.table)

options(future.globals.maxSize = 1e8)
# options(Seurat.object.assay.version = "v5")


csv_to_seurat <- function(ExpressionMarker, MetaData, AssayType='Akoya', FOVName='Image') {
  # return(ExpressionMarker)
  require(Seurat);require(data.table)
  options(Seurat.object.assay.version = "v5")
  # return("test")
  # print("csv to seurat working")
  
  # ExpressionMarker <- fread(ExpressionMarker, data.table = F)
  # MetaData <- fread(MetaData)
  
  ExpressionMarker <- read.csv(ExpressionMarker)
  MetaData <- read.csv(MetaData)
  # return("test")
  
  rownames(ExpressionMarker) <- ExpressionMarker[,1]; ExpressionMarker <- ExpressionMarker[,-1]
  
  
  centroids <- data.frame(x= MetaData[,3],
                          y = MetaData[,2],
                          cell = MetaData[,1])
  # AssayType = "Akoya"
  # print(centroids)
  # print(AssayType)
  # print(head(ExpressionMarker))
  # print(head(as.sparse(t(ExpressionMarker))))
  # print(head(MetaData))
  # print(FOVName)
  # return('done')
  
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
  # return("test")
  obj[[FOVName]] <- coords
  
  return(obj)
}

norm_and_sketch <- function(obj, AssayType = 'Akoya', NCells = 50000) {
  
  #Normalization
  if(AssayType == 'Akoya') {
    obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2)
  } 
  
  if(AssayType == 'Nanostring|Vizgen|Spatial|Xenium') {
    obj <- SCTransform(obj, assay = AssayType)
  }
  
  #Sketching
  obj <- SketchData(
    object = obj,
    ncells = NCells,
    method = "LeverageScore",
    sketched.assay = "sketch", verbose = T
  )
  
  return(obj)
}
