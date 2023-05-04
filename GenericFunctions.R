
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


## QC ##
#..##

## Analysis pipeline #

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

seurat_workflow <- function(obj , NDims = 30, ClusterRes = c(0.3,0.5,1,1.5,2)) {
  
  DefaultAssay(obj) <- "sketch"
  
  if(length(rownames(obj)) < 1000) {
    obj <- ScaleData(obj, features = rownames(obj)) ## using all features to scale if panel is relatively small
    VariableFeatures(obj) <- rownames(obj) # just use all of the features if n < 1000
  } else {
    obj <- FindSpatiallyVariableFeatures(obj) ## Calculates Moran's I for autocorrelation
    obj <- ScaleData(obj, features = SpatiallyVariableFeatures(obj))
  }
  
  obj <- RunPCA(obj, npcs = NDims)
  obj <- RunUMAP(obj, dims = 1:NDims)
  obj <- RunTSNE(obj, dims = 1:NDims)
  obj <- FindNeighbors(obj, dims = 1:NDims)
  obj <- FindClusters(object = obj, resolution = ClusterRes, algorithm = 2)
  
  return(obj)
  
}















































