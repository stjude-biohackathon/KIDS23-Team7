
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

obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2)
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
















































































