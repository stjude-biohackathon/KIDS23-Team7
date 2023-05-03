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


test <- csv_to_seurat(ExpressionMarker='/home/lead/Akoya/BioHackathon_Samples/Sample1_Group2_MarkerExpression.csv',
                      MetaData = '/home/lead/Akoya/BioHackathon_Samples/Sample1_Group2_MetaData.csv',
                      AssayType='Akoya', FOVName='Sample1_Group2')
















