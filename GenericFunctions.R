
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
  #obj <- RunTSNE(obj, dims = 1:NDims)
  obj <- FindNeighbors(obj, dims = 1:NDims)
  obj <- FindClusters(object = obj, resolution = ClusterRes, algorithm = 2)
  
  return(obj)
  
}



## assumes the celltype table is a xlsx file with 3 columns: tissueType,cellName,geneSymbolmore1

annotate_cell_types <- function(obj, celltype_file_path) {
  # load packages
  require(openxlsx);require(dplyr)
  # load gene set preparation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  celltype_table <- read.xlsx(celltype_file_path)
  
  gs_list <- strsplit(celltype_table$geneSymbolmore1, ',')
  names(gs_list) <- celltype_table$cellName
  
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = obj@assays$sketch$scale.data, scaled = TRUE, 
                        gs = gs_list, gs2 = NULL)
  
  # merge by cluster
  md <- na.omit(obj@meta.data)
  cL_resutls = do.call("rbind", lapply(unique(md$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(md[md$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(md$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 2, wt = scores)  
  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/5] = "Unknown"
  print(sctype_scores[,1:3])
  
  obj@meta.data$CelltypePrediction = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    obj@meta.data$CelltypePrediction[obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
  return(obj)
}











































