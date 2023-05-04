
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






## Maycon/Felipe P. - begin
## Seurat object
library(Seurat)
library(dplyr)
library(BPCells)
library(ggplot2)
library(data.table)
library(readxl)
library(install_github)
install.packages(c("Rcpp", "devtools"), dependencies=TRUE)
require(devtools)
install_github("awalker89/openxlsx")
library(openxlsx)

options(future.globals.maxSize = 1e8)
options(Seurat.object.assay.version = "v5")


# Define a function to read in Akoya data and metadata files
inputDir = '/home/lead/Akoya/BioHackathon_Samples'
readSpatialData <- function(inputDir) {
  
  # Load necessary packages
  library(data.table)
  
  # Read in marker expression data
  markerExprFile <- paste0(inputDir, "/Sample1_Group2_MarkerExpression.csv")
  matrix_exp <- fread(markerExprFile, data.table = F)
  rownames(matrix_exp) <- matrix_exp$Row
  matrix_exp <- matrix_exp[, -1]
  
  # Read in metadata
  metaFile <- paste0(inputDir, "/Sample1_Group2_MetaData.csv")
  meta_data <- fread(metaFile)
  
  # Return list of data and metadata
  return(list(matrix_exp = matrix_exp, meta_data = meta_data))
  
}

SpatialData = readSpatialData(inputDir)



matrix_exp = SpatialData$matrix_exp
meta_data = SpatialData$meta_data


# Create seurat object 
create_seurat_obj <- function(matrix_exp, meta_data) {
  
  # Re arrange pixel positions to match tissue
  meta_data$y <- meta_data$Y - min(meta_data$Y) + 1
  meta_data$x <- meta_data$X - min(meta_data$X) + 1
  meta_data$y <- meta_data$y * R * (3/2)
  meta_data$x <- (meta_data$x + 1) / 2
  meta_data$y <- -meta_data$y
  
  centroids <- data.frame(x = meta_data$y, y = meta_data$x, cell = meta_data$Object_ID)
  
  # create Seurat object
  mtx <- t(matrix_exp)
  mtx <- as.sparse(mtx)
  obj <- CreateSeuratObject(
    counts = mtx,
    assay = 'matrix_expoya',
    meta.data = meta_data
  )
  
  # add image to Seurat object
  coords <- suppressWarnings(expr = CreateFOV(
    coords = centroids,
    type = 'centroids',
    key = 'Test',
    assay = 'matrix_expoya'
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

obj = create_seurat_obj(matrix_exp = SpatialData$matrix_exp, meta_data = SpatialData$meta_data)






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



obj_processed@meta.data$seurat_clusters
# Plot 1
dimPlotModule <- function(seu_object, label = TRUE, group_by_feature) {
  
  plot <- DimPlot(seu_object, label = label, group.by = group_by_feature)
  
  return(plot)
  
}

UMAP_plot = dimPlotModule(obj_processed, group_by_feature = "seurat_clusters")

UMAP_plot

DimPlot()


# Plot 2
ImagePlotModule <- function(seu_object, group_by_feature) {
  
  plot <- ImageDimPlot(seu_object, group.by = group_by_feature)
  
  return(plot)
  
}

Image_plot = ImagePlotModule(obj_processed, group_by_feature = "seurat_clusters")
Image_plot






# Cell type annotation

seurat_obj = obj_processed
celltype_file_path = "/home/lead/Akoya/BioHackathon_Samples/CellType_BrCa_Akoya_Table.xlsx"

annotate_cell_types <- function(seurat_obj, celltype_file_path) {
  # load packages
  library(openxlsx)
  # load gene set preparation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  celltype_table <- read_xlsx(celltype_file_path)
  
  gs_list <- strsplit(celltype_table$geneSymbolmore1, ',')
  names(gs_list) <- celltype_table$cellName
  
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = seurat_obj@assays$sketch$scale.data, scaled = TRUE, 
                        gs = gs_list, gs2 = NULL)
  
  # merge by cluster
  md <- na.omit(seurat_obj@meta.data)
  cL_resutls = do.call("rbind", lapply(unique(md$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(md[md$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(md$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 2, wt = scores)  
  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/5] = "Unknown"
  print(sctype_scores[,1:3])
  
  seurat_obj@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    seurat_obj@meta.data$customclassif[seurat_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
  return(seurat_obj)
}

obj_celltype <- annotate_cell_types(obj_processed, "/home/lead/Akoya/BioHackathon_Samples/CellType_BrCa_Akoya_Table.xlsx")



UMAP_plot_celltype = dimPlotModule(obj_celltype, group_by_feature = "customclassif") 
UMAP_plot_celltype

Image_plot_celltype = ImagePlotModule(obj_celltype, group_by_feature = "customclassif")
Image_plot_celltype




### Interactive plotting

# UMAP category labels - interactively plot
dimPlotModule_interactive_cat <- function(seu_object, label = TRUE, group_by_feature) {
  Idents(seu_object) = "customclassif"
  plot <- DimPlot(seu_object, label = label, group.by = group_by_feature)
  plot = HoverLocator(plot = plot, information = FetchData(seu_object, vars = c("seurat_clusters", "customclassif", "nCount_matrix_expoya", "nFeature_matrix_expoya")))
  
  return(plot)
  
}


umap_int_cat = dimPlotModule_interactive_cat(seu_object = obj_celltype, 
                         group_by_feature = "customclassif")
umap_int_cat


# UMAP gene expression - interactively plot
dimPlotModule_interactive_gen <- function(seu_object, label = TRUE, feature_molecule, color_min, color_max) {
  Idents(seu_object) = "customclassif"
  plot <- FeaturePlot(seu_object, label = label, features = feature_molecule, cols = c(color_min, color_max)) 
  plot = HoverLocator(plot = plot, information = FetchData(seu_object, vars = c("seurat_clusters", "customclassif", "nCount_matrix_expoya", "nFeature_matrix_expoya")))
  
  return(plot)
  
}


umap_int_gen = dimPlotModule_interactive_gen(seu_object = obj_celltype,
                                             feature_molecule = rownames(obj_celltype)[1], color_min = "white", color_max = "red")
umap_int_gen




## Maycon/Felipe P. - end








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
















































































