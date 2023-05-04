# Testing Felipe seurat object from /home/lead/Akoya/rds

# Packages
library(Seurat)
library(dplyr)
library(BPCells)
library(ggplot2)
library(data.table)
library(readxl)
#library(install_github)
#install.packages(c("Rcpp", "devtools"), dependencies=TRUE)
#require(devtools)
#install_github("awalker89/openxlsx")
library(openxlsx)
#install.packages("plotly")
library(plotly)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("dittoSeq")
library(dittoSeq)

Sample2_Group2_SeuratObj_sketch <- readRDS("/home/lead/Akoya/rds/Sample2_Group2_SeuratObj_sketch.rds")

obj = Sample2_Group2_SeuratObj_sketch


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
  # read cell type table
  celltype_table <- read_xlsx(celltype_file_path)
  
  # prepare cell type input
  gs_list <- strsplit(celltype_table$geneSymbolmore1, ',')
  names(gs_list) <- celltype_table$cellName
  
  # cell type prediction
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


saveRDS(obj_celltype, file = "/home/mmarcao/Sample2_Group2_SeuratObj_sketch.rds")

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