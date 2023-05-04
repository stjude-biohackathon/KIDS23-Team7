## testing integrative analysis ##
library(Seurat)


obj1 <- readRDS('/mnt/plummergrp/Felipe/Biohackaton2023/Sample1_Group1_SeuratObj_sketch.rds')
DefaultAssay(obj1) <- 'Akoya'
obj1$Sample <- 'Sample1_Group1'

obj2 <- readRDS('/mnt/plummergrp/Felipe/Biohackaton2023/Sample2_Group1_SeuratObj_sketch.rds')
DefaultAssay(obj2) <- 'Akoya'
obj2$Sample <- 'Sample2_Group1'


g1_obj <- merge(obj1, obj2)

g1_obj <- SketchData(
  object = g1_obj,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch", verbose = T
)
g1_obj <- ScaleData(g1_obj, features = rownames(g1_obj))
VariableFeatures(g1_obj) <- rownames(g1_obj)
g1_obj <- RunPCA(g1_obj, npcs = 30)


library(harmony) 
g1_obj <- RunHarmony(g1_obj, group.by.vars = 'Sample')
g1_obj <- RunUMAP(g1_obj, dims = 1:30, reduction = 'harmony')
g1_obj <- FindNeighbors(g1_obj, dims = 1:30, reduction = 'harmony')
g1_obj <- FindClusters(object = g1_obj, resolution = c(1), algorithm = 2)


ImageDimPlot(g1_obj, fov = c('Sample1_Group1', 'Sample2_Group1'), size = 0.01)




##using sc-type
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

celltype_table <- read_xlsx('/mnt/plummergrp/Felipe/Biohackaton2023/CellType_BrCa_Akoya_Table.xlsx')


gs_list <- strsplit(celltype_table$geneSymbolmore1, ',')
names(gs_list) <- celltype_table$cellName



# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = g1_obj@assays$sketch$scale.data, scaled = TRUE, 
                      gs = gs_list, gs2 = NULL)


# merge by cluster
md <- na.omit(g1_obj@meta.data)
cL_resutls = do.call("rbind", lapply(unique(md$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(md[md$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(md$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 2, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/5] = "Unknown"
print(sctype_scores[,1:3])

g1_obj@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  g1_obj@meta.data$customclassif[g1_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
DimPlot(g1_obj, group.by = 'customclassif', label=T)

ImageDimPlot(g1_obj, fov = c('Sample1_Group1', 'Sample2_Group1'), size = 0.1, group.by = 'customclassif', coord.fixed = F,dark.background = F)

saveRDS(g1_obj, file = '/mnt/plummergrp/Felipe/Biohackaton2023/Group1_SeuratObj_integrated.rds')













