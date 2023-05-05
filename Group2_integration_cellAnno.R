## testing integrative analysis ##
library(Seurat)


obj1 <- readRDS('/mnt/plummergrp/Felipe/Biohackaton2023/Sample1_Group2_SeuratObj_sketch.rds')
DefaultAssay(obj1) <- 'Akoya'
obj1$Sample <- 'Sample1_Group2'

obj2 <- readRDS('/mnt/plummergrp/Felipe/Biohackaton2023/Sample2_Group2_SeuratObj_sketch.rds')
DefaultAssay(obj2) <- 'Akoya'
obj2$Sample <- 'Sample2_Group2'


g2_obj <- merge(obj1, obj2)

g2_obj <- SketchData(
  object = g2_obj,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch", verbose = T
)
g2_obj <- ScaleData(g2_obj, features = rownames(g2_obj))
VariableFeatures(g2_obj) <- rownames(g2_obj)
g2_obj <- RunPCA(g2_obj, npcs = 30)


library(harmony) 
g2_obj <- RunHarmony(g2_obj, group.by.vars = 'Sample')
g2_obj <- RunUMAP(g2_obj, dims = 1:30, reduction = 'harmony')
g2_obj <- FindNeighbors(g2_obj, dims = 1:30, reduction = 'harmony')
g2_obj <- FindClusters(object = g2_obj, resolution = c(1), algorithm = 2)


ImageDimPlot(g2_obj, fov = c('Sample1_Group2', 'Sample2_Group2'), size = 0.01)




##using sc-type
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

celltype_table <- read_xlsx('/mnt/plummergrp/Felipe/Biohackaton2023/CellType_BrCa_Akoya_Table.xlsx')


gs_list <- strsplit(celltype_table$geneSymbolmore1, ',')
names(gs_list) <- celltype_table$cellName



# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = g2_obj@assays$sketch$scale.data, scaled = TRUE, 
                      gs = gs_list, gs2 = NULL)


# merge by cluster
md <- na.omit(g2_obj@meta.data)
cL_resutls = do.call("rbind", lapply(unique(md$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(md[md$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(md$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 2, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/5] = "Unknown"
print(sctype_scores[,1:3])

g2_obj@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  g2_obj@meta.data$customclassif[g2_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
DimPlot(g2_obj, group.by = 'customclassif', label=T)

ImageDimPlot(g2_obj, fov = c('Sample1_Group2', 'Sample2_Group2'), size = 0.1, group.by = 'customclassif', coord.fixed = F,dark.background = T)

saveRDS(g2_obj, file = '/mnt/plummergrp/Felipe/Biohackaton2023/Group2_SeuratObj_integrated.rds')













