### Testing loading Akoya sample for biohackaton project ###

library(data.table)
library(Seurat)
options(Seurat.object.assay.version = 'v5')


ak <- fread('/home/lead/Akoya/S18_13562_B1_Scan1/S18_13562_B1_Scan1_Marker_Expression.csv', data.table = F)
rownames(ak) <- ak$Row; ak <- ak[, -1]
ak_meta <- fread('/home/lead/Akoya/S18_13562_B1_Scan1/S18_13562_B1_Scan1_metadata.csv')

##Re arrange pixel positions to match tissue
r <- 1/2
R <- (2 / sqrt(3)) * r
ak_meta$y <- ak_meta$Y - min(ak_meta$Y) + 1
ak_meta$x <- ak_meta$X - min(ak_meta$X) + 1
ak_meta$y <- ak_meta$y * R * (3/2)
ak_meta$x <- (ak_meta$x + 1) / 2
ak_meta$y <- -ak_meta$y

centroids <- data.frame(x= ak_meta$y,
                        y = ak_meta$x,
                        cell = ak_meta$Object_ID)


mtx <- t(ak)
mtx <- as.sparse(mtx)

obj <- CreateSeuratObject(
  counts = mtx,
  assay = 'Akoya',
  meta.data = ak_meta
)

coords <- suppressWarnings(expr = CreateFOV(
  coords = centroids,
  type = 'centroids',
  key = 'Test',
  assay = 'Akoya'
))

suppressWarnings(expr = obj[['Test']] <- coords) # add image to seurat object

obj <- obj[-1, ] ## remove DAPI channel 

obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2)
obj <- ScaleData(obj)
VariableFeatures(obj) <- rownames(obj)  # since the panel is small, treat all features as variable.
obj <- RunPCA(object = obj, npcs = 30)
obj <- RunUMAP(object = obj, dims = 1:30)
obj <- FindNeighbors(object = obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(object = obj, resolution = 1, algorithm = 2)
DimPlot(obj, label = T)
ImageDimPlot(obj, size = 0.1, cols = 'glasbey')


