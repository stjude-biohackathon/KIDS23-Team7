# remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
# remotes::install_github("bnprks/BPCells", quiet = TRUE)
library(ggplot2); library(Seurat); library(dplyr)
source("~/KIDS2/GenericFunctions.R")
options(Seurat.object.assay.version = "v5")


data1 <- fread("/home/lead/Akoya/BioHackathon_Samples/Sample1_Group1_MarkerExpression.csv",data.table = F)[1:10000,]
data2 <- fread("/home/lead/Akoya/BioHackathon_Samples/Sample1_Group1_MetaData.csv", data.table=F)[1:10000,]
write.csv(data1, "./test_MarkerExpression.csv", row.names = F)
write.csv(data2, "./test_MeataData.csv", row.names = F)

# test <- csv_to_seurat(ExpressionMarker="/home/lead/Akoya/BioHackathon_Samples/Sample2_Group2_MarkerExpression.csv",
#                       MetaData = "/home/lead/Akoya/BioHackathon_Samples/Sample2_Group2_MetaData.csv",
#                       AssayType='Akoya', FOVName='Sample2_Group2')

# spare matrix
test@meta.data
test@assays$Akoya$counts
as.data.frame(test[["Akoya"]]$counts)
# test[["Akoya"]]$counts["DAPI",]


# after sketch
test <- readRDS("/home/lead/Akoya/rds/Sample2_Group2_SeuratObj_sketch.rds")
test <- seurat_workflow(test)
test <- annotate_cell_types(test, "/home/lead/Akoya/BioHackathon_Samples/CellType_BrCa_Akoya_Table.xlsx")
# saveRDS(test, "/home/ybae/KIDS2/test_annoated.rds")
test <- readRDS("/home/ybae/KIDS2/test_annoated.rds")

# cell annotated plotly 
celltypes <- unique(test$CelltypePrediction)
plot <- ImageDimPlot(test, fov = c('Sample2_Group2'), size = 0.1, group.by = 'CelltypePrediction', coord.fixed = F, dark.background = F)
p <- ggplotly(plot)

# extract labels from ImageDimPlot
labels_df <- test@meta.data  %>%
  select(CelltypePrediction) %>%
  distinct() 
labels_df$CelltypePrediction
plot$labels
plot$patches
test@

celltype_colors <- Idents(test)$CelltypePredictionColors

# color change 
install.packages("paletteer")
polychrome_pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
 


# create a separate plot with just the labels
ggplot(labels_df, aes(x = 1, y = 10, label = CelltypePrediction)) +
  geom_text(size = 3, color = 'black') +
  theme_void()



Idents(test) <- "CelltypePrediction"
VlnPlot(test, features = c(paste0('nCount_', AssayType),
                           paste0('nFeature_', AssayType)), pt.size = 0)



plotly(RidgePlot(test, features = rownames(test)[1], ncol = 2, layer = 'counts'))
plot <- RidgePlot(test, features = rownames(test)[1], ncol = 2, layer = 'counts')
plot = ggplotly(plot)
plot

library(ggridges)
# create example data
df <- data.frame(
  gene = c(rep("gene1", 100), rep("gene2", 100), rep("gene3", 100)),
  sample = rep(1:100, 3),
  expression = c(rnorm(100, mean = 5, sd = 1), rnorm(100, mean = 10, sd = 1), rnorm(100, mean = 15, sd = 1))
)
df$expression
ggplot(df, aes(x = expression, y = gene, fill = gene)) +
  geom_density_ridges(alpha = 0.8) +
  theme_ridges() +
  xlab("Expression")
# ggplot(test@assays$Akoya@data, aes(x = test@assays$Akoya@data[, "Feature1"], y = rownames(test@assays$Akoya@data))) +
#   geom_density_ridges(alpha = 0.8, scale = 1.5) +
#   theme_ridges() +
#   xlab("Expression")

# plot the ridgeline
typeof(as.data.frame(test[["Akoya"]]$counts))
head(as.data.frame(test[["Akoya"]]$counts))
rownames(test[["Akoya"]]$counts)

test[["Akoya"]]$counts["DAPI",]

ggplot(test[["Akoya"]]$counts, aes(x ="DAPI", y = test[["Akoya"]]$counts["DAPI",])) +
  geom_density_ridges(alpha = 0.8) +
  theme_ridges() 

test[["Akoya"]]$counts["DAPI",]
test[["Akoya"]]$counts


# Convert the sparse matrix to a dense matrix
counts_matrix <- as.matrix(test[["Akoya"]]$counts)
head(as.data.frame(counts_matrix))
show(counts_matrix)
# Create the ggplot
df <- as.data.frame(counts_matrix)
df_long <- melt(df)
dim(df)
df[1,]

df_long[1,]
dim(df_long)
ggplot(data = df_long, aes(x ="DAPI", y = values)) +
  geom_density_ridges(alpha = 0.8) +
  theme_ridges()


## new 
# Convert the dgMatrix to a data.frame
counts_df <- as.data.frame(as.matrix(test[["Akoya"]]$counts))
rownames(counts_df) <- rownames(test)
colnames(counts_df) <- colnames(test)
head(counts_df)

rownames(counts_df)
colnames(counts_df)
# Extract the DAPI column as a vector
dapi <- as.numeric(counts_df["DAPI"])


VlnPlot(test, features="nCount_Akoya", layer='counts')
test[[]]
Layers(test)


Sample1_Group1_SeuratObj_sketch <- readRDS("./test_yj/Sample1_Group1_SeuratObj_sketch.rds")
VlnPlot(Sample1_Group1_SeuratObj_sketch, features="nCount_Akoya", layer='counts')



