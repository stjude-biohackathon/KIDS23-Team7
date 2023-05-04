# remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
# remotes::install_github("bnprks/BPCells", quiet = TRUE)
library(ggplot2); library(Seurat); library(dplyr)
source("csv_to_seurat.R")
options(Seurat.object.assay.version = "v5")


data1 <- fread("/home/lead/Akoya/BioHackathon_Samples/Sample1_Group1_MarkerExpression.csv",data.table = F)[1:10000,]
data2 <- fread("/home/lead/Akoya/BioHackathon_Samples/Sample1_Group1_MetaData.csv", data.table=F)[1:10000,]
write.csv(data1, "./test_MarkerExpression.csv", row.names = F)
write.csv(data2, "./test_MeataData.csv", row.names = F)
# 
data1
data2
# test <- csv_to_seurat(ExpressionMarker="/home/lead/Akoya/BioHackathon_Samples/Sample2_Group2_MarkerExpression.csv",
#                       MetaData = "/home/lead/Akoya/BioHackathon_Samples/Sample2_Group2_MetaData.csv",
#                       AssayType='Akoya', FOVName='Sample2_Group2')

# spare matrix
test
test@meta.data
test@assays$Akoya$counts
test[["Akoya"]]$counts

RidgePlot(test, features = rownames(test)[1:10], ncol = 2, layer = 'counts') & xlab("")
VlnPlot(test, features="nCount_Akoya", layer='counts')
test[[]]
Layers(test)


Sample1_Group1_SeuratObj_sketch <- readRDS("./test_yj/Sample1_Group1_SeuratObj_sketch.rds")
VlnPlot(Sample1_Group1_SeuratObj_sketch, features="nCount_Akoya", layer='counts')
