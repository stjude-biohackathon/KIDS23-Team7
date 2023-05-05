### Statistics - Cell Proportion Between Samples ###

###** GROUP 1 **###

# Load metadata
library(readr)
Group1_df_meta <- read_csv("/home/lead/Akoya/Results/Group1_df_meta.csv")
Group1_df_meta$Sample # it seems columns have been shifted 
Group1_df_meta = data.frame(Group1_df_meta)
head(Group1_df_meta)
dim(Group1_df_meta) #2184767      19

# Fix column names
Group1_df_meta_sub = Group1_df_meta[,c("orig.ident",
                                       "leverage.score.1", 
                                       "cluster_full.score")]
names(Group1_df_meta_sub) = c("orig.ident", "sample", "cluster_and_score")

# Create separated columns for celltype and score
Group1_df_meta_sub = Group1_df_meta_sub %>%
  tidyr::separate(cluster_and_score, into = c("celltype_cat", "celltype_cont"), sep = ",")
names(Group1_df_meta_sub)
head(Group1_df_meta_sub)
dim(Group1_df_meta_sub) #2184767       4

# Subset by samples
table(Group1_df_meta_sub$sample)
Smp_1 = Group1_df_meta_sub[Group1_df_meta_sub$sample %in% "Sample1_Group1", ]
Smp_2 = Group1_df_meta_sub[Group1_df_meta_sub$sample %in% "Sample2_Group1", ]

# Cell proportion - Smp1
Smp_1_cellfreq = data.frame(table(Smp_1$celltype_cat))
colnames(Smp_1_cellfreq) = c("celltype_cat", "ncells_celltype")
Smp_1_cellfreq$sample = "Sample1_Group1"
Smp_1_cellfreq$propCells <- Smp_1_cellfreq$ncells_celltype / sum(Smp_1_cellfreq$ncells_celltype)
sum(Smp_1_cellfreq$propCells) == 1 # TRUE

# Cell proportion - Smp2
Smp_2_cellfreq = data.frame(table(Smp_2$celltype_cat))
colnames(Smp_2_cellfreq) = c("celltype_cat", "ncells_celltype")
Smp_2_cellfreq$sample = "Sample2_Group1"
Smp_2_cellfreq$propCells <- Smp_2_cellfreq$ncells_celltype / sum(Smp_2_cellfreq$ncells_celltype)
sum(Smp_2_cellfreq$propCells) == 1 # TRUE

# bind both dataframes
propCells_df = rbind(Smp_1_cellfreq, Smp_2_cellfreq)
head(propCells_df, 50)

# Plot stack barplot of cell proportions by sample
library(ggplot2); theme_set(theme_classic())
propCells_df = propCells_df[!propCells_df$celltype_cat %in% "B cells", ] # to match with group2
ggplot(propCells_df, aes(x=sample, y=propCells, fill=celltype_cat)) +
  scale_fill_manual(values = jet.colors(length(unique(propCells_df$celltype_cat))), name="Cell types") +
  geom_bar(stat="identity") +
  labs(x="Sample", y="Cell Proportion") +
  theme(axis.title.x = element_text(colour="black", size = 12, vjust=-1), 
        axis.title.y = element_text(colour="black", size = 12), 
        axis.text.x = element_text(size = 12)) + 
  ggtitle("Group1") 


ggsave(file="/home/mmarcao/statistics/Group1_barplot_celltypeProportion.pdf",width = 10, height = 6, dpi = 300)

save(propCells_df, file = "/home/mmarcao/statistics/Group1_obj_to_barplot_stats.rda")
load("/home/mmarcao/statistics/Group1_obj_to_barplot_stats.rda")

# Statistics
celltypes_vec = as.vector(unique(factor(Smp_1$celltype_cat)))
Smp_1$celltype_cont = as.numeric(Smp_1$celltype_cont)
Smp_2$celltype_cont = as.numeric(Smp_2$celltype_cont)

list_results = list()

for(i in 1:length(unique(propCells_df$celltype_cat))) {
  result = t.test(Smp_1[Smp_1$celltype_cat %in% celltypes_vec[i], ]$celltype_cont,
                  Smp_2[Smp_2$celltype_cat %in% celltypes_vec[i], ]$celltype_cont,
                  conf.level = 0.95,
                  alternative = c("two.sided"),
                  paired = FALSE,
                  var.equal = FALSE) 
  
  conf_int = paste0("[", as.character(result$conf.int[1]), "  ", as.character(result$conf.int[2]),"]")
  pvalue = result$p.value 
  
  pvalue_axis = data.frame(comparison = celltypes_vec[i],
                           pvalue = pvalue,
                           conf_int = conf_int)
  
  list_results[[i]] = pvalue_axis
  
}

pvalues_df = do.call('rbind', list_results) 

# T.TEST PROBABLY ISN'T THE BEST WAY. EACH CELLTYPE SCORE HAS A DIFFERENT DISTRIBUTIONS... 
# TAKE A LOOK
i = 2
par(mfrow = c(1,2))
Smp_1[Smp_1$celltype_cat %in% celltypes_vec[i], ]$celltype_cont %>% hist()
Smp_2[Smp_2$celltype_cat %in% celltypes_vec[i], ]$celltype_cont %>% hist()


# Instead, we could ttest directily on the numbers used to calculate FoldChange - aka "cell percentage" itself
i = 5
Smp_1_cellfreq[Smp_1_cellfreq$celltype_cat %in% celltypes_vec[i], ]$propCells
Smp_2_cellfreq[Smp_2_cellfreq$celltype_cat %in% celltypes_vec[i], ]$propCells

# Even though we still need to have more PATIENTS to conduce this kind of statistics. We just have two patients per group. 



# So, we should plot only Percentage (x axis) by median of sc-type score (y axis)

# Sample 1
celltypes_vec = as.vector(unique(factor(Smp_1$celltype_cat)))
df = data.frame(celltype_cat = "test", median_score = NA)
for(i in 1:length(celltypes_vec)) {
  df[i, 1] = celltypes_vec[i]
  df[i, 2] = Smp_1[Smp_1$celltype_cat %in% celltypes_vec[i], ]$celltype_cont %>% median()
  
}

Smp_1_cellfreq = merge(Smp_1_cellfreq, df)


# Sample 2
celltypes_vec = as.vector(unique(factor(Smp_1$celltype_cat)))
df = data.frame(celltype_cat = "test", median_score = NA)
for(i in 1:length(celltypes_vec)) {
  df[i, 1] = celltypes_vec[i]
  df[i, 2] = Smp_2[Smp_2$celltype_cat %in% celltypes_vec[i], ]$celltype_cont %>% median()
  
}

Smp_2_cellfreq = merge(Smp_2_cellfreq, df)

# Scatter plot 

to_plot = rbind(Smp_1_cellfreq, Smp_2_cellfreq)

library(ggplot2); theme_set(theme_classic())
ggplot(to_plot, aes(x = propCells, y = median_score, color = sample, label = celltype_cat, size = ncells_celltype)) + 
  geom_point(alpha = 0.7) +
  #scale_fill_manual(values = jet.colors(2)) +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  scale_size(range = c(1,10)) +
  labs(x = "Celltype Percentage", y = "Celltype Median Score", color = "Sample", size = "Celltype nCells") +
  geom_label_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 100)) +
  ggtitle("Group1")

ggsave(file="/home/mmarcao/statistics/Group1_scatterplot_celltypeProportion_perMedianScore.pdf",width = 10, height = 4, dpi = 300)

save(to_plot, file = "/home/mmarcao/statistics/Group1_obj_to_scatterplot_stats.rda")
load("/home/mmarcao/statistics/Group1_obj_to_scatterplot_stats.rda")


# DON'T RUN THIS CHUNK OF CODE 
# # Boxplot 
# celltypes_vec = as.vector(unique(factor(Smp_1$celltype_cat)))
# df = data.frame(celltype_cat = "test", FoldChange = NA, Orientation = "Smp1/Smp2")
# for(i in 1:length(celltypes_vec)) {
#   
#   propCell_Smp1 = to_plot[to_plot$sample %in% "Sample1_Group1" &
#                             to_plot$celltype_cat %in% celltypes_vec[i], ]$propCells
#   
#   
#   propCell_Smp2 = to_plot[to_plot$sample %in% "Sample2_Group1" &
#                             to_plot$celltype_cat %in% celltypes_vec[i], ]$propCells
#   
#   df[i, 1] = celltypes_vec[i]
#   df[i, 2] = as.numeric(propCell_Smp1/propCell_Smp2)
#   df[i, 3] = "Smp1/Smp2"
#   
# }
# to_plot_2 = df
# 
# library(ggplot2); theme_set(theme_classic())
# ggplot(to_plot_2[,], aes(x=celltype_cat, y=FoldChange), fill = stemness) + 
#   geom_boxplot(outlier.color = NA)+ geom_jitter (alpha=0.01)  +
#   xlab("Cell type") + 
#   ylab("FoldChange - Cell Percentage") + 
#   theme(axis.title.x = element_text(colour="black", size = 12, vjust=-1), 
#         axis.title.y = element_text(colour="black", size = 12), 
#         axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 10)) + 
#   ggtitle("Group 1") 
# 
# ggsave(file="/home/mmarcao/statistics/Group1_boxplot_FoldChange_celltypeProportion.pdf",width = 6, height = 4, dpi = 300)
# 
# save(to_plot_2, file = "/home/mmarcao/statistics/Group1_obj_to_boxplot_stats.rda")
# load("/home/mmarcao/statistics/Group1_obj_to_boxplot_stats.rda")


###** GROUP 2 **###

# Barplot

# Load metadata
library(readr)
Group2_df_meta <- read_csv("/home/lead/Akoya/Results/Group2_df_meta.csv")
Group2_df_meta$Sample # it seems columns have been shifted 
Group2_df_meta = data.frame(Group2_df_meta)
head(Group2_df_meta)
dim(Group2_df_meta) #416587     18

# Fix column names
Group2_df_meta_sub = Group2_df_meta[,c("orig.ident",
                                       "leverage.score.1", 
                                       "cluster_full.score")]
names(Group2_df_meta_sub) = c("orig.ident", "sample", "cluster_and_score")

# Create separated columns for celltype and score
Group2_df_meta_sub = Group2_df_meta_sub %>%
  tidyr::separate(cluster_and_score, into = c("celltype_cat", "celltype_cont"), sep = ",")
names(Group2_df_meta_sub)
head(Group2_df_meta_sub)
dim(Group2_df_meta_sub) #2184767       4

# Subset by samples
table(Group2_df_meta_sub$sample)
Smp_1 = Group2_df_meta_sub[Group2_df_meta_sub$sample %in% "Sample1_Group2", ]
Smp_2 = Group2_df_meta_sub[Group2_df_meta_sub$sample %in% "Sample2_Group2", ]

# Cell proportion - Smp1
Smp_1_cellfreq = data.frame(table(Smp_1$celltype_cat))
colnames(Smp_1_cellfreq) = c("celltype_cat", "ncells_celltype")
Smp_1_cellfreq$sample = "Sample1_Group2"
Smp_1_cellfreq$propCells <- Smp_1_cellfreq$ncells_celltype / sum(Smp_1_cellfreq$ncells_celltype)
sum(Smp_1_cellfreq$propCells) == 1 # TRUE

# Cell proportion - Smp2
Smp_2_cellfreq = data.frame(table(Smp_2$celltype_cat))
colnames(Smp_2_cellfreq) = c("celltype_cat", "ncells_celltype")
Smp_2_cellfreq$sample = "Sample2_Group2"
Smp_2_cellfreq$propCells <- Smp_2_cellfreq$ncells_celltype / sum(Smp_2_cellfreq$ncells_celltype)
sum(Smp_2_cellfreq$propCells) == 1 # TRUE

# bind both dataframes
propCells_df = rbind(Smp_1_cellfreq, Smp_2_cellfreq)
head(propCells_df, 50)

# Plot stack barplot of cell proportions by sample
library(ggplot2); theme_set(theme_classic())
ggplot(propCells_df, aes(x=sample, y=propCells, fill=celltype_cat)) +
  scale_fill_manual(values = jet.colors(length(unique(propCells_df$celltype_cat))), name="Cell types") +
  geom_bar(stat="identity") +
  labs(x="Sample", y="Cell Proportion") +
  theme(axis.title.x = element_text(colour="black", size = 12, vjust=-1), 
        axis.title.y = element_text(colour="black", size = 12), 
        axis.text.x = element_text(size = 12)) + 
  ggtitle("Group2") 


ggsave(file="/home/mmarcao/statistics/Group2_barplot_celltypeProportion.pdf",width = 10, height = 6, dpi = 300)

save(propCells_df, file = "/home/mmarcao/statistics/Group2_obj_to_barplot_stats.rda")
load("/home/mmarcao/statistics/Group2_obj_to_barplot_stats.rda")



# Scatter plot
# Sample 1
celltypes_vec = as.vector(unique(factor(Smp_1$celltype_cat)))
Smp_1$celltype_cont = as.numeric(Smp_1$celltype_cont)
Smp_2$celltype_cont = as.numeric(Smp_2$celltype_cont)
df = data.frame(celltype_cat = "test", median_score = NA)
for(i in 1:length(celltypes_vec)) {
  df[i, 1] = celltypes_vec[i]
  df[i, 2] = Smp_1[Smp_1$celltype_cat %in% celltypes_vec[i], ]$celltype_cont %>% median()
  
}

Smp_1_cellfreq = merge(Smp_1_cellfreq, df)


# Sample 2
df = data.frame(celltype_cat = "test", median_score = NA)
for(i in 1:length(celltypes_vec)) {
  df[i, 1] = celltypes_vec[i]
  df[i, 2] = Smp_2[Smp_2$celltype_cat %in% celltypes_vec[i], ]$celltype_cont %>% median()
  
}

Smp_2_cellfreq = merge(Smp_2_cellfreq, df)

to_plot = rbind(Smp_1_cellfreq, Smp_2_cellfreq)

library(ggplot2); theme_set(theme_classic())
ggplot(to_plot, aes(x = propCells, y = median_score, color = sample, label = celltype_cat, size = ncells_celltype)) + 
  geom_point(alpha = 0.7) +
  #scale_fill_manual(values = jet.colors(2)) +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  scale_size(range = c(1,10)) +
  labs(x = "Celltype Percentage", y = "Celltype Median Score", color = "Sample", size = "Celltype nCells") +
  geom_label_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 100)) +
  ggtitle("Group 2")

ggsave(file="/home/mmarcao/statistics/Group2_scatterplot_celltypeProportion_perMedianScore.pdf",width = 10, height = 4, dpi = 300)

save(to_plot, file = "/home/mmarcao/statistics/Group2_obj_to_scatterplot_stats.rda")
load("/home/mmarcao/statistics/Group2_obj_to_scatterplot_stats.rda")


# DON'T RUN THIS CHUNK OF CODE 
# # Boxplot 
# celltypes_vec = as.vector(unique(factor(Smp_1$celltype_cat)))
# df = data.frame(celltype_cat = "test", FoldChange = NA, Orientation = "Smp1/Smp2")
# for(i in 1:length(celltypes_vec)) {
#   
#   propCell_Smp1 = to_plot[to_plot$sample %in% "Sample1_Group2" &
#                             to_plot$celltype_cat %in% celltypes_vec[i], ]$propCells
#   
#   
#   propCell_Smp2 = to_plot[to_plot$sample %in% "Sample2_Group2" &
#                             to_plot$celltype_cat %in% celltypes_vec[i], ]$propCells
#   
#   df[i, 1] = celltypes_vec[i]
#   df[i, 2] = as.numeric(propCell_Smp1/propCell_Smp2)
#   df[i, 3] = "Smp1/Smp2"
#   
# }
# 
# to_plot_2 = df
# 
# library(ggplot2); theme_set(theme_classic())
# ggplot(to_plot_2[,], aes(x=celltype_cat, y=FoldChange), fill = stemness) + 
#   geom_boxplot(outlier.color = NA)+ geom_jitter (alpha=0.01)  +
#   xlab("Cell type") + 
#   ylab("FoldChange - Cell Percentage") + 
#   theme(axis.title.x = element_text(colour="black", size = 12, vjust=-1), 
#         axis.title.y = element_text(colour="black", size = 12), 
#         axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 10)) + 
#   ggtitle("Group 2") 
# 
# ggsave(file="/home/mmarcao/statistics/Group2_boxplot_FoldChange_celltypeProportion.pdf",width = 6, height = 4, dpi = 300)
# 
# save(to_plot_2, file = "/home/mmarcao/statistics/Group2_obj_to_boxplot_stats.rda")
# load("/home/mmarcao/statistics/Group2_obj_to_boxplot_stats.rda") 




###** GROUP 1  vs  GROUP 2 **###

load("/home/mmarcao/statistics/Group1_obj_to_barplot_stats.rda")
G1_propCells_df = propCells_df

load("/home/mmarcao/statistics/Group2_obj_to_barplot_stats.rda")
G2_propCells_df = propCells_df

# Checking common cell types between groups
length(intersect(as.vector(unique(factor(G1_propCells_df$celltype_cat))),
                 as.vector(unique(factor(G2_propCells_df$celltype_cat))))) # 12 common celltypes

# Filtering G1_propCells_df
G1_propCells_df = G1_propCells_df[G1_propCells_df$celltype_cat %in% as.vector(unique(factor(G2_propCells_df$celltype_cat))), ]


celltypes_vec = as.vector(unique(factor(G1_propCells_df$celltype_cat)))
t.test(G1_propCells_df[G1_propCells_df$celltype_cat %in% celltypes_vec[1], ]$propCells,
       G2_propCells_df[G2_propCells_df$celltype_cat %in% celltypes_vec[1], ]$propCells)


# Pvalue
celltypes_vec = as.vector(unique(factor(G1_propCells_df$celltype_cat)))
list_results = list()

for(i in 1:length(celltypes_vec)) {
  result = t.test(G1_propCells_df[G1_propCells_df$celltype_cat %in% celltypes_vec[i], ]$propCells,
                  G2_propCells_df[G2_propCells_df$celltype_cat %in% celltypes_vec[i], ]$propCells,
                  conf.level = 0.95,
                  alternative = c("two.sided"),
                  paired = FALSE,
                  var.equal = FALSE) 
  
  conf_int = paste0("[", as.character(result$conf.int[1]), "  ", as.character(result$conf.int[2]),"]")
  pvalue = result$p.value 
  
  pvalue_axis = data.frame(comparison = celltypes_vec[i],
                           pvalue = pvalue,
                           conf_int = conf_int)
  
  list_results[[i]] = pvalue_axis
  
}

pvalues_df = do.call('rbind', list_results) 
hist(pvalues_df$pvalue, breaks = 10)


# FoldChange
FC_df = data.frame(comparison = "test", FoldChange = NA, Orientation = "G1/G2")
for(i in 1:length(celltypes_vec)) {
  
  G1 = mean(G1_propCells_df[G1_propCells_df$celltype_cat %in% celltypes_vec[i], ]$propCells)
  
  
  G2 = mean(G2_propCells_df[G2_propCells_df$celltype_cat %in% celltypes_vec[i], ]$propCells)
  
  FC_df[i, 1] = celltypes_vec[i]
  FC_df[i, 2] = as.numeric(G1/G2)
  FC_df[i, 3] = "G1/G2"
  
}


FC_df
pvalues_df$pvalue_log10 = -log10(pvalues_df$pvalue) 
plot_FC_pvalue = merge(pvalues_df, FC_df)

library(ggplot2); theme_set(theme_classic())
ggplot(plot_FC_pvalue, aes(x = FoldChange, y = pvalue_log10, color = comparison, label = comparison)) + 
  geom_point(alpha = 1) +
  #scale_fill_manual(values = jet.colors(2)) +
  scale_color_manual(values = jet.colors(length(as.vector(unique(factor(plot_FC_pvalue$comparison)))))) +
  #scale_size(range = c(1,10)) +
  labs(x = "FoldChange ", y = "-log10(pvalue)", color = "comparison") +
  geom_label_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 100)) +
  ggtitle("<---- Group2         Group1 ---->")

ggsave(file="/home/mmarcao/statistics/Group1_Group2_obj_to_scatterplot_FoldChange.pdf",width = 10, height = 6, dpi = 300)

save(plot_FC_pvalue, file = "/home/mmarcao/statistics/Group1_Group2_obj_to_scatterplot_FoldChange_stats.rda")



#barplot function interactive plot
#load df proportion of cells
load("/home/mmarcao/statistics/Group1_obj_to_barplot_stats.rda")
G1_propCells_df = propCells_df

load("/home/mmarcao/statistics/Group2_obj_to_barplot_stats.rda")
G2_propCells_df = propCells_df


barplot_proportion_interactive <- function(df, title, xname, yname){
  df$propCells <- (round(df$propCells, digits = 4))*100
  plot <- ggplot(df, aes(x=sample, y=propCells, fill=celltype_cat)) +
    scale_fill_manual(values = jet.colors(length(unique(df$celltype_cat))), name="Cell types") +
    geom_bar(stat="identity") +
    labs(x = xname, y = yname) +
    theme(axis.title.x = element_text(colour="black", size = 12, vjust=-1), 
          axis.title.y = element_text(colour="black", size = 12), 
          axis.text.x = element_text(size = 12)) + 
    ggtitle(title) +
    theme_classic()
  #plot <- ggplotly(plot, tooltip = c("sample", "propCells", "celltype_cat"))
  plot <- ggplotly(plot, tooltip = c("propCells", "celltype_cat"))
  return(plot)
}
#ploting for each group
barplot_proportion_interactive(G1_propCells_df, title = "Group1", xname = "Sample", yname = "Cell Proportion (%)")
barplot_proportion_interactive(G2_propCells_df, title = "Group2", xname = "Sample", yname = "Cell Proportion (%)")


#scatter plot FC pvalue interactive
# load plot_FC_pvalue df
load("/home/mmarcao/statistics/Group1_Group2_obj_to_scatterplot_FoldChange_stats.rda")

scatterplot_FC_pval_interactive <- function(plot_FC_pvalue){
 plot <- ggplot(plot_FC_pvalue, aes(x = FoldChange, y = pvalue_log10, color = comparison, label = comparison)) + 
    geom_point(alpha = 1) +
    #scale_fill_manual(values = jet.colors(2)) +
    scale_color_manual(values = jet.colors(length(as.vector(unique(factor(plot_FC_pvalue$comparison)))))) +
    #scale_size(range = c(1,10)) +
    labs(x = "FoldChange ", y = "-log10(pvalue)", color = "comparison") +
    geom_label_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 100)) +
    ggtitle("<---- Group2         Group1 ---->") +
    theme_classic() 
    plot <- ggplotly(plot)
    return(plot)
}

scatterplot_FC_pval_interactive(plot_FC_pvalue)
