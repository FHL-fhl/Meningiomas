#####12_ScRNA-seq quality control and integration######
library(scDblFinder)
library(Seurat)
library(dplyr)
library(cowplot)
library(sctransform)
library(SoupX)
library(glmGamPoi)
library(tidyverse)
library(ggalluvial)
library(seuratObject)
library(harmony)
options(future.globals.maxSize = 10 * 1024^3)  # 设置为 10 GB

# 设置数据目录路径
data_dir <- "/xxx/data"####下载单细胞数据放置在data文件夹中

# 获取数据目录下所有样本文件夹的名称
sample_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)

# 初始化一个空的列表来存储所有scobj对象
scobj_list <- list()


# 循环遍历每个样本文件夹
for (sample_dir in sample_dirs) {
  # 提取样本名称
  sample_name <- basename(sample_dir)
  
  # 读取数据
  sample_data <- Read10X(data.dir = sample_dir)
  
  # 创建scobj对象
  scobj_obj <- CreatescobjObject(counts = sample_data, min.cells = 10, min.features = 200, project = sample_name)
  
  # 识别双重细胞
  dbl <- scDblFinder(scobj_obj@assays$RNA@counts)
  
  # 添加双重细胞逻辑列
  scobj_obj$doublet_logic <- ifelse(dbl$scDblFinder.class == "doublet", TRUE, FALSE)
  
  # 过滤掉双重细胞
  scobj_obj <- subset(scobj_obj, subset = doublet_logic == FALSE)
  
  # 将处理后的scobj对象以样本名称为键添加到列表中
  scobj_list[[sample_name]] <- scobj_obj
}

### 数据质控
###合并多个样本

scobj_list[["S3"]]$gender<-"F"
scobj_list[["S3"]]$pathology<-"Meningothelial"
scobj_list[["S3"]]$WHO<-"WHO1"
scobj_list[["S3"]]$recurrence<-"Primary"
scobj_list[["S3"]]$group<-"Meningioma"

scobj_list[["S4"]]$gender<-"F"
scobj_list[["S4"]]$pathology<-"Transitional"
scobj_list[["S4"]]$WHO<-"WHO1"
scobj_list[["S4"]]$recurrence<-"Primary"
scobj_list[["S4"]]$group<-"Meningioma"

scobj_list[["S5"]]$gender<-"M"
scobj_list[["S5"]]$pathology<-"Meningothelial"
scobj_list[["S5"]]$WHO<-"WHO1"
scobj_list[["S5"]]$recurrence<-"Primary"
scobj_list[["S5"]]$group<-"Meningioma"

scobj_list[["S6"]]$gender<-"F"
scobj_list[["S6"]]$pathology<-"Transitional"
scobj_list[["S6"]]$WHO<-"WHO1"
scobj_list[["S6"]]$recurrence<-"Primary"
scobj_list[["S6"]]$group<-"Meningioma"

scobj_list[["S7"]]$gender<-"F"
scobj_list[["S7"]]$pathology<-"Meningothelial"
scobj_list[["S7"]]$WHO<-"WHO1"
scobj_list[["S7"]]$recurrence<-"Primary"
scobj_list[["S7"]]$group<-"Meningioma"

scobj_list[["S8"]]$gender<-"M"
scobj_list[["S8"]]$pathology<-"Transitional"
scobj_list[["S8"]]$WHO<-"WHO1"
scobj_list[["S8"]]$recurrence<-"Primary"
scobj_list[["S8"]]$group<-"Meningioma"

scobj_list[["S9"]]$gender<-"F"
scobj_list[["S9"]]$pathology<-"Meningothelial"
scobj_list[["S9"]]$WHO<-"WHO1"
scobj_list[["S9"]]$recurrence<-"Primary"
scobj_list[["S9"]]$group<-"Meningioma"

scobj_list[["S10"]]$gender<-"M"
scobj_list[["S10"]]$pathology<-"Meningothelial"
scobj_list[["S10"]]$WHO<-"WHO1"
scobj_list[["S10"]]$recurrence<-"Primary"
scobj_list[["S10"]]$group<-"Meningioma"

scobj_list[["S11"]]$gender<-"F"
scobj_list[["S11"]]$pathology<-"Atypical"
scobj_list[["S11"]]$WHO<-"WHO2"
scobj_list[["S11"]]$recurrence<-"Recurrent"
scobj_list[["S11"]]$group<-"Meningioma"

scobj_list[["S12"]]$gender<-"M"
scobj_list[["S12"]]$pathology<-"Atypical"
scobj_list[["S12"]]$WHO<-"WHO2"
scobj_list[["S12"]]$recurrence<-"Recurrent"
scobj_list[["S12"]]$group<-"Meningioma"

scobj_list[["S15"]]$gender<-"F"
scobj_list[["S15"]]$pathology<-"Anaplastic"
scobj_list[["S15"]]$WHO<-"WHO3"
scobj_list[["S15"]]$recurrence<-"Recurrent"
scobj_list[["S15"]]$group<-"Meningioma"

scobj_list[["S16"]]$gender<-"M"
scobj_list[["S16"]]$pathology<-"Anaplastic"
scobj_list[["S16"]]$WHO<-"WHO3"
scobj_list[["S16"]]$recurrence<-"Recurrent"
scobj_list[["S16"]]$group<-"Meningioma"


save(scobj_list,file="scobj_none_zhuanyi_zhengchang.Rdata")

samples <- basename(sample_dirs)
scobj<- merge(scobj_list[[1]], 
              scobj_list[2:length(scobj_list)],
              add.cell.ids = samples)

#将它们的counts并没有融合成为一个矩阵
#scobj <- JoinLayers(scobj)

scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.ribo"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")

### 正式筛选，筛选的是细胞，最终细胞减少
# 根据具体情况进行筛选阈值的确定
scobj <- subset(scobj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & nCount_RNA > 1000)
scobj <- NormalizeData(scobj)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)
scobj <- ScaleData(scobj, features = rownames(scobj))
scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj),reduction.name = "pca")
ElbowPlot(scobj)

DimPlot(scobj, reduction = "pca")

## 查看数据的复杂程度
ElbowPlot(scobj, reduction = "pca", ndims = 30)
xx <- cumsum(scobj[["pca"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.9) # 可以查看多少PCs解释了90%的方差，假设10%的方差来自于噪声，然后就可以选择相应的PCs
ndim = 30# dims的参数


scobj <- RunHarmony(scobj,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
scobj <- RunUMAP(scobj, reduction = "harmony", dims = 1:30,reduction.name = "harmony_umap")
scobj <- RunTSNE(scobj, reduction = "harmony", dims = 1:30,reduction.name = "harmony_tsne")
DimPlot(scobj, reduction = "harmony_umap",group.by = "orig.ident", label = T)
DimPlot(scobj, reduction = "harmony_tsne",group.by = "orig.ident", label = T)
scobj <- FindNeighbors(scobj, reduction = "harmony", dims = 1:30, verbose = FALSE)



#scobj <- FindClusters(scobj, resolution = 1,verbose=FALSE)
scobj <- FindClusters(scobj, resolution = 0.5,verbose=FALSE)
library(ggsci)
colors <- c(
  "#FF6F61", "#FFD662", "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1",
  "#955251", "#B565A7", "#009B77", "#DD4124", "#45B8AC", "#EFC050",
  "#5B5EA6", "#9B2335", "#DFCFBE", "#BC243C", "#E15D44", "#55B4B0",
  "#98B4D4", "#F7786B", "#DEAB8A", "#BFD641", "#D65076", "#86A59C"
)


DimPlot(scobj, reduction = "harmony_umap", label = TRUE, cols = colors)

scobj.markers <- FindAllMarkers(scobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
markers <-scobj.markers%>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(scobj, features = markers$gene) + NoLegend()


DimPlot(scobj, reduction = "harmony_umap", label = TRUE, cols = colors)
DimPlot(scobj, reduction = "harmony_umap",label = T)
DimPlot(scobj, reduction = "harmony_umap",split.by = "orig.ident", label = T,raster=FALSE)
DimPlot(scobj, reduction = "harmony_umap", split.by = "gender",label = TRUE)
DimPlot(scobj, reduction = "harmony_umap", split.by = "pathology",label = TRUE)
DimPlot(scobj, reduction = "harmony_umap", split.by = "WHO",label = TRUE)
DimPlot(scobj, reduction = "harmony_umap", split.by = "recurrence",label = TRUE)
DimPlot(scobj, reduction = "harmony_umap", split.by = "group",label = TRUE)


###定义细胞类型
new.cluster.ids <- c("Neoplastic_cell","Macrophages","T_cell","Macrophages","Neoplastic_cell" ,
                     "Neoplastic_cell","Endothelial","T_cell","Macrophages","Neoplastic_cell",
                     "Neutrophils","Dendritic Cells","Neoplastic_cell","Neoplastic_cell","Fibroblast",
                     "Neoplastic_cell","Macrophages","Macrophages","B_cell","Macrophages", 
                     "Macrophages","Oligodendrocyte","Macrophages","Neoplastic_cell")


names(new.cluster.ids) <- levels(scobj)              
new.cluster.ids
scobj <- RenameIdents(scobj, new.cluster.ids)
DimPlot(scobj, reduction = "harmony_umap", label = TRUE, pt.size = 0.5) 
save(scobj,file = "/xxx/scobj_none_zhuanyi_zhengchang.Rdata")


##细胞marker
DefaultAssay(scobj) <- "RNA"

DotPlot(scobj, features = c("CLU", "PTN", "LEPR",
                            "CD14","CD163", "CD68","HLA-DRB5", "MS4A6A","FCGR3A", "C1QA",
                            "CD3D", "CD3E", "CD3G", "CD8A", "CD8B",
                            "CD34", "VWF", "CCL14", "PLVAP",
                            "S100A8","S100A9","G0S2","CSF3R",
                            "CD1C","CLEC10A","FCER1A",
                            "ACTA2", "RGS5", 'CD19', 'MS4A1', 'CD79A',
                            "CNP", "MAG", "KLK6", "OLIG2")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_color_gradient(low = "#ADD8E6", high = "#FF4500") +  # 渐变从淡蓝到深红
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))


VlnPlot(scobj,features = c("CLU","CD14","CD3E","PLVAP","CSF3R","CD1C","RGS5","CD79A","CNP"),
        pt.size=0,
        cols = c("#e31a1c","#1f78b4", "#33a02c",  "#ff7f00", "#6a3d9a", "#b15928","#D8709D","#98B4D4","#DFCFBE"  ) )

library(patchwork)

# 生成单独的 FeaturePlot 子图
plots <- FeaturePlot(scobj, 
                     features = c("CLU","CD14","CD3E","PLVAP","RGS5","CD79A","CNP"), 
                     reduction = "harmony_umap", 
                     combine = FALSE)

# 为每个子图应用统一的渐变配色
plots <- lapply(plots, function(p) {
  p + scale_color_gradient(low = "#ADD8E6", high = "#FF4500")
})

combined_plot <- wrap_plots(plots, ncol = 3) 
print(combined_plot)



###计算各组样本不同细胞群比例
{
  table(Idents(scobj))
  Cellratio <- prop.table(table(Idents(scobj)))
  table(Idents(scobj), scobj$orig.ident)
  Cellratio <- prop.table(table(Idents(scobj), scobj$orig.ident), margin = 2)
  Cellratio
  Cellratio <- as.data.frame(Cellratio)
  library(ggplot2)
  # 定义完整的样本顺序 S1-S16（即使数据中部分缺失）
  full_sample_order <- c(paste0("S", 3:12),"S15","S16")
  
  # 将样本列转换为有序因子，强制包含所有 S1-S16 的标签
  Cellratio$Var2 <- factor(
    Cellratio$Var2,
    levels = full_sample_order,  # 强制顺序为 S1-S16
    ordered = TRUE
  )
  
  # 绘图代码
  ggplot(Cellratio) + 
    geom_bar(
      aes(x = Var2, y = Freq, fill = Var1),
      stat = "identity",
      width = 0.7
    ) +
    theme_classic() +
    coord_flip() +  # 翻转坐标系后，x轴变为纵轴
    labs(x = "Sample", y = "Ratio") +
    scale_fill_manual(
      values = c("#e31a1c","#1f78b4", "#33a02c",  "#ff7f00", "#6a3d9a", "#b15928","#D8709D","#98B4D4","#DFCFBE")
    ) +
    scale_x_discrete(
      limits = rev(full_sample_order),  # 强制显示 S1-S16，S1在顶部，S16在底部
      drop = FALSE  # 即使数据中缺失S1、S2、S13、S14，仍显示所有标签
    )
}




######macorphages re-clusters - no harmony again#####

load("/xxx/scobj_none_zhuanyi_zhengchang_cell.Rdata")

mar_scobj= scobj[, Idents(scobj) %in% c("Macrophages")]
mar_scobj <- NormalizeData(mar_scobj, normalization.method = "LogNormalize", scale.factor = 1e4) 
mar_scobj <- FindVariableFeatures(mar_scobj, selection.method = 'vst', nfeatures = 2000)
mar_scobj <- ScaleData(mar_scobj, vars.to.regress = "percent.mt")
mar_scobj <- RunPCA(mar_scobj, features = VariableFeatures(object = mar_scobj)) 


mar_scobj <- FindNeighbors(mar_scobj, dims = 1:30)
mar_scobj <- FindClusters(mar_scobj, resolution = 0.1 )
head(Idents(mar_scobj), 5)
table(mar_scobj$seurat_clusters)
mar_scobj <- RunUMAP(mar_scobj,dims = 1:30,reduction.name = "mar_umap")

mycol22 <- c("#b0d45d","#4c9568","#5066a1",
             "#76afda","#abddff","#fddbc8",
             "#ffc556", "#e8743c","#b20000",
             "#a14462","#cca69c","#7d4444")
p2<-DimPlot(mar_scobj, reduction = "mar_umap",group.by = "orig.ident", label = T, cols = mycol22)
p3<-DimPlot(mar_scobj, reduction = "mar_umap",group.by = "WHO", label = T)
p4<-DimPlot(mar_scobj, reduction = "mar_umap",group.by = "pathology", label = T)
p1<-DimPlot(mar_scobj, reduction = "mar_umap", label = TRUE, cols = mycol22)
FeaturePlot(mar_scobj,reduction ="mar_umap",features = c("CD68","CD80","CD86","TLR2"
                                                         ,"TLR4","CXCL2","IL1A","IL1B"))
FeaturePlot(mar_scobj,reduction ="mar_umap",features = c("CD163","MRC1","STAT6","PPARG"
                                                         ,"IL10","CLEC10A","CLEC7A","IRF4"))
FeaturePlot(mar_scobj,reduction ="mar_umap",features = c("SPP1"))

library(ggplot2)
library(patchwork)
(p1| p2) /(p3 | p4)
save(mar_scobj,file="mar_scobj_noharmony.Rdata")




######macorphages re-clusters - harmony again#####
load("/xxx/scobj_none_zhuanyi_zhengchang_cell.Rdata")
mar_scobj= scobj[, Idents(scobj) %in% c("Macrophages")]


mar_scobj=mar_scobj %>% RunHarmony("orig.ident", plot_convergence = TRUE, theta = 1,    # 原默认2，降低样本特异性惩罚
                                   lambda = 0.5,   # 原默认1，减少矫正强度
                                   reduction.save = "mar_harmony")
mar_scobj <- mar_scobj %>%
  RunUMAP(reduction = "mar_harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "mar_harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.1) %>%
  identity()

mycol22 <- c("#b0d45d","#4c9568","#5066a1",
             "#76afda","#abddff","#fddbc8",
             "#ffc556", "#e8743c","#b20000",
             "#a14462","#cca69c","#7d4444")
p2<-DimPlot(mar_scobj, reduction = "mar_harmony",group.by = "orig.ident", label = T, cols = mycol22)
p3<-DimPlot(mar_scobj, reduction = "mar_harmony",group.by = "WHO", label = T)
p4<-DimPlot(mar_scobj, reduction = "mar_harmony",group.by = "pathology", label = T)
p1<-DimPlot(mar_scobj, reduction = "mar_harmony", label = TRUE, cols = mycol22)
(p1| p2) /(p3 | p4)

save(mar_scobj,file="mar_scobj_harmony.Rdata")

