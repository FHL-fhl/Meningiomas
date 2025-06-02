#####14_Identification of CS1/CS4 macrophages######

library(Seurat)
library(ggplot2)
library(Scissor)
library(magrittr)
library(readr)


####CS4 vs Others####
load("/xxx/GSE212666_tpm.rda")
GSE212666_302sample <- read_csv("/xxx/GSE212666_302sample.csv")
a<-GSE212666_302sample[,c(10,12)]
a$Group <- ifelse(a$Cluster == "CS4", 1, 0)
df<-a[c(1,3)]

df$Group <- as.numeric(as.character(df$Group))
rownames(df) <- df$samID
phenotype <- as.numeric(as.character(df$Group))
names(phenotype) <- df$samID

load("/xxx/mar_scobj_noharmony.Rdata")


trace('Scissor', edit = T, where = asNamespace("Scissor"))
# 源代码 33 行
#dataset1 <- normalize.quantiles(dataset0)
# 修改为
#dataset1 <- normalize.quantiles(as.matrix(dataset0))
tag = c("Others","CS4")
result <- Scissor:::Scissor(bulk_dataset = GSE302_tpm, sc_dataset = mar_scobj,
                            phenotype = phenotype, alpha =0.5,, tag = tag,
                            family = "binomial", Save_file ="Scissor_CS4_inputs.rdata")
#可调整alpha值
result2 <- Scissor(GSE302_tpm, mar_scobj, phenotype, tag = tag, alpha = 0.2,
                   family = "binomial", Load_file = "Scissor_CS4_inputs.rdata")


save(result2,file="CS4_result2_alpha = 0.2.Rdata")
#可视化
Scissor_select <- rep(0, ncol(mar_scobj))
names(Scissor_select) <- colnames(mar_scobj)
Scissor_select[result2$Scissor_pos] <- 1
Scissor_select[result2$Scissor_neg] <- 2
mar_scobj <- AddMetaData(mar_scobj, metadata = Scissor_select, col.name = "scissor2")
DimPlot(mar_scobj, reduction = 'mar_umap', group.by = 'scissor2', cols = c('grey','indianred1','royalblue'), pt.size = 0.5, order = c(2,1))

###计算各组样本不同细胞群比例
mar_scobj_filtered <- mar_scobj[, mar_scobj$scissor2 != 0]
table(mar_scobj_filtered$orig.ident)
Cellratio <- prop.table(table(mar_scobj_filtered$scissor2))
table(mar_scobj_filtered$scissor2, mar_scobj_filtered$orig.ident)
Cellratio <- prop.table(table(mar_scobj_filtered$scissor2, mar_scobj_filtered$orig.ident), margin = 2)
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
  theme_classic() + # 翻转坐标系后，x轴变为纵轴
  labs(x = "Sample", y = "Ratio") +
  scale_fill_manual(
    values = c("#e31a1c","#1f78b4")
  ) +
  scale_x_discrete(
    limits = rev(full_sample_order),  # 强制显示 S1-S16，S1在顶部，S16在底部
    drop = FALSE  # 即使数据中缺失S1、S2、S13、S14，仍显示所有标签
  )



### 使用 FindMarkers() 函数计算两群之间的差异表达基因
DEG_results <- FindMarkers(mar_scobj_filtered, 
                           ident.1 = 1, 
                           ident.2 = 2, 
                           group.by = "scissor2")

DEG_results_filtered <- DEG_results[DEG_results$avg_log2FC > 0.5 & DEG_results$p_val_adj < 0.05, ]
head(DEG_results_filtered)
write.csv(DEG_results,file="CS4_relatived_DEG.csv")


####CS1 vs Others####
load("/xxx/GSE212666_tpm.rda")
GSE212666_302sample <- read_csv("/xxx/GSE212666_302sample.csv")
a<-GSE212666_302sample[,c(10,12)]
a$Group <- ifelse(a$Cluster == "CS1", 1, 0)
df<-a[c(1,3)]

df$Group <- as.numeric(as.character(df$Group))
rownames(df) <- df$samID
phenotype <- as.numeric(as.character(df$Group))
names(phenotype) <- df$samID

load("/xxx/mar_scobj_noharmony.Rdata")

trace(Scissor:::Scissor, edit = T)
# 源代码 33 行
#dataset1 <- normalize.quantiles(dataset0)
# 修改为
#dataset1 <- normalize.quantiles(as.matrix(dataset0))
tag = c("Others","CS1")
result_CS1 <- Scissor:::Scissor(bulk_dataset = GSE302_tpm, sc_dataset = mar_scobj,
                                phenotype = phenotype, alpha =0.2, tag = tag,
                                family = "binomial", Save_file ="Scissor_CS1_inputs.rdata")

save(result_CS1,file="CS1_result2_alpha = 0.2.Rdata")
#可视化
Scissor_select <- rep(0, ncol(mar_scobj))
names(Scissor_select) <- colnames(mar_scobj)
Scissor_select[result_CS1$Scissor_pos] <- 1
Scissor_select[result_CS1$Scissor_neg] <- 2
mar_scobj <- AddMetaData(mar_scobj, metadata = Scissor_select, col.name = "scissor2")
DimPlot(mar_scobj, reduction = 'mar_umap', group.by = 'scissor2', cols = c('grey','#b85ec2','#ffcd5e'), pt.size =0.5, order = c(2,1))


###计算各组样本不同细胞群比例
mar_scobj_filtered <- mar_scobj[, mar_scobj$scissor2 != 0]
table(mar_scobj_filtered$orig.ident)
Cellratio <- prop.table(table(mar_scobj_filtered$scissor2))
table(mar_scobj_filtered$scissor2, mar_scobj_filtered$orig.ident)
Cellratio <- prop.table(table(mar_scobj_filtered$scissor2, mar_scobj_filtered$orig.ident), margin = 2)
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
  theme_classic() + # 翻转坐标系后，x轴变为纵轴
  labs(x = "Sample", y = "Ratio") +
  scale_fill_manual(
    values = c('#b85ec2','#ffcd5e')
  ) +
  scale_x_discrete(
    limits = rev(full_sample_order),  # 强制显示 S1-S16，S1在顶部，S16在底部
    drop = FALSE  # 即使数据中缺失S1、S2、S13、S14，仍显示所有标签
  )



### 使用 FindMarkers() 函数计算两群之间的差异表达基因
DEG_results <- FindMarkers(mar_scobj_filtered, 
                           ident.1 = 1, 
                           ident.2 = 2, 
                           group.by = "scissor2")

DEG_results_filtered <- DEG_results[DEG_results$avg_log2FC > 0.5 & DEG_results$p_val_adj < 0.05, ]
head(DEG_results_filtered)
write.csv(DEG_results,file="CS1_relatived_DEG.csv")




####绘制象限图####
library(ggplot2)
library(ggrepel)
library(readxl)

CS1<- read_excel("/xxx/CS1_CS4.xlsx")
CS4<- read_excel("/xxx/CS1_CS4.xlsx",sheet="Sheet2")
CS1_CS4 <- merge(CS1,CS4,by="Gene")

# 替换为实际文件路径
data<-CS1_CS4
# 检查列名并修改为合适的格式
head(data)

# 添加象限信息
data$Quadrant <- with(data, ifelse(
  CS1_avg_log2FC > 0 & CS4_avg_log2FC > 0, "Q1",
  ifelse(CS1_avg_log2FC < 0 & CS4_avg_log2FC > 0, "Q2",
         ifelse(CS1_avg_log2FC < 0 & CS4_avg_log2FC < 0, "Q3", "Q4")
  )
))

# 定义颜色
quadrant_colors <- c("Q1" ="blue" , "Q2" = "red", "Q3" = "green", "Q4" = "purple")
# 绘制散点图
ggplot(data, aes(x = CS4_avg_log2FC, y = CS1_avg_log2FC, color = Quadrant, label = Gene)) +
  geom_point(size = 3) +  # 绘制点
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # 添加y=0线
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # 添加x=0线
  geom_point(aes(x = 0, y = 0), color = "black", size = 5, shape = 21, fill = "yellow") +  # 加粗显示原点
  geom_text_repel(size = 3, max.overlaps = 20) +  # 避免基因名称重叠
  scale_color_manual(values = quadrant_colors) +  # 设置象限颜色
  labs(
    title = "Scatter Plot with Quadrant Colors",
    x = "CS1_avg_log2FC",
    y = "CS4_avg_log2FC",
    color = "Quadrant"
  ) +
  theme_minimal()

