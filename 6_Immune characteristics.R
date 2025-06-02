#####6_Immune characteristics######


#install.packages("remotes")
#remotes::install_github("omnideconv/immunedeconv")
library(immunedeconv)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(gplots)



load("/xxx/GSE212666_tpm.rda")
load("/xxx/GSE212_cluster5.rda")


an =as.data.frame(a[,c(10)])
rownames(an)<-rownames(a)
colnames(an)<-"cluster"


# 假设 an 是你的数据框
an <- an %>%
  mutate(cluster = recode(cluster, 
                          `1` = "CS1", 
                          `2` = "CS2", 
                          `3` = "CS3", 
                          `4` = "CS4", 
                          `5` = "CS5"))
an <- an %>%
  arrange(cluster)


dat<-GSE302_tpm


#利用estimate算法进行免疫浸润分析
estimate<-immunedeconv::deconvolute(dat, "estimate")
estimate<-as.data.frame(estimate)
estimate$cell_type<-paste0(estimate$cell_type,"_ESTIMATE")
estimate$Immunemethod<-"ESTIMATE"



#cibersort
load("/xxx/Cibersort_RNA_results.Rdata")
cibersort<-t(RNA_results)
cibersort<-as.data.frame(cibersort)
cibersort$cell_type<-rownames(cibersort)
cibersort <- cibersort[, c("cell_type", setdiff(names(cibersort), "cell_type"))]
rownames(cibersort)<-NULL
cibersort<-as.data.frame(cibersort)
cibersort$cell_type<-paste0(cibersort$cell_type,"_CIBESORT")
cibersort$Immunemethod<-"CIBERSORT"



#利用abis算法进行免疫浸润分析
abis<-immunedeconv::deconvolute(dat, "abis")
abis<-as.data.frame(abis)
abis$cell_type<-paste0(abis$cell_type,"_ABIS")
abis$Immunemethod<-"ABIS"


#利用quantiseq算法进行免疫浸润分析
quantiseq<-immunedeconv::deconvolute(dat, "quantiseq")
quantiseq<-as.data.frame(quantiseq)
quantiseq$cell_type<-paste0(quantiseq$cell_type,"_QUANTISEQ")
quantiseq$Immunemethod<-"QUANTISEQ"


#利用xcell算法进行免疫浸润分析
xcell<-immunedeconv::deconvolute(dat, "xcell")
xcell<-as.data.frame(xcell)





#利用epic算法进行免疫浸润分析
epic<-immunedeconv::deconvolute(dat, "epic")
epic<-as.data.frame(epic)
epic$cell_type<-paste0(epic$cell_type,"_EPIC")
epic$Immunemethod<-"EPIC"


#利用mcp_counter算法进行免疫浸润分析
mcp_counter<-immunedeconv::deconvolute(dat,"mcp_counter")
mcp_counter<-as.data.frame(mcp_counter)
mcp_counter$cell_type<-paste0(mcp_counter$cell_type,"_MCPCOUNTER")
mcp_counter$Immunemethod<-"MCPCOUNTER"



#合并免疫细胞矩阵
Immune<-rbind(xcell,cibersort,epic,mcp_counter,abis,quantiseq)
#提取免疫细胞矩阵
final<-select(Immune,-Immunemethod)
#把列名改成行名
final<-column_to_rownames(final,var="cell_type")
rownames(final)<-final$cell_type
final<-final[,-1]
#匹配高低风险组样本信息
Immunesample<-match(rownames(an),colnames(final))
final<-final[,Immunesample]
save(final,file ="final.Rdata")



####巨噬细胞###
load("/xxx/final.Rdata")
a<-as.data.frame(rownames(final))
cibersort<-final[c(8:14),]
cibersort$Immunemethod<-"CIBERSORT"

abis<-final[c(87,88,89,91,92,93,95),]
abis$Immunemethod<-"ABIS"

quantiseq<-final[c(108:110),]
quantiseq$Immunemethod<-"QUANTISEQ"

xcell<-final[c(29:37),]
xcell$Immunemethod<-"XCELL"

epic<-final[c(68:69),]
epic$Immunemethod<-"EPIC"

mcp_counter<-final[c(74:75),]
mcp_counter$Immunemethod<-"MCPCOUNTER"




####T细胞###
load("/xxx/final.Rdata")
a<-as.data.frame(rownames(final))
cibersort<-final[c(17:20),]
cibersort$Immunemethod<-"CIBERSORT"

abis<-final[c(85,97),]
abis$Immunemethod<-"ABIS"

quantiseq<-final[c(103:105),]
quantiseq$Immunemethod<-"QUANTISEQ"

xcell<-final[c(47:49,52),]
xcell$Immunemethod<-"XCELL"

epic<-final[c(71),]
epic$Immunemethod<-"EPIC"

mcp_counter<-final[c(79:80),]
mcp_counter$Immunemethod<-"MCPCOUNTER"



####B细胞###
load("/xxx/final.Rdata")
a<-as.data.frame(rownames(final))

cibersort<-final[c(5:6),]
cibersort$Immunemethod<-"CIBERSORT"

abis<-final[c(90,98,101),]
abis$Immunemethod<-"ABIS"

quantiseq<-final[c(102),]
quantiseq$Immunemethod<-"QUANTISEQ"

xcell<-final[c(28,51,53),]
xcell$Immunemethod<-"XCELL"

epic<-final[c(66),]
epic$Immunemethod<-"EPIC"

mcp_counter<-final[c(78),]
mcp_counter$Immunemethod<-"MCPCOUNTER"




# 定义颜色
#3.绘制免疫细胞热图
methods.col <- brewer.pal(n = length(unique(Immune$Immunemethod)),name = "Paired")
# 创建注释
# 列注释，位于热图顶端
annCol <- data.frame(CLUSTER = an$cluster,
                     # 以上是risk score和risk type两种注释，可以按照这样的格式继续添加更多种类的注释信息，记得在下面的annColors里设置颜色
                     row.names = row.names(an),
                     stringsAsFactors = F)

# 行注释，位于热图左侧
annRow<-data.frame(Methods=factor(Immune$Immunemethod,levels=unique(Immune$Immunemethod)),
                   row.names = Immune$cell_type,
                   stringsAsFactors = F)

# 为各注释信息设置颜色estimate,cibersort,xcell,epic,mcp_counter,abis,quantiseq
annColors <- list(Methods = c("ESTIMATE" = methods.col[1],
                              "XCELL" = methods.col[2],
                              "CIBERSORT"= methods.col[3],
                              "EPIC" = methods.col[4],
                              "MCPCOUNTER" = methods.col[5],
                              "ABIS" = methods.col[6],
                              "QUANTISEQ" = methods.col[7]
),
# 下面是列注释的颜色，可依此设置更多注释的颜色
"CLUSTER" = c("CS1" ="#2EC4B6","CS2" ="#E71D36",
              "CS3" ="#FF9F1C","CS4" ="#BDD5EA",
              "CS5" = "#FFA5AB"))

# 数据标准化
indata <- final # 确保没有富集全为0的细胞
#确定免疫细胞含量范围
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}
plotdata <- standarize.fun(indata,halfwidth = 2)
# pheatmap绘图
pheatmap::pheatmap(mat = as.matrix(plotdata), # 标准化后的数值矩阵
                   border_color = NA, # 无边框色
                   color = bluered(64), # 热图颜色为红蓝
                   cluster_rows = F, # 行不聚类
                   cluster_cols = F, # 列不聚类
                   show_rownames = T, # 显示行名
                   show_colnames = F, # 不显示列名
                   annotation_col = annCol, # 列注释
                   annotation_row = annRow, # 行注释
                   annotation_colors = annColors, # 注释颜色
                   gaps_col = c(table(annCol$CLUSTER)[1],
                                table(annCol$CLUSTER)[1]+table(annCol$CLUSTER)[2],
                                table(annCol$CLUSTER)[1]+table(annCol$CLUSTER)[2]+table(annCol$CLUSTER)[3],
                                table(annCol$CLUSTER)[1]+table(annCol$CLUSTER)[2]+table(annCol$CLUSTER)[3]+table(annCol$CLUSTER)[4],
                                table(annCol$CLUSTER)[1]+table(annCol$CLUSTER)[2]+table(annCol$CLUSTER)[3]+table(annCol$CLUSTER)[4]+table(annCol$CLUSTER)[5]), # 列分割
                   gaps_row = cumsum(table(annRow$Methods)), # 行分割
                   cellwidth = 0.8, # 元素宽度
                   cellheight = 10 ,# 元素高度
                   filename = "immune_heatmap_by_pheatmap_Tcell.pdf")#保存图片








comparisons <- list(
  c("CS1", "CS2"),
  c("CS1", "CS3"),
  c("CS1", "CS4"),
  c("CS1", "CS5"),
  c("CS2", "CS3"),
  c("CS2", "CS4"),
  c("CS2", "CS5"),
  c("CS3", "CS4"),
  c("CS3", "CS5"),
  c("CS4", "CS5")
)

# 使用 Set1 作为配色方案
library(RColorBrewer)
color_palette <- brewer.pal(n = 5, name = "Set1")
re<-t(indata)
b<-merge(an,re,by=0)
rownames(b)<-b$Row.names
b<-b[,-1]
# 假设你已经有数据框 b 和 color_palette，且 comparisons 也已定义

# 创建一个空列表来保存所有图
# 加载必要的包
library(ggplot2)
library(patchwork)
library(ggpubr)
library(rlang)  # 用于 tidy evaluation

# 创建一个空列表来保存所有图
plot_list <- list()

# 遍历列名，从第2列到第23列
for (i in 2:113) {
  col_name <- colnames(b)[i]
  
  # 使用 sym() 将列名转换为符号
  y_axis <- sym(col_name)
  
  # 创建单个 ggplot 图
  p <- ggplot(b, aes(x = cluster, y = !!y_axis, fill = cluster)) +
    geom_boxplot() +
    scale_fill_manual(values = color_palette) +
    stat_compare_means(comparisons = comparisons, 
                       method = "t.test", 
                       label = "p.signif") +
    theme_classic() +
    ggtitle(col_name)  # 添加图表标题，以便识别每个图表对应的列名
  
  # 将图表添加到列表中
  plot_list[[i - 1]] <- p
}


# Number of lists per plot
lists_per_plot <- 4

# Number of columns in each plot
columns_per_plot <- 4

# Split plot_list into chunks of 12
plot_chunks <- split(plot_list, ceiling(seq_along(plot_list)/lists_per_plot))

# Create a list to store the final plots
final_plots <- lapply(plot_chunks, function(chunk) {
  # Use wrap_plots to arrange the plots in a grid of 4 columns
  wrap_plots(chunk, ncol = columns_per_plot)
})

# Print or save each plot
for (i in seq_along(final_plots)) {
  print(final_plots[[i]])
  # You can also save each plot using ggsave
  ggsave(paste0("plot_", i, ".pdf"), final_plots[[i]],dpi = 300,width = 16, height = 10)
}
