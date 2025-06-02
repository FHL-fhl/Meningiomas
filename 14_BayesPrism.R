####BayesPrism#####


library(scDblFinder)
library(Seurat)
library(dplyr)
library(cowplot)
library(sctransform)
library(SoupX)
library(glmGamPoi)
library(CellChat)
library(tidyverse)
library(ggalluvial)
library(SeuratObject)


library(BayesPrism)

load("XXX/scobj_none_zhuanyi_zhengchang_cell.Rdata")
sc <- subset(
  scobj,
  idents = c("Neoplastic_cell", "Endothelial", "Macrophages", "T_cell", 
             "Fibroblast","Dendritic Cells","B_cell","Neutrophils","Oligodendrocyte")
)


# DimPlot(sc.dat, reduction = "umap", label = TRUE, pt.size = 0.5) 
# rm(scobj)

sc.dat <- as.data.frame(t(sc[["RNA"]]@counts))
dim(sc.dat) #[1] 87495 22799
head(rownames(sc.dat))
head(colnames(sc.dat))

load("/project/hlfan/MOVICS_0727/GSE212666_tpm.rda")
T_NO1 <-t(GSE302_tpm)
bk.dat   <- t(GSE302_tpm)     
bk.dat <- as.data.frame(apply(bk.dat, 2, as.numeric)   )         
rownames(bk.dat) <-rownames(T_NO1) 


sc@meta.data$seurat_clusters <- droplevels(sc@meta.data$seurat_clusters)
table(sc@meta.data[["seurat_clusters"]])
cell.state.labels <- sc@meta.data$seurat_clusters




sc@meta.data[["cell_type"]] <- sc@active.ident
cell.type.labels <- sc@meta.data[["cell_type"]]




table(cbind.data.frame(cell.state.labels,cell.type.labels))


##QC of cell type and state labels
plot.cor.phi (input=sc.dat,
              input.labels=cell.state.labels,
              title="cell_state_correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.cs", 
              cexRow=0.2, cexCol=0.2,
              margins=c(2,2))
plot.cor.phi (input=sc.dat, 
              input.labels=cell.type.labels, 
              title="cell type correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.ct",
              cexRow=0.5, cexCol=0.5,
)


sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting. 
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)



bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)




### Filter outlier genes from scRNA-seq data
sc.dat.filtered <- cleanup.genes(input = sc.dat, input.type = "count.matrix", species = "hs",
                                 gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"), exp.cells = 5)
dim(sc.dat.filtered)





#note this function only works for human data. For other species, you are advised to make plots by yourself.

plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)





sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered, gene.type = "protein_coding")

diff.exp.stat <- get.exp.stat(sc.dat=sc.dat[,colSums(sc.dat>0)>3],# filter genes to reduce memory use
                              cell.type.labels=cell.type.labels,
                              cell.state.labels=cell.state.labels,
                              pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                              cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                              n.cores=1 #number of threads
)


sc.dat.filtered.pc.sig <- select.marker (sc.dat=sc.dat.filtered.pc,
                                         stat=diff.exp.stat,
                                         pval.max=0.01,
                                         lfc.min=0.1)


myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key="Neoplastic_cell",# 
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

bp.res <- run.prism(prism = myPrism, n.cores=50)
bp.res#结果
slotNames(bp.res)

#提取细胞类型
theta <- get.fraction(bp=bp.res,
                      which.theta="final",
                      state.or.type="type")

head(theta)
write.csv(theta,file="theta.csv")

#提取变异系数
theta.cv <- bp.res@posterior.theta_f@theta.cv
head(theta.cv)


library(reshape2)
#head(t(Z.tumor[1:5,]))
ratio <- as.data.frame(theta)
ratio <- t(ratio)
ratio <- as.data.frame(ratio)
ratio <- tibble::rownames_to_column(ratio)
ratio <- melt(ratio)
colourCount = length(ratio$rowname)
ggplot(ratio) + 
  geom_bar(aes(x = variable,y = value,fill = rowname),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))


load("/XXX/GSE212_cluster5.rda")


load("/XXX/theta.Rdata")
load("/XXX/GSE212_cluster5.rda")



# 加载必要的包
library(pheatmap)
library(RColorBrewer)
library(dplyr)

## 1. 数据准备 -------------------------------------------------------------

# 假设a是包含cluster信息的数据框
an <- data.frame(cluster = a[, 10], row.names = rownames(a))

# 重命名cluster
an <- an %>%
  mutate(cluster = recode(cluster, 
                          `1` = "CS1", 
                          `2` = "CS2", 
                          `3` = "CS3", 
                          `4` = "CS4", 
                          `5` = "CS5")) %>%
  arrange(cluster)

# 获取细胞类型信息
cell_types <- Idents(scobj)
unique_cell_types <- levels(cell_types)

## 2. 注释数据准备 ---------------------------------------------------------

# 列注释（位于热图顶端）- 按CLUSTER分组
annCol <- data.frame(CLUSTER = an$cluster,
                     row.names = rownames(an),
                     stringsAsFactors = FALSE)

# 行注释（位于热图左侧）- 按Methods分组
annRow <- data.frame(
  Methods = factor(cell_types[rownames(an)], levels = unique_cell_types),
  row.names = rownames(an),
  stringsAsFactors = FALSE
)

# 确保plotdata的行列顺序与注释一致
plotdata <- t(theta[rownames(an), ] ) # 假设theta是表达矩阵

## 3. 颜色设置 -------------------------------------------------------------

# 为Methods设置颜色
methods.col <- brewer.pal(n = length(unique_cell_types), name = "Paired")
names(methods.col) <- unique_cell_types

# 为CLUSTER设置颜色
cluster.col <- c("CS1" = "#2EC4B6", "CS2" = "#E71D36",
                 "CS3" = "#FF9F1C", "CS4" = "#BDD5EA",
                 "CS5" = "#FFA5AB")

annColors <- list(
  Methods = methods.col,
  CLUSTER = cluster.col
)

## 4. 数据标准化 -----------------------------------------------------------

standarize.fun <- function(indata, halfwidth = NULL, centerFlag = TRUE, scaleFlag = TRUE) {  
  outdata <- t(scale(t(indata), center = centerFlag, scale = scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata > halfwidth] <- halfwidth
    outdata[outdata < (-halfwidth)] <- -halfwidth
  }
  return(outdata)
}

plotdata <- standarize.fun(plotdata, halfwidth = 2)

## 5. 计算分割位置 --------------------------------------------------------

# 列分割位置（按CLUSTER）
col_counts <- table(annCol$CLUSTER)
col_gaps <- cumsum(col_counts[-length(col_counts)])

# 行分割位置（按Methods）
row_counts <- table(annRow$Methods)
row_gaps <- cumsum(row_counts[-length(row_counts)])

# 安全验证：确保分割位置不超过矩阵尺寸
col_gaps <- col_gaps[col_gaps < ncol(plotdata)]
row_gaps <- row_gaps[row_gaps < nrow(plotdata)]

## 6. 绘制热图 ------------------------------------------------------------

pheatmap::pheatmap(
  mat = as.matrix(plotdata),
  border_color = NA,
  color = colorRampPalette(c("blue", "white", "red"))(64),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = annCol,
  annotation_row = annRow,
  annotation_colors = annColors,
  gaps_col = col_gaps,  # 横向按CLUSTER分割
  gaps_row = row_gaps,  # 纵向按Methods分割
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Cell Type Expression Heatmap"
)
