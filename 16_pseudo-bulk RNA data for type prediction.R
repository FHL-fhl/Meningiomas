####15_pseudo-bulk RNA data for type prediction####


###转化为伪bulk####

load("/xxx/scobj_none_zhuanyi_zhengchang_cell.Rdata")
counts <- scobj@assays$RNA@counts
metadata<-scobj@meta.data
#将cluster信息加到metadata
metadata$cluster_id<-factor(scobj@active.ident)

#建Vsingle cell experiment对象
sce<-SingleCellExperiment(assays = list(counts = counts),colData = metadata)

#确定用于聚合数据的cluster
groups <- colData(sce)[,c("cluster_id","orig.ident")]

# 提取细胞类型
kids <-purrr::set_names(levels(groups$cluster_id))
#cluster数目
nk<- length(kids)

#提取样本名称
groups$orig.ident <- as.factor(groups$orig.ident)
sids <- purrr::set_names(levels(groups$orig.ident))

# 样本数目
ns<-length(sids)

#确定每个样本细胞数
table(sce$orig.ident)
#转为数字向量
n_cells<-as.numeric(table(sce$orig.ident))
#确定每个样本出现的细胞位置
m<-match(sids,sce$orig.ident)
# 构建样本级别metadata
ei <- data.frame(colData(sce)[m,],
                 n_cells,row.names = NULL)%>%select(-"cluster_id")

#聚合count到样本水平
#细胞信息提取
groups<-colData(sce)[,c("orig.ident")]
#跨样本整合(求和)
library(Matrix.utils)
pb<- aggregate.Matrix(t(counts(sce)),
                      groupings =groups,fun="sum")
pb<-as.data.frame(pb)

scRNA_count<-t(pb)
###将count转化为tpm
#devtools::install_github("IOBR/IOBR")
library(IOBR)
scRNA_tpm=count2tpm(scRNA_count,idType= "SYMBOL")
save(scRNA_tmp,file="scRNA_tpm.Rdata")




####scRNA Sample 进行分类预测####

library("MOVICS")

# run DEA with limma
load("/xxx/mo.data.rda")
load("/xxx/moic_res_list_1000_cluster5_rep_20.rda")
cmoic.brca <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP_cluster5",
                               distance      = "euclidean",
                               linkage       = "average")
load("/xxx/GSE212666_tpm.rda")
tmp<-GSE302_tpm


runDEA(dea.method = "limma",
       expr       = tmp, # normalized expression data
       moic.res   = cmoic.brca,
       prefix     = "meningioma")


marker.up <- runMarker(moic.res      = cmoic.brca,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "meningioma", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 200, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = tmp, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")


yau.ntp.pred <- runNTP(expr       = scRNA_tmp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR SCRNA-seq YAU UP") 

