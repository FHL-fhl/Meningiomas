####7_Differential analysis and functional annotation####

library("MOVICS")



####DGEs####
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

# choose limma result to identify subtype-specific up-regulated biomarkers
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
# choose limma result to identify subtype-specific up-regulated biomarkers
marker.down <- runMarker(moic.res      = cmoic.brca,
                         dea.method    = "limma", # name of DEA method
                         prefix        = "meningioma", # MUST be the same of argument in runDEA()
                         dat.path      = getwd(), # path of DEA files
                         res.path      = getwd(), # path to save marker files
                         p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                         p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                         dirct         = "down", # direction of dysregulation in expression
                         n.marker      = 50, # number of biomarkers for each subtype
                         doplot        = TRUE, # generate diagonal heatmap
                         norm.expr     = tmp, # use normalized expression as heatmap input
                         annCol        = annCol, # sample annotation in heatmap
                         annColors     = annColors, # colors for sample annotation
                         show_rownames = FALSE, # show no rownames (biomarker name)
                         fig.name      = "DownREGULATED BIOMARKER HEATMAP")

#run gene set variation analysis

GSET.FILE <- 
  system.file("extdata", "gene sets of interest.gmt", package = "MOVICS", mustWork = TRUE)

gsva.res <- 
  runGSVA(moic.res      = cmoic.brca,
          norm.expr     = tmp,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          annCol        = annCol,
          annColors     = annColors,
          fig.path      = getwd(),
          fig.name      = "GENE SETS OF INTEREST HEATMAP",
          height        = 5,
          width         = 8)


#####DMPs######
library(dplyr)
library(readr)
library(dml)
library(sesame)
library(SummarizedExperiment)
sesameDataCache() # required at new sesame installation


load("/xxx/GSE212_cluster5.rda")
GSE212666_beta <- readRDS("/xxx/GSE212666_beta.rds")


betas <- sapply(GSE212666_beta, as.numeric)
is.numeric(betas)
betas <- as.data.frame(betas)
rownames(betas)<-rownames(GSE212666_beta)

#设置参考
an =as.data.frame(a[,c(6:8,10)])
rownames(an)<-rownames(a)
colnames(an)<-c("WHO","Sex","Age","Cluster")
an$ID <- rownames(an)

# 假设 an 是你的数据框
an <- an %>%
  mutate(Cluster = recode(Cluster, 
                          `1` = "CS1", 
                          `2` = "CS2", 
                          `3` = "CS3", 
                          `4` = "CS4", 
                          `5` = "CS5"))



an$Sex = relevel(factor(an$Sex),ref="Sex_1")
an$WHO=relevel(factor(an$WHO),ref="WHO_1")


an$Group_CS1 <- ifelse(an$Cluster == "CS1", "group_1", "group_2")
an$Group_CS2 <- ifelse(an$Cluster == "CS2", "group_1", "group_2")
an$Group_CS3 <- ifelse(an$Cluster == "CS3", "group_1", "group_2")
an$Group_CS4 <- ifelse(an$Cluster == "CS4", "group_1", "group_2")
an$Group_CS5 <- ifelse(an$Cluster == "CS5", "group_1", "group_2")

an$Group_CS1=relevel(factor(an$Group_CS1),ref="group_2")
an$Group_CS2=relevel(factor(an$Group_CS2),ref="group_2")
an$Group_CS3=relevel(factor(an$Group_CS3),ref="group_2")
an$Group_CS4=relevel(factor(an$Group_CS4),ref="group_2")
an$Group_CS5=relevel(factor(an$Group_CS5),ref="group_2")



# 假设 Group_CS1 至 Group_CS5 是你想要比较的分组变量名称
# 假设 an 是包含 WHO、Sex 和分组信息的数据框

# 循环计算每个分组的 smry 并保存
for (i in 1:5) {
  group_name <- paste0("Group_CS", i)  # 生成分组名称
  formula_str <- paste0("~", group_name, "+WHO+Sex") # 创建 R 公式字符串
  
  # 数据清洗，过滤 NA
  str(an)
  ok1 <- checkLevels(betas, an[, group_name])  # 假设分组变量在 an 数据框中
  sum(ok1)
  betas <- betas[ok1, ]
  dim(betas)
  
  # 将模拟处理组织和性别的 DNA 甲基化变异作为协变量
  smry <- DML(betas, formula(formula_str), meta=an)
  smry
  smry[[1]]
  
  # 保存结果
  saveRDS(smry, file = paste0(group_name, "_VS_other.rds"))
}



#####以CS5为例####
CS5_VS_other <- readRDS("/xxx/CS5_VS_other.rds")
test_result = summaryExtractTest(CS5_VS_other)
colnames(test_result) # the column names, show four groups of statistics
rownames(test_result)<-rownames(betas)
test_result$Probe_ID<-rownames(test_result)


###特异性探针选取
library(dplyr)
library(tidyr)
test_result %>% dplyr::filter(FPval_Group_CS5< 0.05, Eff_Group_CS5> 0.25) %>%
  dplyr::select(FPval_Group_CS5,Eff_Group_CS5)
#使用 0.1 作为效应大小阈值。
#这意味着DNA 甲基化差异低于 0.1 （10%） 被认为是非生物学差异的

test_result %>%
  mutate(Sex_specific =
           ifelse(FPval_Sex < 0.05 & Eff_Sex > 0.1, TRUE, FALSE)) %>%
  mutate(WHO_specific =
           ifelse(FPval_WHO < 0.05 & Eff_WHO > 0.1, TRUE, FALSE)) %>%
  mutate(Group_CS5_specific =
           ifelse(FPval_Group_CS5 < 0.05 & Eff_Group_CS5> 0.1, TRUE, FALSE)) %>%
  dplyr::select(Sex_specific, WHO_specific, Group_CS5_specific) %>% table



DMP_CS5<-test_result %>% dplyr::filter(Pval_Group_CS5group_1 < 0.05, abs(Est_Group_CS5group_1) > 0.25) %>%
  dplyr::select(Probe_ID, Pval_Group_CS5group_1,Est_Group_CS5group_1)


##差异甲基化位点基因注释
load("/xxx/probe_anno_epic.Rdata")

DMP_CS5anno<-merge(DMP_CS5,probe.features,by.x="Probe_ID",by.y=0)
DMP_CS5anno <- DMP_CS5anno[,c(1:4,8,9)]



DMP_CS5betas <- merge(DMP_CS5anno,betas,by.x="Probe_ID",by.y=0)
DMP_CS5betas <- DMP_CS5betas[,-c(1:14)]

#转化为数值型
DMP_CS5betas<-as.data.frame(DMP_CS5betas)
DMP_CS5betas<-sapply(DMP_CS5betas, as.numeric)
is.numeric(DMP_CS5betas) #TRUE





