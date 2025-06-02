######11_Drug sensitivity######
library("MOVICS")
library(readxl)

load("/xxx/GSE212666_tpm.rda")
tpm<-GSE302_tpm
load("/xxx/moic_res_list_1000_cluster5_rep_20.rda")
cmoic.brca <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP_cluster5",
                               distance      = "euclidean",
                               linkage       = "average")


load("/xxx/possibleDrugs2016.rda")

# 假设 `a` 中包含所有药物的名称
a <- unique(drugData2016$Drug.name)

# 创建一个空列表来存储结果
results_list <- list()

# 循环遍历每个药物名称
for (drug in a) {
  # 动态生成变量名称
  result_name <- paste0("drug.brca_", drug)
  
  # 调用 compDrugsen 函数
  results_list[[result_name]] <- compDrugsen(moic.res    = cmoic.brca,
                                             norm.expr   = tmp,
                                             drugs       = c(drug), # 替换药物名称
                                             tissueType  = "all", 
                                             test.method = "nonparametric", 
                                             prefix      = "BOXVIOLIN OF ESTIMATED IC50")
}

save(results_list,file = "Drug_result_list.rda")
