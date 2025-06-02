######3_PCA#####


library("MOVICS")
library(readxl)

####加载beta矩阵和配对TPM表达矩阵###
GSE212666_beta <- readRDS("/xxx/GSE212666_beta.rds")
load("/xxx/GSE212666_tpm.rda")
tpm<-GSE302_tpm


###加载302个样本信息####
GSE212666_302 <- read_excel("/xxx/GSE212666-302.xlsx")
surv.info<- as.data.frame(GSE212666_302)
rownames(surv.info)<-surv.info$Name


# #######表达矩阵tpmPCA_TOP1000提取#######
pca_result_tpm <- prcomp(t(tpm), center = TRUE, scale. = FALSE)

#提取pca_result_tpm前12个主成分top1000
top_12_pcs <- pca_result_tpm$rotation[, 1:12]
num_sites <- 1000
top_Gene_sites <- list()

# 遍历每个主成分
for (i in 1:12) {
  # 计算绝对权重并排序，选择前num_sites个甲基化位点
  abs_weights <- abs(top_12_pcs[, i])
  sorted_indices <- order(-abs_weights)[1:num_sites]
  
  # 存储甲基化位点信息
  top_Gene_sites[[paste0("PC", i)]] <- data.frame(
    Gene_Site = rownames(top_12_pcs)[sorted_indices],
    Weight = top_12_pcs[sorted_indices, i]
  )
}

# 合并所有PC的甲基化位点并去重
merged_sites <- unlist(lapply(top_Gene_sites, function(x) x$Gene_Site))
unique_sites <- unique(merged_sites) 
# 匹配TPM数据并保存
matched_tpm <- tpm[rownames(tpm) %in% unique_sites, ]####2312genes
#save(matched_tpm, file = "matched_tpm_1000.rda")






# #######甲基化矩阵beta_PCA_TOP500提取#######
meth.beta<- GSE212666_beta 
pca_result_beta <- prcomp(t(meth.beta), center = TRUE, scale. = FALSE)
summary(pca_result_beta )

# 计算 PCA
top_12_pcs <- pca_result_beta$rotation[, 1:12]
merged_column <- list()
# 初始化一个列表来存储每个 PC 的前500个甲基化位点
top_1000_methylation_sites <- list()
# 循环处理每个主成分
for (i in 1:12) {
  # 计算绝对权重并排序，选择前1000个甲基化位点
  abs_weights <- abs(top_12_pcs[, i])
  sorted_indices <- order(-abs_weights)[1:500]
  # 创建数据框并存储在列表中
  top_1000_methylation_sites[[paste("PC", i, sep="")]] <- data.frame(
    Methylation_Site = rownames(top_12_pcs)[sorted_indices],
    Weight = top_12_pcs[sorted_indices, i]
  )
  # 提取当前 PC 的前1000个甲基化位点并合并到 merged_column 列表中
  pc <- top_1000_methylation_sites[[paste("PC", i, sep="")]]$Methylation_Site
  merged_column <- c(merged_column, pc)
}
# 将 merged_column 列表转换为向量
merged_column <- unlist(merged_column)
# 获取唯一值
unique_merged_column <- unique(merged_column)####5334 probes
matched_meth.beta<- meth.beta[row.names(meth.beta) %in% unique_merged_column, ]
matched_meth.beta<-as.data.frame(matched_meth.beta)

mo.data <- list(omics1 = matched_tpm,
                omics2 = matched_meth.beta)

save(mo.data, file = "/xxx/mo.data.rda")


