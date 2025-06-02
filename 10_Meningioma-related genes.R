####10_Meningioma-related genes####
library(dplyr)
library(ggplot2)
library(ggpubr)


GSE212666_beta <- readRDS("/xxx/GSE212666_beta.rds")


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



###NF2####

NF2_expression <- GSE302_tpm["NF2", ]
# 将表达数据与样本分组信息合并
expression_data <- data.frame(
  SampleID = colnames(GSE302_tpm),
  NF2_expression = as.numeric(NF2_expression)
)
expression_data <- merge(expression_data, an, by.x = "SampleID", by.y = 0)

# 画箱线图
pdf(file = "NF2.pdf",height = 5,width = 3.5)
ggplot(expression_data, aes(x = factor(cluster), y = NF2_expression, fill = cluster)) +
  geom_boxplot() +
  labs(title = "NF2 ") +
  scale_fill_manual(values = c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB")) +
  stat_compare_means(comparisons = comparisons, 
                     method = "t.test", 
                     label = "p.signif") +
  theme_classic()
dev.off()



####SMARCB1###
SMARCB1_expression <- GSE302_tpm["SMARCB1", ]
# 将表达数据与样本分组信息合并
expression_data <- data.frame(
  SampleID = colnames(GSE302_tpm),
  SMARCB1_expression = as.numeric(SMARCB1_expression)
)
expression_data <- merge(expression_data, an, by.x = "SampleID", by.y = 0)

# 画箱线图
pdf(file = "SMARCB1.pdf",height = 5,width = 3.5)
ggplot(expression_data, aes(x = factor(cluster), y = SMARCB1_expression, fill = cluster)) +
  geom_boxplot() +
  labs(title = "SMARCB1 ") +
  scale_fill_manual(values = c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB")) +
  stat_compare_means(comparisons = comparisons, 
                     method = "t.test", 
                     label = "p.signif") +
  theme_classic()
dev.off()


