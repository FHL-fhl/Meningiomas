#####9_Estimation of CNV######


s<-readRDS( "s.rds")
segs <- lapply(s, cnSegmentation)
#saveRDS(segs, file = "/project/hlfan/meningiomas_methylation/GSE212666/segs.rds")


#########5）比较基因组改变的分数########
#segs <- readRDS("/project/hlfan/meningiomas_methylation/GSE212666/segs.rds")
segment <- data.frame()
for (sample_name in names(segs)) {
  # 获取当前样本的 seg 数据
  seg_data <- segs[[sample_name]][["seg.signals"]]
  # 删除不需要的列
  seg_data <- seg_data[,-c(7:12)]
  # 添加样本名称列
  seg_data$ID <- sample_name
  # 将当前样本的数据添加到合并后的数据框中
  segment <- rbind(segment, seg_data)
}
segment<-segment[,-5]
colnames(segment) <- c("sample","chrom","start","end","value")
segment$ID<-sub("_[^_]+$", "", segment$sample)
a<-merge(segment,GSE212666_302,by.x="ID",by.y="Meth_ID")
a$sample<-a$Name
segment<-a[,c(2:6)]
rm(a)
# compare FGA, FGG, and FGL
fga.brca <- compFGA(moic.res     = cmoic.brca,
                    segment      = segment,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.2, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "BARPLOT OF FGA")
head(fga.brca$summary)
#save(fga.brca,file="fga.brca.rda")
df<-fga.brca[["summary"]]


library(ggplot2)
library(ggpubr)
library(RColorBrewer)
# 指定两两比较的组
comparisons <- list(c("CS1", "CS2"),
                    c("CS1", "CS3"),
                    c("CS1", "CS4"),
                    c("CS1", "CS5"),
                    c("CS2", "CS3"),
                    c("CS2", "CS4"),
                    c("CS2", "CS5"),
                    c("CS3", "CS4"),
                    c("CS3", "CS5"),
                    c("CS4", "CS5"))

# 使用 Set1 作为配色方案
color_palette <- brewer.pal(n = 5, name = "Set1")

# FGA的箱线图
p1 <- ggplot(df, aes(x = Subtype, y = FGA, fill = Subtype)) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette) +
  stat_compare_means(comparisons = comparisons, 
                     method = "t.test", 
                     label = "p.signif") +
  theme_classic() +
  labs(title = "FGA across Subtypes", y = "FGA")

# FGG的箱线图
p2 <- ggplot(df, aes(x = Subtype, y = FGG, fill = Subtype)) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette) +
  stat_compare_means(comparisons = comparisons, 
                     method = "t.test", 
                     label = "p.signif") +
  theme_classic() +
  labs(title = "FGG across Subtypes", y = "FGG")

# FGL的箱线图
p3 <- ggplot(df, aes(x = Subtype, y = FGL, fill = Subtype)) +
  geom_boxplot() +
  scale_fill_manual(values = color_palette) +
  stat_compare_means(comparisons = comparisons, 
                     method = "t.test", 
                     label = "p.signif") +
  theme_classic() +
  labs(title = "FGL across Subtypes", y = "FGL")

# 使用ggpubr将所有图组合到一起
ggarrange(p1, p2, p3, ncol = 3, nrow = 1)




####IGV可视化文件整理####
library(readxl)
segs <- readRDS("/xxx/segs.rds")

GSE212666_302 <- read_excel("/xxx/GSE212666-302.xlsx")


segs <- readRDS("/xxx/segs.rds")
names(segs) <- substr(names(segs), 1, nchar(names(segs)) - 7)



load("/xxx/GSE212_cluster5.rda")

an =as.data.frame(a[,c(5,10)])

rownames(an)<-NULL
colnames(an)<-c("ID","cluster")
library(dplyr)

# 假设 an 是你的数据框
an <- an %>%
  mutate(cluster = recode(cluster, 
                          `1` = "CS1", 
                          `2` = "CS2", 
                          `3` = "CS3", 
                          `4` = "CS4", 
                          `5` = "CS5"))

####删除XY和13，14，15，21，22的p臂
merged_data_list <- list()
# 顶心染色体的p臂范围（假设参考的是人类基因组 hg38 版本）
acrocentric_p_arm_ranges <- list(
  "13" = c(0, 16000000),   # 13号染色体 p臂估计范围
  "14" = c(0, 16000000),   # 14号染色体 p臂估计范围
  "15" = c(0, 17083673),   # 15号染色体 p臂估计范围
  "21" = c(0, 10864560),   # 21号染色体 p臂估计范围
  "22" = c(0, 12954788)    # 22号染色体 p臂估计范围
)

# 遍历每个样本分组
for (group in unique(an$cluster)) {
  # 获取当前分组的样本列表
  samples <- an$ID[an$cluster == group]
  
  # 创建一个空的数据框来存储当前分组合并后的数据
  merged_data <- data.frame()
  
  # 遍历当前分组的每个样本
  for (sample_name in samples) {
    # 获取当前样本的 seg 数据
    seg_data <- segs[[sample_name]][["seg.signals"]]
    
    # 删除不需要的列
    seg_data <- seg_data[,-c(7:12)]
    
    # 截取染色体名称的前三个字符
    seg_data$chrom <- substr(seg_data$chrom, 4, nchar(seg_data$chrom))
    
    # 更改列名
    colnames(seg_data) <- c("Sample", "chromosome", "Start", "End", "Num_Probe", "Segment_mean")
    
    # 添加样本名称列
    seg_data$Sample <- sample_name
    
    # 去除X和Y染色体
    seg_data <- subset(seg_data, !(chromosome %in% c("X", "Y")))
    
    # 去除顶心染色体的p臂范围
    for (chrom in names(acrocentric_p_arm_ranges)) {
      p_arm_range <- acrocentric_p_arm_ranges[[chrom]]
      seg_data <- seg_data[!(seg_data$chromosome == chrom & seg_data$End <= p_arm_range[2]), ]
    }
    
    # 将当前样本的数据添加到当前分组合并后的数据框中
    merged_data <- rbind(merged_data, seg_data)
  }
  
  # 将当前分组合并后的数据添加到列表中
  merged_data_list[[paste("group_", group, "_merged_seg_data", sep = "")]] <- merged_data
}

# 将每个分组的合并后的数据写入到文件中
for (group_name in names(merged_data_list)) {
  write.table(merged_data_list[[group_name]], paste(group_name, ".seg", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
}



# 定义要查找NF2范围的染色体和范围
target_chromosome <- "22"
target_start <- 0  # 你可以根据需要调整范围
target_end <- Inf 
# 创建一个列表来存储每组中符合条件的变化
filtered_changes_list <- list()
# 遍历每个分组的seg数据
for (group_name in names(merged_data_list)) {
  seg_data <- merged_data_list[[group_name]]
  
  # 筛选出目标范围内的变化
  filtered_changes <- subset(seg_data, chromosome == target_chromosome & 
                               Start <= target_end & 
                               End >= target_start& 
                               Segment_mean < -0.1)
  
  # 将筛选出的变化添加到列表中
  filtered_changes_list[[paste("group_", group_name, "CHR22", sep = "")]] <- filtered_changes
}
# 将每组筛选出的变化写入到文件中
for (group_name in names(filtered_changes_list)) {
  write.table(filtered_changes_list[[group_name]], paste(group_name, "_CHR22.seg", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
}


# 定义要查找的染色体22和数量
target_chromosome <- "22"
target_start <- 0  # 你可以根据需要调整范围
target_end <- Inf  # 设置为Inf表示不限结束位置
#平均强度值小于 -0.1 的片段定义为拷贝数丢失
# 创建一个列表来存储每个组中具有缺失的样本数量
group_deletions_count <- list()

# 遍历每个分组的seg数据
for (group_name in names(merged_data_list)) {
  seg_data <- merged_data_list[[group_name]]
  
  # 筛选出目标范围内的缺失变化（假设Segment_mean < 0表示缺失）
  filtered_changes <- subset(seg_data, chromosome == target_chromosome & 
                               Start <= target_end & 
                               End >= target_start & 
                               Segment_mean < -0.1)
  
  # 获取具有缺失的样本
  samples_with_deletions <- unique(filtered_changes$Sample)
  
  # 统计当前组中具有缺失的样本数量
  group_deletions_count[[group_name]] <- length(samples_with_deletions)
}

# 打印每个组中具有22号染色体缺失的样本数量
for (group_name in names(group_deletions_count)) {
  cat("组", group_name, "中具有22号染色体缺失的样本个数:", group_deletions_count[[group_name]], "\n")
}
