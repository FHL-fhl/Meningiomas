####13_ScRNA-seq CNV#####
library(infercnv)
library(Seurat)
library(gplots)
library(ggplot2)
library(harmony)
library(dplyr)


load("/xxx/scobj_none_zhuanyi_zhengchang.Rdata")

##去掉Oligodendrocyte
sce= scobj[, Idents(scobj) %in% c( "Neoplastic_cell","T_cell","Macrophages","B_cell","Neutrophils","Dendritic Cells")]
dat <- as.matrix(GetAssayData(sce, slot='counts'))
dim(dat)


# 文件制作2：样本的描述
groupinfo_1 <- data.frame(v1 = colnames(dat),
                          v2 = sce@meta.data$orig.ident,
                          v3 = sce@active.ident)
groupinfo_1$v4 <- paste0(groupinfo_1$v2, "-", groupinfo_1$v3)
groupinfo<-data.frame(v1 = groupinfo_1$v1,
                      v2 = groupinfo_1$v4)C
rownames(groupinfo)<-groupinfo_1$v1

 
# 文件制作3：基因在染色体中的坐标
library(AnnoProbe)
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
# 去除性染色体（chrX 和 chrY）
non_sex_chromosomes <- geneInfor[!(geneInfor$chr %in% c("chrX", "chrY","chrM")), ]
# 自定义染色体排序（从 chr1 到 chr22）
chromosome_order <- c(paste0("chr", 1:22))
# 将 Chr 列转换为因子并按自定义顺序排序
non_sex_chromosomes$chr <- factor(non_sex_chromosomes$chr, levels = chromosome_order)
# 排序非性染色体
non_sex_chromosomes <- non_sex_chromosomes[order(non_sex_chromosomes$chr), ]
geneInfor<-non_sex_chromosomes[,c(1:4)]
dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),]
dim(dat)

expFile='expFile.txt'
#write.table(dat,file = expFile,sep = '\t',quote = F)
library(data.table)
fwrite(dat, file = expFile,  row.names = T,sep = "\t", quote = FALSE)


#groupFiles:是一个变量，存储写入的文件的文件名或路径。在这里文件名是groupFiles.txt。
groupFiles <- 'groupFiles.txt'
write.table(groupinfo,file = groupFiles, sep = '\t',
            quote = F, col.names = F, row.names = F)

#geneFile:是一个变量，存储写入的文件的文件名或路径。在这里文件名是geneFile.txt。
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)



#########所有细胞以免疫细胞作为参考###########
##注意作图时Myeloid_cells(including Macrophages,Neutrophil,DCs)
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = expFile,
                                     annotations_file = groupFiles,
                                     delim = "\t",
                                     gene_order_file = geneFile,
                                     ref_group_names = c(
                                                         "S3-Myeloid_cells","S3-T_cell","S3-B_cell",
                                                         "S4-Myeloid_cells","S4-T_cell","S4-B_cell",
                                                         "S5-Myeloid_cells","S5-T_cell","S5-B_cell",
                                                         "S6-Myeloid_cells","S6-T_cell","S6-B_cell",
                                                         "S7-Myeloid_cells","S7-T_cell","S7-B_cell",
                                                         "S8-Myeloid_cells","S8-T_cell","S8-B_cell",
                                                         "S9-Myeloid_cells","S9-T_cell","S9-B_cell",
                                                         "S10-Myeloid_cells","S10-T_cell","S10-B_cell",
                                                         "S11-Myeloid_cellss","S11-T_cell","S11-B_cell",
                                                         "S12-Myeloid_cells","S12-T_cell","S12-B_cell",
                                                         "S15-Myeloid_cells","S15-T_cell","S15-B_cell",
                                                         "S16-Myeloid_cells","S16-T_cell","S16-B_cell"))
#save(infercnv_obj,file = "infercnv_output_Neoplastic_and_immune_refimmune_infercnv_obj.Rdata")

#load("/project/hlfan/NO-meningiomas/picture_none_zhuanyi_zhengchang/infercnv_output_Neoplastic_and_immune_refimmune_infercnv_obj.Rdata")
infercnv_obj_run <- infercnv::run(infercnv_obj,
                                  cutoff =  0.1, #smart-seq选择1,10X选择0.1
                                  out_dir = "infercnv_output_Neoplastic_and_immune_refimmune", # dir is auto
                                  # 是否根据细胞注释文件的分组
                                  # 对肿瘤细胞进行分组
                                  # 影响read.dendrogram, 如果有多个细胞类型，且设置为TRUE，
                                  # 后续的read.dendrogram无法执行
                                  cluster_by_groups =  T, #是否根据患者类型(由细胞注释文件中定义)
                                  hclust_method = "ward.D2",# ward.D2 方法进行层次聚类
                                  analysis_mode = "samples", # 默认是samples，推荐是subclusters
                                  denoise = TRUE, # 去噪音
                                  HMM = F,  ##特别耗时间,是否要去背景噪音
                                  plot_steps = F, #不在每个步骤后生成图形。
                                  leiden_resolution = "auto", #可以手动调参数
                                  num_threads = 10 ,
                                  output_format = "pdf",
                                  plot_probabilities=FALSE)
