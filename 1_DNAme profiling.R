###### 1_DNAme profiling######

library(sesame)
library(parallel)
library(DNAcopy)
sesameDataCache()

###sesame数据预处理得到clean_betas
{
  
  #读入SigDF对象并保存QC结果
  s<-mclapply(searchIDATprefixes("."), readIDATpair)
  saveRDS(s, file = "s.rds")
  qcs = openSesame(s, prep="", func=sesameQC_calcStats, funs="detection")
  head(do.call(rbind, lapply(qcs, as.data.frame)))
  
  
  
  #对一组 IDAT 文件进行数据处理，
  #包括读取、预处理、背景校正、染料偏差校正和计算 β 值等步骤，
  #并将处理后的 β 值按列合并成一个矩阵
  idat_dir="/xxx/GSE212666"
  
  betas = do.call(cbind, BiocParallel::bplapply(
    searchIDATprefixes(idat_dir), function(pfx) {
      getBetas(noob(pOOBAH(dyeBiasNL(inferInfiniumIChannel(qualityMask(
        readIDATpair(pfx)))))))
    }, BPPARAM = BiocParallel::MulticoreParam(2)))
  
  #saveRDS(betas, file = "betas.rds")

  sum(rowSums(is.na(betas)) > 0)  #248190 探针被过滤
  # 找到不包含NA值的行的索引
  clean_betas <- betas[complete.cases(betas), ]
  nrow(clean_betas) #618363探针被保留
  colnames(clean_betas) <- substr(colnames(clean_betas), 1, 10)
  saveRDS(clean_betas, file =  "GSE212666_beta.rds")
}




#cd /project/hlfan/meningiomas_methylation/sesame_code
#nohup Rscript segs.R &



s<-readRDS( "s.rds")
segs <- lapply(s, cnSegmentation)
saveRDS(segs, file = "/xxx/segs.rds")

