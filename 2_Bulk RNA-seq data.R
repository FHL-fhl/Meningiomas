#####2_Bulk RNA-seq data#####


GSE212666_tpm <- as.data.frame(read_excel("/xxx/GSE212666_tpm.xlsx"))
rownames(GSE212666_tpm)<-GSE212666_tpm$target_id
GSE212666_tpm<-GSE212666_tpm[,-1]
colnames(GSE212666_tpm) <- sub(".*_", "", colnames(GSE212666_tpm))
save(GSE212666_tpm, file = "GSE212666_tpm.Rdata")

