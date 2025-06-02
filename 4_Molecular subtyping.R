#######4_Molecular subtyping######



library("MOVICS")
library(readxl)
GSE212666_302 <- read_excel("/xxx/GSE212666-302.xlsx")
surv.info<- as.data.frame(GSE212666_302)
rownames(surv.info)<-surv.info$Name


###########3）获取聚类的最佳数量#######################
# identify optimal clustering number (may take a while)
optk.brca <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F), # note: the 4th data is somatic mutation which is a binary matrix
                         try.N.clust = 2:8, # try cluster number from 2 to 8
                         fig.path = getwd(),
                         fig.name    = "CLUSTER NUMBER")


iClusterBayes.res <- getiClusterBayes(data        = mo.data,
                                      N.clust     = 5,
                                      type        = c("gaussian","gaussian"),
                                      n.burnin    = 1000,
                                      n.draw      = 1000,
                                      prior.gamma = c(0.5, 0.5),
                                      sdev        = 0.05,
                                      thin        = 3)
############5）一次从多个算法中获取结果##################
# perform multi-omics integrative clustering with the rest of 10 algorithms
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"),
                         N.clust     = 5,
                         type        = c("gaussian", "gaussian"))
#时间很久
# attach iClusterBayes.res as a list using append() to moic.res.list with 9 results already
moic.res.list <- append(moic.res.list, 
                        list("iClusterBayes" = iClusterBayes.res))
# save moic.res.list to local path
save(moic.res.list, file = "moic_res_list_1000_cluster5_rep_20.rda")


#############6）从不同的算法中获得共识###################
cmoic.brca <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP_cluster5",
                               distance      = "euclidean",
                               linkage       = "average")


#############7） 使用 Silhoutte 量化相似度#################
getSilhouette(sil      = cmoic.brca$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "SILHOUETTE_cluster5",
              height   = 5.5,
              width    = 5)

#############8）基于聚类结果获取多组学热图###############

indata <- mo.data
indata$omics2  <- log2(indata$omics2 / (1 - indata$omics2))
# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2), # no truncation for mutation
                     centerFlag = c(T,T), # no center for mutation
                     scaleFlag  = c(T,T)) # no scale for mutation

feat   <- moic.res.list$iClusterBayes$feat.res
feat1  <- feat[which(feat$dataset == "omics1"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "omics2"),][1:10,"feature"]
annRow <- list(feat1, feat2)
# set color for each omics data
# if no color list specified all subheatmaps will be unified to green and red color pattern

omics1   <- c("#6699CC", "white"  , "#FF3C38")
omics2   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
col.list   <- list(omics1, omics2)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("omics1","omics2"),
             is.binary     = c(F,F), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.TMP","M value"),
             clust.res     = moic.res.list$ConsensusClustering$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")
# extract PAM50, pathologic stage and age for sample annotation
annCol    <- surv.info[,c("WHO", "Year", "Sex")]

# generate corresponding colors for sample annotation
annColors <- list(Year    = circlize::colorRamp2(breaks = c(min(annCol$Year),
                                                            median(annCol$Year),
                                                            max(annCol$Year)), 
                                                 colors = c("#0000AA", "#555555", "#AAAA00")),
                  WHO  = c("WHO_1" = "#A3C9D5",
                           "WHO_2"   = "#F6C63C",
                           "WHO_3"   = "red"),
                  Sex = c("Sex_1"    = "#E26844",
                          "Sex_2"    = "#66BC98"
                  ))

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA.TM","M value"),
             is.binary     = c(F,F), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.TMP","M value"),
             clust.res     = cmoic.brca$clust.res, # consensusMOIC results
             clust.dend    = NULL, # show no dendrogram for samples
             show.rownames = c(F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(F,F), # show no dendrogram for features
             annRow        = NULL, # no selected features
             color         = col.list,
             annCol        = annCol, # annotation for samples
             annColors     = annColors, # annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF CONSENSUSMOIC")



##########8）与其他亚型的一致性比较#########
# customize the factor level for pstage

a<-merge(surv.info,cmoic.brca[["clust.res"]],by=0)
row.names(a)<-a$Row.names
a$WHO <- factor(a$WHO, levels = c("WHO_1","WHO_2","WHO_3"))
a$clust <- factor(a$clust, levels = c("1","2","3","4","5"))
save(a,file="GSE212_cluster5.rda")

# agreement comparison (support up to 6 classifications include current subtype)
agree.brca <- compAgree(moic.res  = cmoic.brca,
                        subt2comp = a[,c(,"WHO")],
                        doPlot    = TRUE,
                        box.width = 0.2,
                        fig.name  = "AGREEMENT OF CONSENSUSMOIC WITH PAM50 AND PSTAGE")
#> --all samples matched.
print(agree.brca)


