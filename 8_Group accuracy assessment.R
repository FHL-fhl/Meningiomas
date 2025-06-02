#####8_Group accuracy assessment#######

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


####检测准确定####
# predict subtype  using NTP
ntp.pred <- runNTP(expr      = tmp,
                        templates = marker.up$templates,
                        doPlot    = FALSE) 


# predict subtype cohort using PAM
pam.pred <- runPAM(train.expr  = tmp,
                        moic.res    = cmoic.brca,
                        test.expr   = tmp)


# check consistency between current and NTP-predicted subtype 
runKappa(subt1     = cmoic.brca$clust.res$clust,
         subt2     = ntp.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "NTP")

# check consistency between current and PAM-predicted subtype 
runKappa(subt1     = cmoic.brca$clust.res$clust,
         subt2     = pam.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "PAM")


