####16_Target sensitivity analysis######
library(Seurat)
library(scRank)

###全部样本和细胞
load("/xxx/scobj_none_zhuanyi_zhengchang_cell.Rdata")

scobj@meta.data$cell_type<-scobj@active.ident


s16= scobj[, scobj@meta.data$orig.ident=="S16"]
obj_SPP1_s16<- CreateScRank(input = s16,
                            species = "human",
                            cell_type = "cell_type",
                            target=c("SPP1"))
obj_SPP1_s16<- scRank::Constr_net(obj_SPP1_s16, n.core = 30)
obj_SPP1_s16<- scRank::rank_celltype(obj_SPP1_s16)
obj_SPP1_s16 <- init_mod(obj_SPP1_s16, min_ModuleSize = 10)
plot_dim(obj_SPP1_s16,reductions ="harmony_umap",box_padding=0.1)
plot_net(obj_SPP1_s16, mode = "network", celltype ="Macrophages",vertex_label_cex = 0.5,charge = 0.1)#4.5*4.5
plot_net(obj_SPP1_s16, mode = "heatmap", celltype ="Macrophages")
drug_function_s16 <-plot_drug_function(obj_SPP1_s16,
                                       celltype ="Macrophages",
                                       category = "C5",#C5 : GO gene sets  #C7 : immunologic signatures
                                       top_number = 5,
                                       show_leading_edge = T)
drug_function_s16



s8= scobj[, scobj@meta.data$orig.ident=="S8"]
obj_SPP1_s8<- CreateScRank(input = s8,
                           species = "human",
                           cell_type = "cell_type",
                           target=c("SPP1"))
obj_SPP1_s8<- scRank::Constr_net(obj_SPP1_s8, n.core = 30)
obj_SPP1_s8<- scRank::rank_celltype(obj_SPP1_s8)
obj_SPP1_s8 <- init_mod(obj_SPP1_s8, min_ModuleSize = 10)
plot_dim(obj_SPP1_s8,reductions ="harmony_umap",box_padding=0.1)
plot_net(obj_SPP1_s8, mode = "network", celltype ="Macrophages",vertex_label_cex = 0.5,charge = 0.1)#4.5*4.5
plot_net(obj_SPP1_s8, mode = "heatmap", celltype ="Macrophages")

drug_function_s8 <-plot_drug_function(drug_function_s8,
                                      celltype ="Macrophages",
                                      category = "C5",#C5 : GO gene sets  #C7 : immunologic signatures
                                      top_number = 5,
                                      show_leading_edge = T)
drug_function_s8
