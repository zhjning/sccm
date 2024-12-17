## load data
dataDir = "/Volumes/ZJN_backup/BWH_backup220923/MCAO/public_datasets/result.sc.20220725/GSE174574/combined/"
so = loadRDS(file.path(dataDir,grep("sos",list.files(dataDir),value=T)))
so = Seurat::UpdateSeuratObject(so)
## umap-based analysis
seurat_labels = named_vec(c(0:20),c("Endo","Mig1","aEndo","EPC1","SMC","Olig","Ast","Mac","Mono","Neu","vEndo","Peri","capEndo","EPC2","Mig2","T","Fib1","Fib2","Neuron","Peri","RBC"))
so$seurat_labels = seurat_labels[so$seurat_clusters %>% as.character]
p = DimPlot(so, group.by = "seurat_labels", label = T) + NoLegend()
## extract the cellchat input files from seurat v3 object
so <- subset(so, seurat_labels %in% c("Endo","vEndo","aEndo","Mig1","Mig2","Mac","Mono","SMC","Peri","Neu","Ast","capEndo","Fib1","Fib2"))
so$ct.major <- so$seurat_labels
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
labels <- so$ct.major
meta <- data.frame(ct = labels, row.names = names(labels), sample = so$orig.ident)
ident_id = unique(so$orig.ident)[1]


## Functional setting
require(msigdbr)
funcDB = msigdbr::msigdbr(species = "Mus musculus",category = "C5")
funcToView = funcDB %>%
  filter( gs_name %in% (grep("PROLIFERATION",gs_name,value=T) %>% unique)) %>%
  filter( gs_name %in% (grep("ENDO",gs_name,value=T) %>% unique))
pfuncToUse = c("GOBP_ENDOTHELIAL_CELL_PROLIFERATION",
               "GOBP_VASCULAR_ENDOTHELIAL_CELL_PROLIFERATION",
               "GOBP_POSITIVE_REGULATION_OF_ENDOTHELIAL_CELL_PROLIFERATION" ,
               "GOBP_POSITIVE_REGULATION_OF_BLOOD_VESSEL_ENDOTHELIAL_CELL_PROLIFERATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS",
               "GOBP_POSITIVE_REGULATION_OF_VASCULAR_ENDOTHELIAL_CELL_PROLIFERATION")
funcToUse = "GOBP_ENDOTHELIAL_CELL_PROLIFERATION"
funcGeneSets = funcDB %>% filter(gs_name == funcToUse) %>% pull(gene_symbol) %>% unique

oDir <- file.path(dataDir, "cellchat.out")
# Cell relationships setting
feature.name = "ct.major"
# ctlist_to_check = c("Mono,Mig1,aEndo","Mono,Mig2,aEndo","Mono,Mac,aEndo","Mono,Mac,SMC","Mono,aEndo,SMC","Mono,SMC,aEndo",
#                     "Mono,aEndo,Mig1","Mono,aEndo,Mig2","Mono,aEndo,Mac","Mono,aEndo,SMC")
# ctlist_to_check = combn(so$ct.major %>% unique,3) %>% .[,which(.[1,] == "Mono")] %>% apply(.,2,function(x) paste0(x, collapse=","))
#' cell type list to check
#' generate customed cell types to check

