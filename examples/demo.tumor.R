## known reference markers for cell types
ref.ctlist = list(Cancer = c("EPCAM","KRT7","KRT18"), #"MKI67"
                  Myeloid = c("CD68","LYZ","AIF1"),
                  T_cell = c("CD3D","CD3E","CD3G"),
                  Mast_cell = c("MS4A2","TPSAB1","CPA3"),
                  Mast = c("MS4A2","TPSAB1","CPA3"),
                  B_cell = c("CD79A","CD79B"),
                  DC = c("CD1C","CD207","CLEC9A","LILRA4","CCL17"),
                  Fibroblast = c("COL1A1","BGN","DCN"),
                  Alveolar = c("CLDN18","SFTPA1","SFTPA2","SFTPC"),
                  EC = c("CLDN5","PECAM1","VWF"),
                  Enteric_glia = c("S100B","PLP1"),
                  Erythroblast = c("HBB","HBA1","HBA2","HBG2"),
                  Epithelial.CRC = c("MT1E","MT1G","ITLN1","ZG16"),
                  Epithelial.LC = c("CAPS", "TPP3"),
                  Pericyte = c("RGS5")) ## different markers of the same cell type from different tissues were markered by ".tissue"


## markers used to define the subclusters in the paper, unaligned and aligned
## For ECs: unaligned: 13, aligned 9
## F
ref.subctlist.unaligned = list(EC=list(
  tipEC_ESM1 = c("ESM1","NID2"),
  vEC_ACKR1 = c("ACKR1","SELP"),
  capEC_CA4 = c("CA4","CD36"),
  artEC = c("FBLN5","GJA5"),
  lymEC = c("PROX1","PDPN"),
  cEC_1 = c("HSPG2") )
)
ref.subctlist.aligned = list()


## Create seurat objects, remove low quality features and cells -----
## pretreatment for Lung cancer
dataDir = "data/E-MTAB-6149/LungCancer/LC_counts/"
dataType = "10x"
metafile = "data/E-MTAB-6149/LungCancer/2097-Lungcancer_metadata.csv"

treatedDir = "data/E-MTAB-6149/LungCancer/LC_treated/"
createDir(treatedDir)

so = load_data(dataDir = dataDir, dataType = dataType, metaPath = metafile)
so = Seurat::UpdateSeuratObject(so)
solist = separate_so(so, batch = "orig.ident")

ct.colorlist = get_random_colorList(so$CellType %>% unique)
saveRDS(ct.colorlist, file.path(treatedDir,"ct.colorlist.rds"))

for (i in 1:length(solist)){
  run_seurat(solist[[i]], label = names(solist)[i], features_to_view = c("CellFromTumor","PatientNumber","TumorType","TumorSite","CellType"),
             save = T, so_dir = treatedDir, treat.cellname = T, mtlabel = "MT-", min_avg = 0.01, assay = "RNA", addNamedColorList = ct.colorlist,
             scale_factor = 1e+06, min_bin = 50, min_cutoff = c(0.01, Inf), dimnum = "auto", dim_sd = 1.25, dim_pval = 0.01, clustRes = 0.8)  
}
## pretreatment for CRC
dataDir = "data/E-MTAB-6149/Colerectalcancer/CRC_counts/"
dataType = "10x"
metafile = "data/E-MTAB-6149/Colerectalcancer/2099-Colorectalcancer_metadata.csv"

treatedDir = "data/E-MTAB-6149/Colerectalcancer/CRC_treated/"
createDir(treatedDir)

so = load_data(dataDir = dataDir, dataType = dataType, metaPath = metafile)
so = Seurat::UpdateSeuratObject(so)
solist = separate_so(so, batch = "orig.ident")

ct.colorlist = c(ct.colorlist, get_random_colorList(setdiff(so$CellType %>% unique, ct.colorlist)))
ct.colorlist = ct.colorlist[!is.na(names(ct.colorlist))]
saveRDS(ct.colorlist, file.path(treatedDir,"ct.colorlist.rds"))

for (i in 1:length(solist)){
  run_seurat(solist[[i]], label = names(solist)[i], features_to_view = c("CellFromTumor","PatientNumber","TumorType","TumorSite","CellType"),
             save = T, so_dir = treatedDir, treat.cellname = T, mtlabel = "MT-", min_avg = 0.01, assay = "RNA", addNamedColorList = ct.colorlist,
             scale_factor = 1e+06, min_bin = 50, min_cutoff = c(0.01, Inf), dimnum = "auto", dim_sd = 1.25, dim_pval = 0.01, clustRes = 0.8)  
}
## pretreatment for OvC
dataDir = "data/E-MTAB-6149/Ovariancancer/OvC_counts/"
dataType = "10x"
metafile = "data/E-MTAB-6149/Ovariancancer/2101-Ovariancancer_metadata.csv"

treatedDir = "data/E-MTAB-6149/Ovariancancer/OvC_treated/"
createDir(treatedDir)

so = load_data(dataDir = dataDir, dataType = dataType, metaPath = metafile)
so = Seurat::UpdateSeuratObject(so)
solist = separate_so(so, batch = "orig.ident")

ct.colorlist = c(ct.colorlist, get_random_colorList(setdiff(so$CellType %>% unique, ct.colorlist)))
saveRDS(ct.colorlist, file.path(treatedDir,"ct.colorlist.rds"))

for (i in 1:length(solist)){
  run_seurat(solist[[i]], label = names(solist)[i], features_to_view = c("CellFromTumor","PatientNumber","TumorType","TumorSite","CellType"),
             save = T, so_dir = treatedDir, treat.cellname = T, mtlabel = "MT-", min_avg = 0.01, assay = "RNA", addNamedColorList = ct.colorlist,
             scale_factor = 1e+06, min_bin = 50, min_cutoff = c(0.01, Inf), dimnum = "auto", dim_sd = 1.25, dim_pval = 0.01, clustRes = 0.8)  
}
## pretreatment for Breast cancer
dataDir = "data/E-MTAB-6149/Breastcancer/BC_counts/"
dataType = "10x"
metafile = "data/E-MTAB-6149/Breastcancer/2103-Breastcancer_metadata.csv"

treatedDir = "data/E-MTAB-6149/Breastcancer/BrC_treated/"
createDir(treatedDir)

meta = load_file(metafile)
so = load_data(dataDir = dataDir, dataType = dataType, metaPath = NULL)
so@meta.data = cbind(so@meta.data, meta)
solist = separate_so(so, batch = "orig.ident")

ct.colorlist = c(ct.colorlist, get_random_colorList(setdiff(so$CellType %>% unique, ct.colorlist)))
saveRDS(ct.colorlist, file.path(treatedDir,"ct.colorlist.rds"))

for (i in 1:length(solist)){
  run_seurat(solist[[i]], label = names(solist)[i], features_to_view = c("CellFromTumor","PatientNumber","TumorType","TumorSite","CellType"),
             save = T, so_dir = treatedDir, treat.cellname = T, mtlabel = "MT-", min_avg = 0.01, assay = "RNA", addNamedColorList = ct.colorlist,
             scale_factor = 1e+06, min_bin = 50, min_cutoff = c(0.01, Inf), dimnum = "auto", dim_sd = 1.25, dim_pval = 0.01, clustRes = 0.8)  
}



## Update seurat objects with given annotation cell types and predicted by seurat clusters ----
## default clustering resolution used 0.8
treatedDirs = named_vec(c("LC","CRC","OvC", "BrC"), 
                        c("data/E-MTAB-6149/LungCancer/LC_treated",
                          "data/E-MTAB-6149/Colerectalcancer/CRC_treated",
                          "data/E-MTAB-6149/Ovariancancer/OvC_treated",
                          "data/E-MTAB-6149/Breastcancer/BrC_treated"))
for (test_tissue in c("LC","CRC","OvC","BrC")) { #NULL, LC
  test_dirs = list.dirs(treatedDirs[test_tissue],recursive = F)
  for (test_dir in test_dirs){
    test = load_so(so_dir = test_dir)
    test_annotation = estimate_pretreated_so_by_refMarkers(so = test, 
                                                           refList = ref.ctlist, 
                                                           tissue = test_tissue, 
                                                           clusterMarkerPath = file.path(test_dir, "t_markers.cluster.csv"),
                                                           annotateName = "CellType",
                                                           clusterName = "seurat_clusters")
    test_annotation = test_annotation %>% dplyr::arrange(as.numeric(cluster))
    ## check mislabel and conflicts manually
    dup_clusts = test_annotation$cluster %>% table %>% .[. > 1] %>% names %>% unique
    if (length(dup_clusts) > 0){
      test$CellTypeStatus = ""
      test$CellTypeM = ""
      test$CellTypeMarker = ""
      for (j in 1:nrow(test_annotation)){
        test@meta.data[which(test$CellType %in% (test_annotation[j, "ct"] %>% unlist %>% strsplit(",") %>% unlist) & test$seurat_clusters == test_annotation[j, "cluster"]),"CellTypeStatus"] = test_annotation$status[j]
        test@meta.data[which(test$CellType %in% (test_annotation[j, "ct"] %>% unlist %>% strsplit(",") %>% unlist) & test$seurat_clusters == test_annotation[j, "cluster"]),"CellTypeM"] = test_annotation$ct[j]
        test@meta.data[which(test$CellType %in% (test_annotation[j, "ct"] %>% unlist %>% strsplit(",") %>% unlist) & test$seurat_clusters == test_annotation[j, "cluster"]),"CellTypeMarker"] = test_annotation$newMark[j]
      } 
    } else {
      test$CellTypeStatus = named_vec(test_annotation$cluster,test_annotation$status)[test$seurat_clusters %>% as.character]
      test$CellTypeM = named_vec(test_annotation$cluster,test_annotation$ct)[test$seurat_clusters %>% as.character]
      test$CellTypeMarker = named_vec(test_annotation$cluster,test_annotation$newMark)[test$seurat_clusters %>% as.character]
    }
    update_so(so = test, so_dir = test_dir, meta.only = T)
  }
}

##################################################################################
##################################################################################
################################ Main Analysis 1 #################################
##################################################################################
require(R.filesets)
require(jzhang.utils)
require(scCellMediator)
require(dplyr)
`%+%` <- paste0 
setwd("~/Documents/workspace/softwares/iRPackages/scCellMediator/")
treatedDirs = named_vec(c("LC","CRC","OvC", "BrC"), 
                        c("data/E-MTAB-6149/Lungcancer/LC_treated",
                          "data/E-MTAB-6149/Colerectalcancer/CRC_treated",
                          "data/E-MTAB-6149/Ovariancancer/OvC_treated",
                          "data/E-MTAB-6149/Breastcancer/BrC_treated"))
for (test_tissue in c("LC","CRC","OvC","BrC")) { #NULL, LC test_tissue = "CRC"
  test_dir = treatedDirs[test_tissue]
  run_combined_seurat(test_dir, 
                      save_dir = test_dir %+% "_integrated/", 
                      ifReturn = FALSE, 
                      filterRules = 'CellType %in% c("EC","T_cell","Cancer") & stringr::str_detect(.$CellTypeStatus, "CT_confirmed") & nCount_RNA > 500 & nFeature_RNA > 500', 
                      acceptRules = 'CellType,>,50', 
                      dimnum = 50, k.anchor = 5, k.filter = 50, k.weight = 20)
  ## accepted samples for CRC:
  ## c("filtered_scrEXT001","filtered_scrEXT002","filtered_scrEXT009","filtered_scrEXT0012","filtered_scrEXT0013","filtered_scrEXT0018","filtered_scrEXT0019",
  ##   "filtered_scrEXT0020","filtered_scrEXT0023","filtered_scrEXT0025","filtered_scrEXT0027","filtered_scrEXT0029")
}
## analysis for CRC ## ----
## ONLY check Cancer, T_cell and EC
test_tissue = "CRC"
test_dir = treatedDirs[test_tissue]
merged_dir = file.path(test_dir %+% "_integrated")
merged = load_so(so_dir = merged_dir)
## Generate UMAP plots for merged seurat object
features_to_view = c("CellFromTumor","PatientNumber","TumorSite","CellType","integrated_snn_res.1","CellTypeStatus","CellTypeM","CellTypeMarker")
addNamedColorList = loadRDS(file.path(test_dir,"ct.colorlist.rds"))
for (ifeature in features_to_view) {
  print(ifeature)
  colorlist = get_random_colorlist(merged@meta.data[,ifeature] %>% unique)
  labelon = F
  fnum = uniqlen(merged@meta.data[,ifeature])
  if (fnum > 1 & fnum < 80) {
    if (fnum < 20){ labelon = T}
    if (ifeature %in% c("CellType","CellTypeM")){
      colorlist = c(addNamedColorList, colorlist)
    }
    p = UMAPPlot(merged, group.by = ifeature, pt.size = 1.2, label = labelon) +
      ggplot2::scale_color_manual(values = colorlist) +
      NoAxes()
    # p = AugmentPlot(p, width = 6.5, height = 5, dpi = 100)
    if (ifeature %in% c("CellTypeMarker")){
      savePlot(pObj = p, pType = "umap", pName = ifeature, pDir = merged_dir, pFmt = "pdf", pWd = 12, pHt = 5)
    } else {
      savePlot(pObj = p, pType = "umap", pName = ifeature, pDir = merged_dir, pFmt = "pdf", pWd = 6.5, pHt = 5)
    }
  }
}
## Total 6991 cells accepted
## Separate cancer-originated (6276 cells) and normal-originated (715 cells)
merged_so4cellchat = separate_so(merged, batch = "CellFromTumor")
names(merged_so4cellchat) = c("CRC","CTRL")
### For CRC
so_label = "CRC"
so = merged_so4cellchat[[so_label]]
so$ct.major <- so$CellType
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
labels <- so$ct.major
meta <- data.frame(ct = labels, row.names = names(labels), sample = so$orig.ident)
#### run cell chat analysis
cellchatDir = file.path(merged_dir,"tmp." %+% so_label %+% ".cellchat.out")
createDir(cellchatDir)
run_cellchat(so, feature.name = "ct.major", batch.name = "orig.ident",
             output = cellchatDir, 
             minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10,
             suggestedRunTimes = F, species = "human")
#### run causal mediation analysis for Endothelial cells
ctlist_to_check = generate_triple_ct_model(so$ct.major %>% unique, cells_to_fix = named_vec(c(3),c("T_cell")), keep_or_remove="keep")
#### this is to avoid investigate the relationships between two subtypes from the same major cell types (like explore the relationships between aEndo, vEndo and capEndo)
for (cts in ctlist_to_check){
  endo_sum = stringr::str_match_all(cts,"EC") %>% unlist %>% length
  if (endo_sum > 1){
    ctlist_to_check = ctlist_to_check[!(ctlist_to_check == cts)]
  }
}
#### start cma analysis ##
inpath = cellchatDir
mat4cma_filepath = file.path(dirname(cellchatDir),"tmp." %+% so_label %+% ".CausalMediationOutput","mat4cma.Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
#### check functions to T cell exhaustion
get_functions_to_view(c("exhaust","t"), mod = "inspect", overlapKeywords = TRUE, species = "human", cat = "C7") %>% unique # , cat = "C5"
func_list = get_functions_to_view(c("GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP"), mod = "catch", overlapKeywords = FALSE, species = "human", cat = "C7")
funcToUse = func_list
#### generate matrix for causal mediation analysis ## the function is not available for multiple function test yet.
DefaultAssay(so) = "RNA"
mat4cma = generate_matrix_for_CMA(so,
                                  dataDir = inpath,
                                  feature.name = "ct.major",
                                  ctlist_to_check = ctlist_to_check,
                                  funcToUse = funcToUse,
                                  species = "human",
                                  batch = "orig.ident",
                                  ifSave = TRUE,
                                  ifReturn = TRUE,
                                  outputFile = mat4cma_filepath)
#### generate causal mediation relationships
pvalue = 0.05
cma_path = file.path(dirname(dataDir),"tmp." %+% so_label %+% ".CausalMediationOutput","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
# mat4cma = loadRDS(mat4cma_filepath)
run_CausalMediation(mat4mca, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "CRC_CancerECTcell_TcellExhaustion", saveDir = cma_path) #mod = "spc",
### For normal
so_label = "CTRL"
so = merged_so4cellchat[[so_label]]
so$ct.major <- so$CellType
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
labels <- so$ct.major
meta <- data.frame(ct = labels, row.names = names(labels), sample = so$orig.ident)
#### run cell chat analysis
cellchatDir = file.path(merged_dir,"tmp." %+% so_label %+% ".cellchat.out")
createDir(cellchatDir)
run_cellchat(so, feature.name = "ct.major", batch.name = "orig.ident",
             output = cellchatDir, 
             minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10,
             suggestedRunTimes = F, species = "human")
#### run causal mediation analysis for Endothelial cells
ctlist_to_check = generate_triple_ct_model(so$ct.major %>% unique, cells_to_fix = named_vec(c(3),c("T_cell")), keep_or_remove="keep")
#### this is to avoid investigate the relationships between two subtypes from the same major cell types (like explore the relationships between aEndo, vEndo and capEndo)
for (cts in ctlist_to_check){
  endo_sum = stringr::str_match_all(cts,"EC") %>% unlist %>% length
  if (endo_sum > 1){
    ctlist_to_check = ctlist_to_check[!(ctlist_to_check == cts)]
  }
}
#### start cma analysis ##
inpath = cellchatDir
mat4cma_filepath = file.path(dirname(cellchatDir),"tmp.CausalMediationOutput","mat4cma.Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
#### check functions to T cell exhaustion
get_functions_to_view(c("exhaust","t"), mod = "inspect", overlapKeywords = TRUE, species = "human", cat = "C7") %>% unique # , cat = "C5"
func_list = get_functions_to_view(c("GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP"), mod = "catch", overlapKeywords = FALSE, species = "human", cat = "C7")
funcToUse = func_list
#### generate matrix for causal mediation analysis ## the function is not available for multiple function test yet.
DefaultAssay(so) = "RNA"
mat4cma = generate_matrix_for_CMA(so,
                                  dataDir = inpath,
                                  feature.name = "ct.major",
                                  ctlist_to_check = ctlist_to_check,
                                  funcToUse = funcToUse,
                                  species = "human",
                                  batch = "orig.ident",
                                  ifSave = TRUE,
                                  ifReturn = TRUE,
                                  outputFile = mat4cma_filepath)
#### generate causal mediation relationships
pvalue = 0.05
cma_path = file.path(dirname(dataDir),"tmp.CausalMediationOutput","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
# mat4cma = loadRDS(mat4cma_filepath)
run_CausalMediation(mat4mca, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "CRC_CancerECTcell_TcellExhaustion", saveDir = cma_path) #mod = "spc",



##################################################################################
##################################################################################
################################ Main Analysis 2 #################################
##################################################################################
require(R.filesets)
require(jzhang.utils)
require(scCellMediator)
require(dplyr)
`%+%` <- paste0 
setwd("~/Documents/workspace/softwares/iRPackages/scCellMediator/")
treatedDirs = named_vec(c("LC","CRC","OvC", "BrC"), 
                        c("data/E-MTAB-6149/Lungcancer/LC_treated",
                          "data/E-MTAB-6149/Colerectalcancer/CRC_treated",
                          "data/E-MTAB-6149/Ovariancancer/OvC_treated",
                          "data/E-MTAB-6149/Breastcancer/BrC_treated"))
test_tissue = "CRC"
test_dir = treatedDirs[test_tissue]
run_combined_seurat(test_dir, 
                    save_dir = test_dir %+% "_integrated_ECMyeloidTcell/", 
                    ifReturn = FALSE, 
                    filterRules = 'CellType %in% c("EC","T_cell","Myeloid") & stringr::str_detect(.$CellTypeStatus, "CT_confirmed") & nCount_RNA > 500 & nFeature_RNA > 500', 
                    acceptRules = 'CellType,>,50', 
                    dimnum = 50, k.anchor = 5, k.filter = 50, k.weight = 20)
## accepted samples for CRC:
## c("filtered_scrEXT001","filtered_scrEXT002","filtered_scrEXT003","filtered_scrEXT009","filtered_scrEXT0010",
##   "filtered_scrEXT0011","filtered_scrEXT0012","filtered_scrEXT0013","filtered_scrEXT0018","filtered_scrEXT0019",
##   "filtered_scrEXT0020","filtered_scrEXT0021","filtered_scrEXT0025","filtered_scrEXT0027","filtered_scrEXT0029")
## analysis for CRC ## ----
## ONLY check Myeloid, T_cell and EC
test_tissue = "CRC"
test_dir = treatedDirs[test_tissue]
merged_dir = file.path(test_dir %+% "_integrated_ECMyeloidTcell/")
merged = load_so(so_dir = merged_dir)
## Purify cell types
markers_for_hsNeutrophil = c("ITGAM","FCGR3A","PTPRC","LY6H","FUT4") # ITGAM -> CD11B, "FCGR3A" -> CD16, PTPRC -> CD45, LY6H -> Ly6C, FUT4 -> CD15
p = DotPlot(merged, features = markers_for_hsNeutrophil, assay = "RNA", group.by = "integrated_snn_res.1") 
savePlot(pObj = p, pType = "dotplot", pName = "neutrophil_markers", pDir = merged_dir, pFmt = "pdf", pWd = 6.5, pHt = 4.5)
clusters_for_hsNeutrophil = p$data %>% filter(avg.exp > 0, pct.exp > 0, features.plot %in% c("ITGAM", "FCGR3A"), avg.exp.scaled >= 1) %>% pull(id) %>% unique %>% as.character %>% as.numeric
merged$CellType.old = merged$CellType
merged$CellType[merged$seurat_clusters %in% clusters_for_hsNeutrophil] = "Neutrophil"
cells_filtered_based_on_umap = c(Cells(merged)[merged@reductions$umap@cell.embeddings[,1] < -5 & merged$CellType %in% c("Neutrophil")],
                                 Cells(merged)[merged@reductions$umap@cell.embeddings[,1] > -5 & merged@reductions$umap@cell.embeddings[,2] > 5 & merged$CellType %in% c("EC")],
                                 Cells(merged)[merged@reductions$umap@cell.embeddings[,1] > -5 & merged@reductions$umap@cell.embeddings[,2] < 5 & merged$CellType %in% c("T_cell")])
merged = subset(merged, cells = cells_filtered_based_on_umap)
p = UMAPPlot(merged, group.by = "CellType", pt.size = 1.2, label = TRUE) +
  ggplot2::scale_color_manual(values=c(addNamedColorList,named_vec("Neutrophil","blue"))) +
  NoAxes()
savePlot(pObj = p, pType = "umap", pName = "new_celltype", pDir = merged_dir, pFmt = "pdf", pWd = 6.5, pHt = 5)
p = UMAPPlot(merged, group.by = "CellFromTumor", pt.size = 1.2, split.by = "CellFromTumor") +
  ggsci::scale_color_aaas() +
  NoAxes()
savePlot(pObj = p, pType = "umap", pName = "split_cellFromTumor", pDir = merged_dir, pFmt = "pdf", pWd = 8.4, pHt = 4)
update_so(merged, so_dir = merged_dir, meta.only = F)
## Generate UMAP plots for merged seurat object
features_to_view = c("CellFromTumor","PatientNumber","TumorSite","CellType","integrated_snn_res.1","CellTypeStatus","CellTypeM","CellTypeMarker")
addNamedColorList = loadRDS(file.path(test_dir,"ct.colorlist.rds"))
for (ifeature in features_to_view) {
  print(ifeature)
  colorlist = get_random_colorlist(merged@meta.data[,ifeature] %>% unique)
  labelon = F
  fnum = uniqlen(merged@meta.data[,ifeature])
  if (fnum > 1 & fnum < 80) {
    if (fnum < 20){ labelon = T}
    if (ifeature %in% c("CellType","CellTypeM")){
      colorlist = c(addNamedColorList, colorlist)
    }
    p = UMAPPlot(merged, group.by = ifeature, pt.size = 1.2, label = labelon) +
      ggplot2::scale_color_manual(values = colorlist) +
      NoAxes()
    # p = AugmentPlot(p, width = 6.5, height = 5, dpi = 100)
    if (ifeature %in% c("CellTypeMarker")){
      savePlot(pObj = p, pType = "umap", pName = ifeature, pDir = merged_dir, pFmt = "pdf", pWd = 12, pHt = 5)
    } else {
      savePlot(pObj = p, pType = "umap", pName = ifeature, pDir = merged_dir, pFmt = "pdf", pWd = 6.5, pHt = 5)
    }
  }
}
## Total 8451 cells accepted
## Separate cancer-originated (6678 cells) and normal-originated (1773 cells)
merged_so4cellchat = separate_so(merged, batch = "CellFromTumor")
names(merged_so4cellchat) = c("CRC","CTRL")
### For CRC
so_label = "CRC"
so = merged_so4cellchat[[so_label]]
so$ct.major <- so$CellType
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
labels <- so$ct.major
meta <- data.frame(ct = labels, row.names = names(labels), sample = so$orig.ident)
#### run cell chat analysis
cellchatDir = file.path(merged_dir,"tmp." %+% so_label %+% ".cellchat.out")
createDir(cellchatDir)
run_cellchat(so, feature.name = "ct.major", batch.name = "orig.ident",
             output = cellchatDir, 
             minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10,
             suggestedRunTimes = F, species = "human")
#### run causal mediation analysis for Endothelial cells
ctlist_to_check = generate_triple_ct_model(so$ct.major %>% unique, cells_to_fix = named_vec(c(3),c("T_cell")), keep_or_remove="keep")
#### this is to avoid investigate the relationships between two subtypes from the same major cell types (like explore the relationships between aEndo, vEndo and capEndo)
for (cts in ctlist_to_check){
  endo_sum = stringr::str_match_all(cts,"EC") %>% unlist %>% length
  if (endo_sum > 1){
    ctlist_to_check = ctlist_to_check[!(ctlist_to_check == cts)]
  }
}
#### start cma analysis ##
inpath = cellchatDir
mat4cma_filepath = file.path(dirname(cellchatDir),"tmp." %+% so_label %+% ".CausalMediationOutput","mat4cma.Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
#### check functions to T cell exhaustion
get_functions_to_view(c("exhaust","t"), mod = "inspect", overlapKeywords = TRUE, species = "human", cat = "C7") %>% unique # , cat = "C5"
func_list = get_functions_to_view(c("GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP"), mod = "catch", overlapKeywords = FALSE, species = "human", cat = "C7")
funcToUse = func_list
#### generate matrix for causal mediation analysis ## the function is not available for multiple function test yet.
DefaultAssay(so) = "RNA"
mat4cma = generate_matrix_for_CMA(so,
                                  dataDir = inpath,
                                  feature.name = "ct.major",
                                  ctlist_to_check = ctlist_to_check,
                                  funcToUse = funcToUse,
                                  species = "human",
                                  batch = "orig.ident",
                                  ifSave = TRUE,
                                  ifReturn = TRUE,
                                  outputFile = mat4cma_filepath)
#### generate causal mediation relationships
pvalue = 0.05
cma_path = file.path(dirname(dataDir),"tmp." %+% so_label %+% ".CausalMediationOutput","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
# mat4cma = loadRDS(mat4cma_filepath)
run_CausalMediation(mat4mca, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "CRC_CancerECTcell_TcellExhaustion", saveDir = cma_path) #mod = "spc",
### For normal
so_label = "CTRL"
so = merged_so4cellchat[[so_label]]
so$ct.major <- so$CellType
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
labels <- so$ct.major
meta <- data.frame(ct = labels, row.names = names(labels), sample = so$orig.ident)
#### run cell chat analysis
cellchatDir = file.path(merged_dir,"tmp." %+% so_label %+% ".cellchat.out")
createDir(cellchatDir)
run_cellchat(so, feature.name = "ct.major", batch.name = "orig.ident",
             output = cellchatDir, 
             minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10,
             suggestedRunTimes = F, species = "human")
#### run causal mediation analysis for Endothelial cells
ctlist_to_check = generate_triple_ct_model(so$ct.major %>% unique, cells_to_fix = named_vec(c(3),c("T_cell")), keep_or_remove="keep")
#### this is to avoid investigate the relationships between two subtypes from the same major cell types (like explore the relationships between aEndo, vEndo and capEndo)
for (cts in ctlist_to_check){
  endo_sum = stringr::str_match_all(cts,"EC") %>% unlist %>% length
  if (endo_sum > 1){
    ctlist_to_check = ctlist_to_check[!(ctlist_to_check == cts)]
  }
}
#### start cma analysis ##
inpath = cellchatDir
mat4cma_filepath = file.path(dirname(cellchatDir),"tmp." %+% so_label %+% ".CausalMediationOutput","mat4cma.Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
#### check functions to T cell exhaustion
get_functions_to_view(c("exhaust","t"), mod = "inspect", overlapKeywords = TRUE, species = "human", cat = "C7") %>% unique # , cat = "C5"
func_list = get_functions_to_view(c("GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP"), mod = "catch", overlapKeywords = FALSE, species = "human", cat = "C7")
funcToUse = func_list
#### generate matrix for causal mediation analysis ## the function is not available for multiple function test yet.
DefaultAssay(so) = "RNA"
mat4cma = generate_matrix_for_CMA(so,
                                  dataDir = inpath,
                                  feature.name = "ct.major",
                                  ctlist_to_check = ctlist_to_check,
                                  funcToUse = funcToUse,
                                  species = "human",
                                  batch = "orig.ident",
                                  ifSave = TRUE,
                                  ifReturn = TRUE,
                                  outputFile = mat4cma_filepath)
#### generate causal mediation relationships
pvalue = 0.05
cma_path = file.path(dirname(dataDir),"tmp." %+% so_label %+% ".CausalMediationOutput","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
# mat4cma = loadRDS(mat4cma_filepath)
run_CausalMediation(mat4mca, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "CRC_CancerECTcell_TcellExhaustion", saveDir = cma_path) #mod = "spc",



##################################################################################
##################################################################################
################################ Main Analysis 3 #################################
##################################################################################
require(R.filesets)
require(jzhang.utils)
require(scCellMediator)
require(dplyr)
`%+%` <- paste0 
setwd("~/Documents/workspace/softwares/iRPackages/scCellMediator/")
treatedDirs = named_vec(c("LC","CRC","OvC", "BrC"), 
                        c("data/E-MTAB-6149/Lungcancer/LC_treated",
                          "data/E-MTAB-6149/Colerectalcancer/CRC_treated",
                          "data/E-MTAB-6149/Ovariancancer/OvC_treated",
                          "data/E-MTAB-6149/Breastcancer/BrC_treated"))
test_tissue = "LC"
test_dir = treatedDirs[test_tissue]
run_combined_seurat(test_dir, 
                    save_dir = test_dir %+% "_integrated_ECTcellCancer/", 
                    ifReturn = FALSE, 
                    filterRules = 'CellType %in% c("EC","T_cell","Cancer") & stringr::str_detect(.$CellTypeStatus, "CT_confirmed") & nCount_RNA > 500 & nFeature_RNA > 500', 
                    acceptRules = 'CellType,>,50', 
                    dimnum = 50, k.anchor = 5, k.filter = 50, k.weight = 20)
## accepted samples for LC:
## list.dirs(test_dir %+% "_integrated_ECTcellCancer/")
## analysis for LC ## ----
## ONLY check Cancer, T_cell and EC
test_tissue = "LC"
test_dir = treatedDirs[test_tissue]
merged_dir = file.path(test_dir %+% "_integrated_ECTcellCancer")
merged = load_so(so_dir = merged_dir)
## Generate UMAP plots for merged seurat object
features_to_view = c("CellFromTumor","PatientNumber","TumorSite","CellType","integrated_snn_res.1","CellTypeStatus","CellTypeM","CellTypeMarker")
addNamedColorList = loadRDS(file.path(test_dir,"ct.colorlist.rds"))
for (ifeature in features_to_view) {
  print(ifeature)
  colorlist = get_random_colorlist(merged@meta.data[,ifeature] %>% unique)
  labelon = F
  fnum = uniqlen(merged@meta.data[,ifeature])
  if (fnum > 1 & fnum < 80) {
    if (fnum < 20){ labelon = T}
    if (ifeature %in% c("CellType","CellTypeM")){
      colorlist = c(addNamedColorList, colorlist)
    }
    p = UMAPPlot(merged, group.by = ifeature, pt.size = 1.2, label = labelon) +
      ggplot2::scale_color_manual(values = colorlist) +
      NoAxes()
    # p = AugmentPlot(p, width = 6.5, height = 5, dpi = 100)
    if (ifeature %in% c("CellTypeMarker")){
      savePlot(pObj = p, pType = "umap", pName = ifeature, pDir = merged_dir, pFmt = "pdf", pWd = 12, pHt = 5)
    } else {
      savePlot(pObj = p, pType = "umap", pName = ifeature, pDir = merged_dir, pFmt = "pdf", pWd = 6.5, pHt = 5)
    }
  }
}
## Total 10764 cells accepted
## Separate cancer-originated (9020 cells) and normal-originated (1744 cells)
merged_so4cellchat = separate_so(merged, batch = "CellFromTumor")
names(merged_so4cellchat) = c("LC","CTRL")
### For LC
so_label = "LC"
so = merged_so4cellchat[[so_label]]
so$ct.major <- so$CellType
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
labels <- so$ct.major
meta <- data.frame(ct = labels, row.names = names(labels), sample = so$orig.ident)
#### run cell chat analysis
cellchatDir = file.path(merged_dir,"tmp." %+% so_label %+% ".cellchat.out")
createDir(cellchatDir)
run_cellchat(so, feature.name = "ct.major", batch.name = "orig.ident",
             output = cellchatDir, 
             minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10,
             suggestedRunTimes = F, species = "human")
#### run causal mediation analysis for Endothelial cells
ctlist_to_check = generate_triple_ct_model(so$ct.major %>% unique, cells_to_fix = named_vec(c(3),c("T_cell")), keep_or_remove="keep")
#### this is to avoid investigate the relationships between two subtypes from the same major cell types (like explore the relationships between aEndo, vEndo and capEndo)
for (cts in ctlist_to_check){
  endo_sum = stringr::str_match_all(cts,"EC") %>% unlist %>% length
  if (endo_sum > 1){
    ctlist_to_check = ctlist_to_check[!(ctlist_to_check == cts)]
  }
}
#### start cma analysis ##
inpath = cellchatDir
mat4cma_filepath = file.path(dirname(cellchatDir),"tmp." %+% so_label %+% ".CausalMediationOutput","mat4cma.Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
#### check functions to T cell exhaustion
get_functions_to_view(c("exhaust","t"), mod = "inspect", overlapKeywords = TRUE, species = "human", cat = "C7") %>% unique # , cat = "C5"
func_list = get_functions_to_view(c("GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP"), mod = "catch", overlapKeywords = FALSE, species = "human", cat = "C7")
funcToUse = func_list
#### generate matrix for causal mediation analysis ## the function is not available for multiple function test yet.
DefaultAssay(so) = "RNA"
mat4cma = generate_matrix_for_CMA(so,
                                  dataDir = inpath,
                                  feature.name = "ct.major",
                                  ctlist_to_check = ctlist_to_check,
                                  funcToUse = funcToUse,
                                  species = "human",
                                  batch = "orig.ident",
                                  ifSave = TRUE,
                                  ifReturn = TRUE,
                                  outputFile = mat4cma_filepath)
#### generate causal mediation relationships
pvalue = 0.05
cma_path = file.path(dirname(cellchatDir),"tmp." %+% so_label %+% ".CausalMediationOutput","mat4cma.Tcell_receiver.rds")
createDir(cma_path)
# mat4cma = loadRDS(mat4cma_filepath)
run_CausalMediation(mat4mca, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "LC_ECMyeloidTcell_TcellExhaustion", saveDir = dirname(cma_path)) #mod = "spc",
### For normal
so_label = "CTRL"
so = merged_so4cellchat[[so_label]]
so$ct.major <- so$CellType
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
labels <- so$ct.major
meta <- data.frame(ct = labels, row.names = names(labels), sample = so$orig.ident)
#### run cell chat analysis
cellchatDir = file.path(merged_dir,"tmp." %+% so_label %+% ".cellchat.out")
createDir(cellchatDir)
run_cellchat(so, feature.name = "ct.major", batch.name = "orig.ident",
             output = cellchatDir, 
             minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10,
             suggestedRunTimes = F, species = "human")
#### run causal mediation analysis for Endothelial cells
ctlist_to_check = generate_triple_ct_model(so$ct.major %>% unique, cells_to_fix = named_vec(c(3),c("T_cell")), keep_or_remove="keep")
#### this is to avoid investigate the relationships between two subtypes from the same major cell types (like explore the relationships between aEndo, vEndo and capEndo)
for (cts in ctlist_to_check){
  endo_sum = stringr::str_match_all(cts,"EC") %>% unlist %>% length
  if (endo_sum > 1){
    ctlist_to_check = ctlist_to_check[!(ctlist_to_check == cts)]
  }
}
#### start cma analysis ##
inpath = cellchatDir
mat4cma_filepath = file.path(dirname(cellchatDir),"tmp.CausalMediationOutput","mat4cma.Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
#### check functions to T cell exhaustion
get_functions_to_view(c("exhaust","t"), mod = "inspect", overlapKeywords = TRUE, species = "human", cat = "C7") %>% unique # , cat = "C5"
func_list = get_functions_to_view(c("GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP"), mod = "catch", overlapKeywords = FALSE, species = "human", cat = "C7")
funcToUse = func_list
#### generate matrix for causal mediation analysis ## the function is not available for multiple function test yet.
DefaultAssay(so) = "RNA"
mat4cma = generate_matrix_for_CMA(so,
                                  dataDir = inpath,
                                  feature.name = "ct.major",
                                  ctlist_to_check = ctlist_to_check,
                                  funcToUse = funcToUse,
                                  species = "human",
                                  batch = "orig.ident",
                                  ifSave = TRUE,
                                  ifReturn = TRUE,
                                  outputFile = mat4cma_filepath)
#### generate causal mediation relationships
pvalue = 0.05
cma_path = file.path(dirname(dataDir),"tmp.CausalMediationOutput","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
# mat4cma = loadRDS(mat4cma_filepath)
run_CausalMediation(mat4mca, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "CRC_CancerECTcell_TcellExhaustion", saveDir = cma_path) #mod = "spc",



##################################################################################
##################################################################################
################################ Main Analysis 4 #################################
##################################################################################
require(R.filesets)
require(jzhang.utils)
require(scCellMediator)
require(dplyr)
`%+%` <- paste0 
setwd("~/Documents/workspace/softwares/iRPackages/scCellMediator/")
treatedDirs = named_vec(c("LC","CRC","OvC", "BrC"), 
                        c("data/E-MTAB-6149/Lungcancer/LC_treated",
                          "data/E-MTAB-6149/Colerectalcancer/CRC_treated",
                          "data/E-MTAB-6149/Ovariancancer/OvC_treated",
                          "data/E-MTAB-6149/Breastcancer/BrC_treated"))
test_tissue = "LC"
test_dir = treatedDirs[test_tissue]
run_combined_seurat(test_dir, 
                    save_dir = test_dir %+% "_integrated_ECMyeloidTcell/", 
                    ifReturn = FALSE, 
                    filterRules = 'CellType %in% c("EC","T_cell","Myeloid") & stringr::str_detect(.$CellTypeStatus, "CT_confirmed") & nCount_RNA > 500 & nFeature_RNA > 500', 
                    acceptRules = 'CellType,>,50', 
                    dimnum = 50, k.anchor = 5, k.filter = 50, k.weight = 20)
## accepted samples for LC:
## c("filtered_scrEXT001","filtered_scrEXT002","filtered_scrEXT003","filtered_scrEXT009","filtered_scrEXT0010",
##   "filtered_scrEXT0011","filtered_scrEXT0012","filtered_scrEXT0013","filtered_scrEXT0018","filtered_scrEXT0019",
##   "filtered_scrEXT0020","filtered_scrEXT0021","filtered_scrEXT0025","filtered_scrEXT0027","filtered_scrEXT0029")
## analysis for CRC ## ----
## ONLY check Myeloid, T_cell and EC
test_tissue = "LC"
test_dir = treatedDirs[test_tissue]
merged_dir = file.path(test_dir %+% "_integrated_ECMyeloidTcell/")
merged = load_so(so_dir = merged_dir)
## Purify cell types
markers_for_hsNeutrophil = c("ITGAM","FCGR3A","PTPRC","LY6H","FUT4") # ITGAM -> CD11B, "FCGR3A" -> CD16, PTPRC -> CD45, LY6H -> Ly6C, FUT4 -> CD15
p = DotPlot(merged, features = markers_for_hsNeutrophil, assay = "RNA", group.by = "integrated_snn_res.1") 
savePlot(pObj = p, pType = "dotplot", pName = "neutrophil_markers", pDir = merged_dir, pFmt = "pdf", pWd = 6.5, pHt = 4.5)
clusters_for_hsNeutrophil = p$data %>% filter(avg.exp > 0, pct.exp > 0, features.plot %in% c("ITGAM", "FCGR3A"), avg.exp.scaled >= 1) %>% pull(id) %>% unique %>% as.character %>% as.numeric
merged$CellType.old = merged$CellType
merged$CellType[merged$seurat_clusters %in% clusters_for_hsNeutrophil] = "Neutrophil"
## Need to pay attention to the order of the cell types
cells_filtered_based_on_umap = c(Cells(merged)[merged@reductions$umap@cell.embeddings[,1] < -5 & merged$CellType %in% c("EC")],
                                 Cells(merged)[merged@reductions$umap@cell.embeddings[,1] > -5 & merged@reductions$umap@cell.embeddings[,2] > 0 & merged$CellType %in% c("T_cell")],
                                 Cells(merged)[merged@reductions$umap@cell.embeddings[,1] > -5 & merged@reductions$umap@cell.embeddings[,2] < 0 & merged$CellType %in% c("Neutrophil")])
merged = subset(merged, cells = cells_filtered_based_on_umap)
p = UMAPPlot(merged, group.by = "CellType", pt.size = 1.2, label = TRUE) +
  ggplot2::scale_color_manual(values=c(addNamedColorList,named_vec("Neutrophil","blue"))) +
  NoAxes()
savePlot(pObj = p, pType = "umap", pName = "new_celltype", pDir = merged_dir, pFmt = "pdf", pWd = 6.5, pHt = 5)
p = UMAPPlot(merged, group.by = "CellFromTumor", pt.size = 1.2, split.by = "CellFromTumor") +
  ggsci::scale_color_aaas() +
  NoAxes()
savePlot(pObj = p, pType = "umap", pName = "split_cellFromTumor", pDir = merged_dir, pFmt = "pdf", pWd = 8.4, pHt = 4)
update_so(merged, so_dir = merged_dir, meta.only = F)
## Generate UMAP plots for merged seurat object
features_to_view = c("CellFromTumor","PatientNumber","TumorSite","CellType","integrated_snn_res.1","CellTypeStatus","CellTypeM","CellTypeMarker")
addNamedColorList = loadRDS(file.path(test_dir,"ct.colorlist.rds"))
for (ifeature in features_to_view) {
  print(ifeature)
  colorlist = get_random_colorlist(merged@meta.data[,ifeature] %>% unique)
  labelon = F
  fnum = uniqlen(merged@meta.data[,ifeature])
  if (fnum > 1 & fnum < 80) {
    if (fnum < 20){ labelon = T}
    if (ifeature %in% c("CellType","CellTypeM")){
      colorlist = c(addNamedColorList, colorlist)
    }
    p = UMAPPlot(merged, group.by = ifeature, pt.size = 1.2, label = labelon) +
      ggplot2::scale_color_manual(values = colorlist) +
      NoAxes()
    # p = AugmentPlot(p, width = 6.5, height = 5, dpi = 100)
    if (ifeature %in% c("CellTypeMarker")){
      savePlot(pObj = p, pType = "umap", pName = ifeature, pDir = merged_dir, pFmt = "pdf", pWd = 12, pHt = 5)
    } else {
      savePlot(pObj = p, pType = "umap", pName = ifeature, pDir = merged_dir, pFmt = "pdf", pWd = 6.5, pHt = 5)
    }
  }
}
## Total 15943 cells accepted
## Separate cancer-originated (10308 cells) and normal-originated (5635 cells)
merged_so4cellchat = separate_so(merged, batch = "CellFromTumor")
names(merged_so4cellchat) = c("LC","CTRL")
### For LC
so_label = "LC"
so = merged_so4cellchat[[so_label]]
so$ct.major <- so$CellType
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
labels <- so$ct.major
meta <- data.frame(ct = labels, row.names = names(labels), sample = so$orig.ident)
#### run cell chat analysis
cellchatDir = file.path(merged_dir,"tmp." %+% so_label %+% ".cellchat.out")
createDir(cellchatDir)
run_cellchat(so, feature.name = "ct.major", batch.name = "orig.ident",
             output = cellchatDir, 
             minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10,
             suggestedRunTimes = F, species = "human")
#### run causal mediation analysis for Endothelial cells
ctlist_to_check = generate_triple_ct_model(so$ct.major %>% unique, cells_to_fix = named_vec(c(3),c("T_cell")), keep_or_remove="keep")
#### this is to avoid investigate the relationships between two subtypes from the same major cell types (like explore the relationships between aEndo, vEndo and capEndo)
for (cts in ctlist_to_check){
  endo_sum = stringr::str_match_all(cts,"EC") %>% unlist %>% length
  if (endo_sum > 1){
    ctlist_to_check = ctlist_to_check[!(ctlist_to_check == cts)]
  }
}
#### start cma analysis ##
inpath = cellchatDir
mat4cma_filepath = file.path(dirname(cellchatDir),"tmp." %+% so_label %+% ".CausalMediationOutput","mat4cma.Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
#### check functions to T cell exhaustion
get_functions_to_view(c("exhaust","t"), mod = "inspect", overlapKeywords = TRUE, species = "human", cat = "C7") %>% unique # , cat = "C5"
func_list = get_functions_to_view(c("GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP"), mod = "catch", overlapKeywords = FALSE, species = "human", cat = "C7")
funcToUse = func_list
#### generate matrix for causal mediation analysis ## the function is not available for multiple function test yet.
DefaultAssay(so) = "RNA"
mat4cma = generate_matrix_for_CMA(so,
                                  dataDir = inpath,
                                  feature.name = "ct.major",
                                  ctlist_to_check = ctlist_to_check,
                                  funcToUse = funcToUse,
                                  species = "human",
                                  batch = "orig.ident",
                                  ifSave = TRUE,
                                  ifReturn = TRUE,
                                  outputFile = mat4cma_filepath)
#### generate causal mediation relationships
pvalue = 0.05
cma_path = file.path(dirname(cellchatDir),"tmp." %+% so_label %+% ".CausalMediationOutput","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
# mat4cma = loadRDS(mat4cma_filepath)
run_CausalMediation(mat4mca, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "LC_ECNeutrophilTcell_TcellExhaustion", saveDir = cma_path) #mod = "spc",
### For normal
so_label = "CTRL"
so = merged_so4cellchat[[so_label]]
so$ct.major <- so$CellType
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
labels <- so$ct.major
meta <- data.frame(ct = labels, row.names = names(labels), sample = so$orig.ident)
#### run cell chat analysis
cellchatDir = file.path(merged_dir,"tmp." %+% so_label %+% ".cellchat.out")
createDir(cellchatDir)
run_cellchat(so, feature.name = "ct.major", batch.name = "orig.ident",
             output = cellchatDir, 
             minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10,
             suggestedRunTimes = F, species = "human")
#### run causal mediation analysis for Endothelial cells
ctlist_to_check = generate_triple_ct_model(so$ct.major %>% unique, cells_to_fix = named_vec(c(3),c("T_cell")), keep_or_remove="keep")
#### this is to avoid investigate the relationships between two subtypes from the same major cell types (like explore the relationships between aEndo, vEndo and capEndo)
for (cts in ctlist_to_check){
  endo_sum = stringr::str_match_all(cts,"EC") %>% unlist %>% length
  if (endo_sum > 1){
    ctlist_to_check = ctlist_to_check[!(ctlist_to_check == cts)]
  }
}
#### start cma analysis ##
inpath = cellchatDir
mat4cma_filepath = file.path(dirname(cellchatDir),"tmp." %+% so_label %+% ".CausalMediationOutput","mat4cma.Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
#### check functions to T cell exhaustion
get_functions_to_view(c("exhaust","t"), mod = "inspect", overlapKeywords = TRUE, species = "human", cat = "C7") %>% unique # , cat = "C5"
func_list = get_functions_to_view(c("GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP"), mod = "catch", overlapKeywords = FALSE, species = "human", cat = "C7")
funcToUse = func_list
#### generate matrix for causal mediation analysis ## the function is not available for multiple function test yet.
DefaultAssay(so) = "RNA"
mat4cma = generate_matrix_for_CMA(so,
                                  dataDir = inpath,
                                  feature.name = "ct.major",
                                  ctlist_to_check = ctlist_to_check,
                                  funcToUse = funcToUse,
                                  species = "human",
                                  batch = "orig.ident",
                                  ifSave = TRUE,
                                  ifReturn = TRUE,
                                  outputFile = mat4cma_filepath)
#### generate causal mediation relationships
pvalue = 0.05
cma_path = file.path(dirname(dataDir),"tmp." %+% so_label %+% ".CausalMediationOutput","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
# mat4cma = loadRDS(mat4cma_filepath)
run_CausalMediation(mat4mca, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "LC_ECNeutrophilTcell_TcellExhaustion", saveDir = cma_path) #mod = "spc",



##################################################################################
##################################################################################
################################ Main Analysis 5 #################################
##################################################################################
require(R.filesets)
require(jzhang.utils)
require(scCellMediator)
require(dplyr)
`%+%` <- paste0 
setwd("~/Documents/workspace/softwares/iRPackages/scCellMediator/")
treatedDirs = named_vec(c("LC","CRC","OvC", "BrC"), 
                        c("data/E-MTAB-6149/Lungcancer/LC_treated",
                          "data/E-MTAB-6149/Colerectalcancer/CRC_treated",
                          "data/E-MTAB-6149/Ovariancancer/OvC_treated",
                          "data/E-MTAB-6149/Breastcancer/BrC_treated"))
test_tissue = "OvC"
test_dir = treatedDirs[test_tissue]
run_combined_seurat(test_dir, 
                    save_dir = test_dir %+% "_integrated_ECMyeloidTcell/", 
                    ifReturn = FALSE, 
                    filterRules = 'CellType %in% c("EC","T_cell","Myeloid") & stringr::str_detect(.$CellTypeStatus, "CT_confirmed") & nCount_RNA > 500 & nFeature_RNA > 500', 
                    acceptRules = 'CellType,>,50', 
                    dimnum = 50, k.anchor = 5, k.filter = 50, k.weight = 20)
## accepted samples for OvC:
## analysis for OvC ## ----
## ONLY check Myeloid, T_cell and EC, 8207
test_tissue = "OvC"
test_dir = treatedDirs[test_tissue]
merged_dir = file.path(test_dir %+% "_integrated_ECMyeloidTcell/")
merged = load_so(so_dir = merged_dir)
## Purify cell types,=
markers_for_hsNeutrophil = c("ITGAM","FCGR3A","PTPRC","LY6H","FUT4") # ITGAM -> CD11B, "FCGR3A" -> CD16, PTPRC -> CD45, LY6H -> Ly6C, FUT4 -> CD15
p = DotPlot(merged, features = markers_for_hsNeutrophil, assay = "RNA", group.by = "integrated_snn_res.1") 
savePlot(pObj = p, pType = "dotplot", pName = "neutrophil_markers", pDir = merged_dir, pFmt = "pdf", pWd = 6.5, pHt = 4.5)
clusters_for_hsNeutrophil = p$data %>% filter(avg.exp > 0, pct.exp > 0, features.plot %in% c("ITGAM", "FCGR3A"), avg.exp.scaled >= 1) %>% pull(id) %>% unique %>% as.character %>% as.numeric
merged$CellType.old = merged$CellType
merged$CellType[merged$seurat_clusters %in% clusters_for_hsNeutrophil] = "Neutrophil"
## !!! Need to pay attention to the order of the cell types
cells_filtered_based_on_umap = c(Cells(merged)[merged@reductions$umap@cell.embeddings[,1] > 5 & merged$CellType %in% c("EC")],
                                 Cells(merged)[merged@reductions$umap@cell.embeddings[,1] < 5 & merged@reductions$umap@cell.embeddings[,2] > 0 & merged$CellType %in% c("Neutrophil")],
                                 Cells(merged)[merged@reductions$umap@cell.embeddings[,1] < 5 & merged@reductions$umap@cell.embeddings[,2] < 0 & merged$CellType %in% c("T_cell")])
merged = subset(merged, cells = cells_filtered_based_on_umap) # total 6337 cells
p = UMAPPlot(merged, group.by = "CellType", pt.size = 1.2, label = TRUE) +
  ggplot2::scale_color_manual(values=c(addNamedColorList,named_vec("Neutrophil","blue"))) +
  NoAxes()
savePlot(pObj = p, pType = "umap", pName = "new_celltype", pDir = merged_dir, pFmt = "pdf", pWd = 6.5, pHt = 5)
p = UMAPPlot(merged, group.by = "CellFromTumor", pt.size = 1.2, split.by = "CellFromTumor") +
  ggsci::scale_color_aaas() +
  NoAxes()
savePlot(pObj = p, pType = "umap", pName = "split_cellFromTumor", pDir = merged_dir, pFmt = "pdf", pWd = 8.4, pHt = 4)
update_so(merged, so_dir = merged_dir, meta.only = F)
## Generate UMAP plots for merged seurat object
features_to_view = c("CellFromTumor","PatientNumber","TumorSite","CellType","integrated_snn_res.1","CellTypeStatus","CellTypeM","CellTypeMarker")
addNamedColorList = loadRDS(file.path(test_dir,"ct.colorlist.rds"))
for (ifeature in features_to_view) {
  print(ifeature)
  colorlist = get_random_colorlist(merged@meta.data[,ifeature] %>% unique)
  labelon = F
  fnum = uniqlen(merged@meta.data[,ifeature])
  if (fnum > 1 & fnum < 80) {
    if (fnum < 20){ labelon = T}
    if (ifeature %in% c("CellType","CellTypeM")){
      colorlist = c(addNamedColorList, colorlist)
    }
    p = UMAPPlot(merged, group.by = ifeature, pt.size = 1.2, label = labelon) +
      ggplot2::scale_color_manual(values = colorlist) +
      NoAxes()
    # p = AugmentPlot(p, width = 6.5, height = 5, dpi = 100)
    if (ifeature %in% c("CellTypeMarker")){
      savePlot(pObj = p, pType = "umap", pName = ifeature, pDir = merged_dir, pFmt = "pdf", pWd = 12, pHt = 5)
    } else {
      savePlot(pObj = p, pType = "umap", pName = ifeature, pDir = merged_dir, pFmt = "pdf", pWd = 6.5, pHt = 5)
    }
  }
}
## Total 5635 cells accepted
## Separate cancer-originated (6607 cells) and normal-originated (5635 cells)
merged_so4cellchat = separate_so(merged, batch = "CellFromTumor")
names(merged_so4cellchat) = c("OvC","CTRL")
### For OvC
so_label = "OvC"
so = merged_so4cellchat[[so_label]]
so$ct.major <- so$CellType
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
labels <- so$ct.major
meta <- data.frame(ct = labels, row.names = names(labels), sample = so$orig.ident)
#### run cell chat analysis
cellchatDir = file.path(merged_dir,"tmp." %+% so_label %+% ".cellchat.out")
createDir(cellchatDir)
run_cellchat(so, feature.name = "ct.major", batch.name = "orig.ident",
             output = cellchatDir, 
             minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10,
             suggestedRunTimes = F, species = "human")
#### run causal mediation analysis for Endothelial cells
ctlist_to_check = generate_triple_ct_model(so$ct.major %>% unique, cells_to_fix = named_vec(c(3),c("T_cell")), keep_or_remove="keep")
#### this is to avoid investigate the relationships between two subtypes from the same major cell types (like explore the relationships between aEndo, vEndo and capEndo)
for (cts in ctlist_to_check){
  endo_sum = stringr::str_match_all(cts,"EC") %>% unlist %>% length
  if (endo_sum > 1){
    ctlist_to_check = ctlist_to_check[!(ctlist_to_check == cts)]
  }
}
#### start cma analysis ##
inpath = cellchatDir
mat4cma_filepath = file.path(dirname(cellchatDir),"tmp." %+% so_label %+% ".CausalMediationOutput","mat4cma.Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
#### check functions to T cell exhaustion
get_functions_to_view(c("exhaust","t"), mod = "inspect", overlapKeywords = TRUE, species = "human", cat = "C7") %>% unique # , cat = "C5"
func_list = get_functions_to_view(c("GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP"), mod = "catch", overlapKeywords = FALSE, species = "human", cat = "C7")
funcToUse = func_list
#### generate matrix for causal mediation analysis ## the function is not available for multiple function test yet.
DefaultAssay(so) = "RNA"
mat4cma = generate_matrix_for_CMA(so,
                                  dataDir = inpath,
                                  feature.name = "ct.major",
                                  ctlist_to_check = ctlist_to_check,
                                  funcToUse = funcToUse,
                                  species = "human",
                                  batch = "orig.ident",
                                  ifSave = TRUE,
                                  ifReturn = TRUE,
                                  outputFile = mat4cma_filepath)
#### generate causal mediation relationships
pvalue = 0.05
cma_path = file.path(dirname(cellchatDir),"tmp." %+% so_label %+% ".CausalMediationOutput","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
# mat4cma = loadRDS(mat4cma_filepath)
run_CausalMediation(mat4mca, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "OvC_ECNeutrophilTcell_TcellExhaustion", saveDir = cma_path) #mod = "spc",
### For normal
so_label = "CTRL"
so = merged_so4cellchat[[so_label]]
so$ct.major <- so$CellType
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
labels <- so$ct.major
meta <- data.frame(ct = labels, row.names = names(labels), sample = so$orig.ident)
#### run cell chat analysis
cellchatDir = file.path(merged_dir,"tmp." %+% so_label %+% ".cellchat.out")
createDir(cellchatDir)
run_cellchat(so, feature.name = "ct.major", batch.name = "orig.ident",
             output = cellchatDir, 
             minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10,
             suggestedRunTimes = F, species = "human")
#### run causal mediation analysis for Endothelial cells
ctlist_to_check = generate_triple_ct_model(so$ct.major %>% unique, cells_to_fix = named_vec(c(3),c("T_cell")), keep_or_remove="keep")
#### this is to avoid investigate the relationships between two subtypes from the same major cell types (like explore the relationships between aEndo, vEndo and capEndo)
for (cts in ctlist_to_check){
  endo_sum = stringr::str_match_all(cts,"EC") %>% unlist %>% length
  if (endo_sum > 1){
    ctlist_to_check = ctlist_to_check[!(ctlist_to_check == cts)]
  }
}
#### start cma analysis ##
inpath = cellchatDir
mat4cma_filepath = file.path(dirname(cellchatDir),"tmp." %+% so_label %+% ".CausalMediationOutput","mat4cma.Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
#### check functions to T cell exhaustion
get_functions_to_view(c("exhaust","t"), mod = "inspect", overlapKeywords = TRUE, species = "human", cat = "C7") %>% unique # , cat = "C5"
func_list = get_functions_to_view(c("GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP"), mod = "catch", overlapKeywords = FALSE, species = "human", cat = "C7")
funcToUse = func_list
#### generate matrix for causal mediation analysis ## the function is not available for multiple function test yet.
DefaultAssay(so) = "RNA"
mat4cma = generate_matrix_for_CMA(so,
                                  dataDir = inpath,
                                  feature.name = "ct.major",
                                  ctlist_to_check = ctlist_to_check,
                                  funcToUse = funcToUse,
                                  species = "human",
                                  batch = "orig.ident",
                                  ifSave = TRUE,
                                  ifReturn = TRUE,
                                  outputFile = mat4cma_filepath)
#### generate causal mediation relationships
pvalue = 0.05
cma_path = file.path(dirname(cellchatDir),"tmp." %+% so_label %+% ".CausalMediationOutput","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
# mat4cma = loadRDS(mat4cma_filepath)
run_CausalMediation(mat4mca, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "OvC_ECNeutrophilTcell_TcellExhaustion", saveDir = cma_path) #mod = "spc",



##################################################################################
##################################################################################
################################ Main Analysis 6 #################################
##################################################################################
require(R.filesets)
require(jzhang.utils)
require(scCellMediator)
require(dplyr)
`%+%` <- paste0 
setwd("~/Documents/workspace/softwares/iRPackages/scCellMediator/")
treatedDirs = named_vec(c("LC","CRC","OvC", "BrC"), 
                        c("data/E-MTAB-6149/Lungcancer/LC_treated",
                          "data/E-MTAB-6149/Colerectalcancer/CRC_treated",
                          "data/E-MTAB-6149/Ovariancancer/OvC_treated",
                          "data/E-MTAB-6149/Breastcancer/BrC_treated"))
test_tissue = "BrC"
test_dir = treatedDirs[test_tissue]
run_combined_seurat(test_dir, 
                    save_dir = test_dir %+% "_integrated_ECMyeloidTcell/", 
                    ifReturn = FALSE, 
                    filterRules = 'CellType %in% c("EC","T_cell","Myeloid") & stringr::str_detect(.$CellTypeStatus, "CT_confirmed") & nCount_RNA > 500 & nFeature_RNA > 500', 
                    acceptRules = 'CellType,>,50', 
                    dimnum = 50, k.anchor = 5, k.filter = 50, k.weight = 20)
## accepted samples for BrC:
## analysis for BrC ## ----
## ONLY check Myeloid, T_cell and EC
test_tissue = "BrC"
test_dir = treatedDirs[test_tissue]
merged_dir = file.path(test_dir %+% "_integrated_ECMyeloidTcell/")
merged = load_so(so_dir = merged_dir)
## Purify cell types,=
markers_for_hsNeutrophil = c("ITGAM","FCGR3A","PTPRC","LY6H","FUT4") # ITGAM -> CD11B, "FCGR3A" -> CD16, PTPRC -> CD45, LY6H -> Ly6C, FUT4 -> CD15
p = DotPlot(merged, features = markers_for_hsNeutrophil, assay = "RNA", group.by = "integrated_snn_res.1") 
savePlot(pObj = p, pType = "dotplot", pName = "neutrophil_markers", pDir = merged_dir, pFmt = "pdf", pWd = 6.5, pHt = 4.5)
clusters_for_hsNeutrophil = p$data %>% filter(avg.exp > 0, pct.exp > 0, features.plot %in% c("ITGAM", "FCGR3A"), avg.exp.scaled >= 1) %>% pull(id) %>% unique %>% as.character %>% as.numeric
merged$CellType.old = merged$CellType
merged$CellType[merged$seurat_clusters %in% clusters_for_hsNeutrophil] = "Neutrophil"
## !!! Need to pay attention to the order of the cell types
cells_filtered_based_on_umap = c(Cells(merged)[merged@reductions$umap@cell.embeddings[,1] > 0 & merged@reductions$umap@cell.embeddings[,2] > 0 & merged$CellType %in% c("EC")],
                                 Cells(merged)[merged@reductions$umap@cell.embeddings[,1] < 0 & merged@reductions$umap@cell.embeddings[,2] > 0 & merged$CellType %in% c("Neutrophil")],
                                 Cells(merged)[merged@reductions$umap@cell.embeddings[,2] < 0 & merged$CellType %in% c("T_cell")])
merged = subset(merged, cells = cells_filtered_based_on_umap) # total 8801 cells
p = UMAPPlot(merged, group.by = "CellType", pt.size = 1.2, label = TRUE) +
  ggplot2::scale_color_manual(values=c(addNamedColorList,named_vec("Neutrophil","blue"))) +
  NoAxes()
savePlot(pObj = p, pType = "umap", pName = "new_celltype", pDir = merged_dir, pFmt = "pdf", pWd = 6.5, pHt = 5)
p = UMAPPlot(merged, group.by = "CellFromTumor", pt.size = 1.2, split.by = "CellFromTumor") +
  ggsci::scale_color_aaas() +
  NoAxes()
savePlot(pObj = p, pType = "umap", pName = "split_cellFromTumor", pDir = merged_dir, pFmt = "pdf", pWd = 8.4, pHt = 4)
update_so(merged, so_dir = merged_dir, meta.only = F)
## Generate UMAP plots for merged seurat object
features_to_view = c("CellFromTumor","PatientNumber","TumorSite","CellType","integrated_snn_res.1","CellTypeStatus","CellTypeM","CellTypeMarker")
addNamedColorList = loadRDS(file.path(test_dir,"ct.colorlist.rds"))
for (ifeature in features_to_view) {
  print(ifeature)
  colorlist = get_random_colorlist(merged@meta.data[,ifeature] %>% unique)
  labelon = F
  fnum = uniqlen(merged@meta.data[,ifeature])
  if (fnum > 1 & fnum < 80) {
    if (fnum < 20){ labelon = T}
    if (ifeature %in% c("CellType","CellTypeM")){
      colorlist = c(addNamedColorList, colorlist)
    }
    p = UMAPPlot(merged, group.by = ifeature, pt.size = 1.2, label = labelon) +
      ggplot2::scale_color_manual(values = colorlist) +
      NoAxes()
    # p = AugmentPlot(p, width = 6.5, height = 5, dpi = 100)
    if (ifeature %in% c("CellTypeMarker")){
      savePlot(pObj = p, pType = "umap", pName = ifeature, pDir = merged_dir, pFmt = "pdf", pWd = 12, pHt = 5)
    } else {
      savePlot(pObj = p, pType = "umap", pName = ifeature, pDir = merged_dir, pFmt = "pdf", pWd = 6.5, pHt = 5)
    }
  }
}
## Total 8801 cells accepted, all cancer-originated
merged_so4cellchat = separate_so(merged, batch = "CellFromTumor")
names(merged_so4cellchat) = c("BrC")
### For BrC
so_label = "BrC"
so = merged_so4cellchat[[so_label]]
so$ct.major <- so$CellType
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
labels <- so$ct.major
meta <- data.frame(ct = labels, row.names = names(labels), sample = so$orig.ident)
#### run cell chat analysis
cellchatDir = file.path(merged_dir,"tmp." %+% so_label %+% ".cellchat.out")
createDir(cellchatDir)
run_cellchat(so, feature.name = "ct.major", batch.name = "orig.ident",
             output = cellchatDir, 
             minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10,
             suggestedRunTimes = F, species = "human")
#### run causal mediation analysis for Endothelial cells
ctlist_to_check = generate_triple_ct_model(so$ct.major %>% unique, cells_to_fix = named_vec(c(3),c("T_cell")), keep_or_remove="keep")
#### this is to avoid investigate the relationships between two subtypes from the same major cell types (like explore the relationships between aEndo, vEndo and capEndo)
for (cts in ctlist_to_check){
  endo_sum = stringr::str_match_all(cts,"EC") %>% unlist %>% length
  if (endo_sum > 1){
    ctlist_to_check = ctlist_to_check[!(ctlist_to_check == cts)]
  }
}
#### start cma analysis ##
inpath = cellchatDir
mat4cma_filepath = file.path(dirname(cellchatDir),"tmp." %+% so_label %+% ".CausalMediationOutput","mat4cma.Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
#### check functions to T cell exhaustion
get_functions_to_view(c("exhaust","t"), mod = "inspect", overlapKeywords = TRUE, species = "human", cat = "C7") %>% unique # , cat = "C5"
func_list = get_functions_to_view(c("GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP"), mod = "catch", overlapKeywords = FALSE, species = "human", cat = "C7")
funcToUse = func_list
#### generate matrix for causal mediation analysis ## the function is not available for multiple function test yet.
DefaultAssay(so) = "RNA"
mat4cma = generate_matrix_for_CMA(so,
                                  dataDir = inpath,
                                  feature.name = "ct.major",
                                  ctlist_to_check = ctlist_to_check,
                                  funcToUse = funcToUse,
                                  species = "human",
                                  batch = "orig.ident",
                                  ifSave = TRUE,
                                  ifReturn = TRUE,
                                  outputFile = mat4cma_filepath)
#### generate causal mediation relationships
pvalue = 0.05
cma_path = file.path(dirname(cellchatDir),"tmp." %+% so_label %+% ".CausalMediationOutput","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
# mat4cma = loadRDS(mat4cma_filepath)
run_CausalMediation(mat4mca, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "BrC_ECNeutrophilTcell_TcellExhaustion", saveDir = cma_path) #mod = "spc",
### No normal cells

