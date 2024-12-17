
#########################################
## A simple case to run scCellMediator ##
#########################################

require(Seurat)
require(dplyr)
require(tibble)
require(stringr)
require(R.filesets)
require(scCellMediator)
require(CellChat)
require(iTALK)

workDir = "~/Desktop/scCellMediator/demo/"
dataDir = "~/Desktop/scCellMediator/demo/data"
scriptDir = "~/Desktop/scCellMediator/demo/scripts"

sapply(c(workDir, dataDir, scriptDir), createDir)
setwd(workDir)

## import raw 10x dataset and meta
dsname = "PBMC" # dataset's name
meta.raw = read.delim(file.path(dataDir, dsname, "meta.txt"))[-1,]
mtx = readMM(file.path(dataDir, dsname, "counts.umi.txt.gz"))
cells = read.delim(file.path(dataDir, dsname, "cells.umi.txt"), header=F)
features = read.delim(file.path(dataDir, dsname, "genes.umi.txt"), header=F)
colnames(mtx) = cells[,1]
rownames(mtx) = features[,1]
so = CreateSeuratObject(mtx, assay = "RNA")
meta.raw = meta.raw %>% 
  filter(NAME %in% Cells(so)) %>% 
  column_to_rownames("NAME")
## pretreating dataset and calling seurat to generate the seurat object for following analysis
so = subset(so, cells = rownames(meta.raw))
so = AddMetaData(so, meta.raw[Cells(so),], colnames(meta.raw)) # so@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt")] = NULL
solist = separate_so(so = so, batch = "orig.ident")
ct.colorlist = get_random_colorlist(purrr::reduce(lapply(solist, function(tmp.so) tmp.so$CellType %>% unique), union))
saveRDS(ct.colorlist, file.path(treatedDir, dsname, "ct.colorlist.rds"))
## separating single cells according to the original dataset sources and generating seurat objects for single cells per source
treatedDir = file.path(workDir, "output")
createDir(treatedDir)
for (i in 1:length(solist)){
  sname = names(solist)[i]
  tmp.so = solist[[i]]
  run_seurat(tmp.so, label = sname, features_to_view = c("orig.ident", "CellType"),
             save = T, so_dir = file.path(treatedDir, dsname), treat.cellname = F,
             mtlabel = "-MT-", min_avg = 0.01, assay = "RNA", addNamedColorList = ct.colorlist,
             scale_factor = 1e6, min_bin = 50, min_cutoff = c(0.01, Inf), calcClustMarker = F,
             dimnum = "auto", dim_sd = 1.25, dim_pval = 0.01, clustRes = 0.8)
}

## Study the relationship between CD4+ T cells, CD16+ Monocytes and CD14+ Monocytes ----
## generating a combined seurat object for triple-cell causal mediation analysis only including the three cell types to study
run_combined_seurat(solist_dir = file.path(treatedDir, dsname), ## this folder only allows to store separated datasets rather than the integrated ones
                    save_dir = paste0(file.path(treatedDir, dsname), "_integrated_CD14Mono_CD16Mono_CD4T/"), 
                    ifReturn = FALSE, 
                    filterRules = 'CellType %in% c("CD14+ monocyte","CD4+ T cell","CD16+ monocyte") & nGene > 500 & nUMI > 1000', 
                    acceptRules = 'CellType,>,20', 
                    dimnum = 50, k.anchor = 5, k.filter = 50, k.weight = 20)
## pretreating dataset and calling seurat to generate the seurat object for following analysis
so = load_so(file.path(treatedDir, paste0(dsname, "_integrated_CD14Mono_CD16Mono_CD4T/")))
p = UMAPPlot(so, group.by = "CellType", cols = ct.colorlist)
savePlot(p, file.path(file.path(treatedDir, paste0(dsname, "_integrated_CD14Mono_CD16Mono_CD4T/"))), "umap", "celltype", pWd = 6.3, pHt = 3, "pdf")
## running cellchat to estimate the interactions
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
meta <- data.frame(ct = so$CellType, row.names = Cells(so), sample = so$orig.ident)
ctlist = unique(so$CellType)
## generate potential regulatory relationships
ctlist_to_check <- generate_triple_ct_model(ctlist,
                                            cells_to_fix = named_vec(c(3), c("CD16+ monocyte")),
                                            keep_or_remove = "keep")
## search functions about the cell activation and differentiation of monocytes 
get_functions_to_view("monocyte", mod = "inspect", overlapKeywords = TRUE,  species = "human", cat = "C5") %>% unique
func_list = sapply(c("GOBP_MONOCYTE_DIFFERENTIATION","GOBP_MONOCYTE_ACTIVATION"),
                   function(x) get_functions_to_view(x, mod = "catch", overlapKeywords = FALSE, species = "human", cat = "C5"))
funcToUse = func_list %>% unlist %>% unique 
## running cell chat analysis
cellchatDir = file.path(treatedDir, paste0(dsname, "_integrated_CD14Mono_CD16Mono_CD4T/"), "cellchat.out")
iTALKDir = file.path(treatedDir, paste0(dsname, "_integrated_CD14Mono_CD16Mono_CD4T/"), "iTALK.out")
createDir(cellchatDir)
createDir(iTALKDir)
customedSymbolList = named_vec(rownames(data.input), 
                               rownames(data.input) %>% str_replace("ENSG[0-9]{11}-","")) #change EnsemblID+Symbol to Symbol
run_cellchat(so, feature.name = "CellType", batch.name = "orig.ident", output = cellchatDir, 
             minCNUM.tot = 5, sampleProp = 0.6, minCNUM.new = 5, runs = 10, species = "human", 
             storeResampledSO = F, suggestedRunTimes = F, customedSymbolList = NULL)
run_iTALK(so, feature.name = "CellType", batch.name = "orig.ident", output = iTALKDir, 
          useStoredResampledSO = F, storedResampledSODir = NULL, 
          minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10, 
          storeResampledSO = FALSE)

## running cma analysis
funcToUse = func_list %>% unlist %>% unique 
mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput", "mat4cma.CD16Mono_receiver.rds")
createDir(dirname(mat4cma_filepath))
mat4cma <- generate_matrix_for_CMA(so, 
                                   dataDir = cellchatDir, 
                                   feature.name = "CellType", 
                                   int.type = "cellchat.raw", 
                                   ctlist_to_check = ctlist_to_check,
                                   funcToUse = list(CD16Mono=funcToUse),
                                   species = "human",
                                   batch = "orig.ident",
                                   customedSymbolList = customedSymbolList,
                                   ifSave = T,
                                   ifReturn = T,
                                   outputFile = mat4cma_filepath)

funcToUse = func_list %>% unlist %>% unique 
mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput2", "mat4cma.CD16Mono_receiver.rds")
createDir(dirname(mat4cma_filepath))
mat4cma <- generate_matrix_for_CMA(so, 
                                   dataDir = cellchatDir, 
                                   feature.name = "CellType", 
                                   int.type = "cellchat.p", 
                                   ctlist_to_check = ctlist_to_check,
                                   funcToUse = list(CD16Mono=funcToUse),
                                   species = "human",
                                   batch = "orig.ident",
                                   customedSymbolList = customedSymbolList,
                                   ifSave = T,
                                   ifReturn = T,
                                   outputFile = mat4cma_filepath)

funcToUse = func_list %>% unlist %>% unique 
mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput3", "mat4cma.CD16Mono_receiver.rds")
createDir(dirname(mat4cma_filepath))
mat4cma <- generate_matrix_for_CMA(so, 
                                   dataDir = iTALKDir, 
                                   feature.name = "CellType", 
                                   int.type = "iTALK", 
                                   ctlist_to_check = ctlist_to_check,
                                   funcToUse = list(CD16Mono=funcToUse),
                                   species = "human",
                                   batch = "orig.ident",
                                   customedSymbolList = customedSymbolList,
                                   ifSave = T,
                                   ifReturn = T,
                                   outputFile = mat4cma_filepath)

#### generate causal mediation relationships
mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput", "mat4cma.CD16Mono_receiver.rds")
mat4cma = loadRDS(mat4cma_filepath)
pvalue = 0.05#1e-5#1e-3
cma_path = file.path(dirname(cellchatDir), "CMOutput","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
run_CausalMediation(mat4cma, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "CD14Mono_CD4Tcell_CD16Monofunctions", saveDir = cma_path) #mod = "spc",

mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput2", "mat4cma.CD16Mono_receiver.rds")
mat4cma = loadRDS(mat4cma_filepath)
pvalue = 0.05#1e-5#1e-3
cma_path = file.path(dirname(cellchatDir), "CMOutput2","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
run_CausalMediation(mat4cma, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "CD14Mono_CD4Tcell_CD16Monofunctions", saveDir = cma_path) #mod = "spc",

mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput3", "mat4cma.CD16Mono_receiver.rds")
mat4cma = loadRDS(mat4cma_filepath)
pvalue = 0.05#1e-5#1e-3
cma_path = file.path(dirname(cellchatDir), "CMOutput3","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
run_CausalMediation(mat4cma, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "CD14Mono_CD4Tcell_CD16Monofunctions", saveDir = cma_path) #mod = "spc",

## check iso.cellchat@data.signaling
# update_so(so,file.path(treatedDir, dsname), meta.only = F)
## Study the relationship between CD4+ T cells, CD14+ Monocytes and Dentritic cells ----
## generating a combined seurat object for triple-cell causal mediation analysis only including the three cell types to study
run_combined_seurat(solist_dir = file.path(treatedDir, dsname), ## this folder only allows to store separated datasets rather than the integrated ones
                    save_dir = paste0(file.path(treatedDir, dsname), "_integrated_CD14Mono_DC_CD4T/"), 
                    ifReturn = FALSE, 
                    filterRules = 'CellType %in% c("CD14+ monocyte","Dendritic cell","CD4+ T cell") & nGene > 500 & nUMI > 1000', 
                    acceptRules = 'CellType,>,20', 
                    dimnum = 50, k.anchor = 5, k.filter = 50, k.weight = 20)
## pretreating dataset and calling seurat to generate the seurat object for following analysis
so = load_so(file.path(treatedDir, paste0(dsname, "_integrated_CD14Mono_DC_CD4T/")))
p = UMAPPlot(so, group.by = "CellType", cols = ct.colorlist)
savePlot(p, file.path(file.path(treatedDir, paste0(dsname, "_integrated_CD14Mono_DC_CD4T/"))), "umap", "celltype", pWd = 6.3, pHt = 3, "pdf")
## running cellchat to estimate the interactions
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
meta <- data.frame(ct = so$CellType, row.names = Cells(so), sample = so$orig.ident)
ctlist = unique(so$CellType)
## generate potential regulatory relationships
ctlist_to_check <- generate_triple_ct_model(ctlist,
                                            cells_to_fix = named_vec(c(3), c("CD4+ T cell")),
                                            keep_or_remove = "keep")
## search functions about the cell activation and differentiation of monocytes 
get_functions_to_view("CD4", mod = "inspect", overlapKeywords = TRUE,  species = "human", cat = "C5") %>% unique
func_list = sapply(c("GOBP_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION",
                     "GOBP_CD4_POSITIVE_ALPHA_BETA_T_CELL_CYTOKINE_PRODUCTION"),
                   function(x) get_functions_to_view(x, mod = "catch", overlapKeywords = FALSE, species = "human", cat = "C5"))
funcToUse = func_list %>% unlist %>% unique 
## running cell chat analysis
cellchatDir = file.path(treatedDir, paste0(dsname, "_integrated_DC_CD14Mono_CD4T/"), "cellchat.out")
iTALKDir = file.path(treatedDir, paste0(dsname, "_integrated_DC_CD14Mono_CD4T/"), "iTALK.out")
createDir(cellchatDir)
createDir(iTALKDir)
customedSymbolList = named_vec(rownames(data.input), 
                               rownames(data.input) %>% str_replace("ENSG[0-9]{11}-","")) #change EnsemblID+Symbol to Symbol
run_cellchat(so, feature.name = "CellType", batch.name = "orig.ident", output = cellchatDir, 
             minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10, species = "human", 
             storeResampledSO = F, suggestedRunTimes = F, customedSymbolList = customedSymbolList)
run_iTALK(so, feature.name = "CellType", batch.name = "orig.ident", output = iTALKDir, 
          useStoredResampledSO = F, storedResampledSODir = NULL, customedSymbolList = customedSymbolList,
          minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10, 
          storeResampledSO = FALSE)

## running cma analysis
funcToUse = func_list %>% unlist %>% unique
mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput", "mat4cma.CD4Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
mat4cma <- generate_matrix_for_CMA(so, 
                                   dataDir = cellchatDir, 
                                   feature.name = "CellType", 
                                   int.type = "cellchat.raw", 
                                   ctlist_to_check = ctlist_to_check,
                                   funcToUse = list(CD4Tcell=funcToUse),
                                   species = "human",
                                   batch = "orig.ident",
                                   customedSymbolList = customedSymbolList,
                                   ifSave = T,
                                   ifReturn = T,
                                   outputFile = mat4cma_filepath)

funcToUse = func_list %>% unlist %>% unique 
mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput2", "mat4cma.CD4Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
mat4cma <- generate_matrix_for_CMA(so, 
                                   dataDir = cellchatDir, 
                                   feature.name = "CellType", 
                                   int.type = "cellchat.p", 
                                   ctlist_to_check = ctlist_to_check,
                                   funcToUse = list(CD4Tcell=funcToUse),
                                   species = "human",
                                   batch = "orig.ident",
                                   customedSymbolList = customedSymbolList,
                                   ifSave = T,
                                   ifReturn = T,
                                   outputFile = mat4cma_filepath)

funcToUse = func_list %>% unlist %>% unique 
mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput3", "mat4cma.CD4Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
mat4cma <- generate_matrix_for_CMA(so, 
                                   dataDir = iTALKDir, 
                                   feature.name = "CellType", 
                                   int.type = "iTALK", 
                                   ctlist_to_check = ctlist_to_check,
                                   funcToUse = list(CD4Tcell=funcToUse),
                                   species = "human",
                                   batch = "orig.ident",
                                   customedSymbolList = customedSymbolList,
                                   ifSave = T,
                                   ifReturn = T,
                                   outputFile = mat4cma_filepath)

#### generate causal mediation relationships
mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput", "mat4cma.CD4Tcell_receiver.rds")
mat4cma = loadRDS(mat4cma_filepath)
pvalue = 0.05#1e-5#1e-3
cma_path = file.path(dirname(cellchatDir), "CMOutput","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
run_CausalMediation(mat4cma, pvalue=pvalue, useColor = TRUE, colList = NULL,
                    ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "DC_CD14Mono_CD4Tcellfunctions", saveDir = cma_path) #mod = "spc",

mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput2", "mat4cma.CD4Tcell_receiver.rds")
mat4cma = loadRDS(mat4cma_filepath)
pvalue = 0.05#1e-5#1e-3
cma_path = file.path(dirname(cellchatDir), "CMOutput2","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
run_CausalMediation(mat4cma, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "DC_CD14Mono_CD4Tcellfunctions", saveDir = cma_path) #mod = "spc",

mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput3", "mat4cma.CD4Tcell_receiver.rds")
mat4cma = loadRDS(mat4cma_filepath)
pvalue = 0.05#1e-5#1e-3
cma_path = file.path(dirname(cellchatDir), "CMOutput3","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
run_CausalMediation(mat4cma, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "DC_CD14Mono_CD4Tcellfunctions", saveDir = cma_path) #mod = "spc",




## Study the relationship between CD4+ T cells, CD16+ Monocytes and Dentritic cells ----
## generating a combined seurat object for triple-cell causal mediation analysis only including the three cell types to study
run_combined_seurat(solist_dir = file.path(treatedDir, dsname), ## this folder only allows to store separated datasets rather than the integrated ones
                    save_dir = paste0(file.path(treatedDir, dsname), "_integrated_CD16Mono_DC_CD4T/"), 
                    ifReturn = FALSE, 
                    filterRules = 'CellType %in% c("CD16+ monocyte","Dendritic cell","CD4+ T cell") & nGene > 500 & nUMI > 1000', 
                    acceptRules = 'CellType,>,20', 
                    dimnum = 50, k.anchor = 5, k.filter = 50, k.weight = 20)
## pretreating dataset and calling seurat to generate the seurat object for following analysis
so = load_so(file.path(treatedDir, paste0(dsname, "_integrated_CD16Mono_DC_CD4T/")))
p = UMAPPlot(so, group.by = "CellType", cols = ct.colorlist)
savePlot(p, file.path(file.path(treatedDir, paste0(dsname, "_integrated_CD16Mono_DC_CD4T/"))), "umap", "celltype", pWd = 6.3, pHt = 3, "pdf")
## running cellchat to estimate the interactions
data.input <- GetAssayData(so, assay = "RNA", slot = "data")
meta <- data.frame(ct = so$CellType, row.names = Cells(so), sample = so$orig.ident)
ctlist = unique(so$CellType)
## generate potential regulatory relationships
ctlist_to_check <- generate_triple_ct_model(ctlist,
                                            cells_to_fix = named_vec(c(3), c("CD4+ T cell")),
                                            keep_or_remove = "keep")
## search functions about the cell activation and differentiation of monocytes 
get_functions_to_view("CD4", mod = "inspect", overlapKeywords = TRUE,  species = "human", cat = "C5") %>% unique
func_list = sapply(c("GOBP_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION",
                     "GOBP_CD4_POSITIVE_ALPHA_BETA_T_CELL_CYTOKINE_PRODUCTION"),
                   function(x) get_functions_to_view(x, mod = "catch", overlapKeywords = FALSE, species = "human", cat = "C5"))
funcToUse = func_list %>% unlist %>% unique 
## running cell chat analysis
cellchatDir = file.path(treatedDir, paste0(dsname, "_integrated_DC_CD16Mono_CD4T/"), "cellchat.out")
iTALKDir = file.path(treatedDir, paste0(dsname, "_integrated_DC_CD16Mono_CD4T/"), "iTALK.out")
createDir(cellchatDir)
createDir(iTALKDir)
customedSymbolList = named_vec(rownames(data.input), 
                               rownames(data.input) %>% str_replace("ENSG[0-9]{11}-","")) #change EnsemblID+Symbol to Symbol
run_cellchat(so, feature.name = "CellType", batch.name = "orig.ident", output = cellchatDir, 
             minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10, species = "human", 
             storeResampledSO = F, suggestedRunTimes = F, customedSymbolList = customedSymbolList)
run_iTALK(so, feature.name = "CellType", batch.name = "orig.ident", output = iTALKDir, 
          useStoredResampledSO = F, storedResampledSODir = NULL, customedSymbolList = customedSymbolList,
          minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10, 
          storeResampledSO = FALSE)

## running cma analysis
funcToUse = func_list %>% unlist %>% unique 
mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput", "mat4cma.CD4Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
mat4cma <- generate_matrix_for_CMA(so, 
                                   dataDir = cellchatDir, 
                                   feature.name = "CellType", 
                                   int.type = "cellchat.raw", 
                                   ctlist_to_check = ctlist_to_check,
                                   funcToUse = list(CD4Tcell=funcToUse),
                                   species = "human",
                                   batch = "orig.ident",
                                   customedSymbolList = customedSymbolList,
                                   ifSave = T,
                                   ifReturn = T,
                                   outputFile = mat4cma_filepath)

funcToUse = func_list %>% unlist %>% unique 
mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput2", "mat4cma.CD4Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
mat4cma <- generate_matrix_for_CMA(so, 
                                   dataDir = cellchatDir, 
                                   feature.name = "CellType", 
                                   int.type = "cellchat.p", 
                                   ctlist_to_check = ctlist_to_check,
                                   funcToUse = list(CD4Tcell=funcToUse),
                                   species = "human",
                                   batch = "orig.ident",
                                   customedSymbolList = customedSymbolList,
                                   ifSave = T,
                                   ifReturn = T,
                                   outputFile = mat4cma_filepath)

funcToUse = func_list %>% unlist %>% unique 
mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput3", "mat4cma.CD4Tcell_receiver.rds")
createDir(dirname(mat4cma_filepath))
mat4cma <- generate_matrix_for_CMA(so, 
                                   dataDir = iTALKDir, 
                                   feature.name = "CellType", 
                                   int.type = "iTALK", 
                                   ctlist_to_check = ctlist_to_check,
                                   funcToUse = list(CD4Tcell=funcToUse),
                                   species = "human",
                                   batch = "orig.ident",
                                   customedSymbolList = customedSymbolList,
                                   ifSave = T,
                                   ifReturn = T,
                                   outputFile = mat4cma_filepath)

#### generate causal mediation relationships
mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput", "mat4cma.CD4Tcell_receiver.rds")
mat4cma = loadRDS(mat4cma_filepath)
pvalue = 0.05#1e-5#1e-3
cma_path = file.path(dirname(cellchatDir), "CMOutput","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
run_CausalMediation(mat4cma, pvalue=pvalue, useColor = TRUE, colList = NULL,
                    ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "DC_CD16Mono_CD4Tcellfunctions", saveDir = cma_path) #mod = "spc",

mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput2", "mat4cma.CD4Tcell_receiver.rds")
mat4cma = loadRDS(mat4cma_filepath)
pvalue = 0.05#1e-5#1e-3
cma_path = file.path(dirname(cellchatDir), "CMOutput2","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
run_CausalMediation(mat4cma, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "DC_CD16Mono_CD4Tcellfunctions", saveDir = cma_path) #mod = "spc",

mat4cma_filepath = file.path(dirname(cellchatDir), "CMOutput3", "mat4cma.CD4Tcell_receiver.rds")
mat4cma = loadRDS(mat4cma_filepath)
pvalue = 0.05#1e-5#1e-3
cma_path = file.path(dirname(cellchatDir), "CMOutput3","output_" %+% pvalue %+% "_bothdict")
createDir(cma_path)
run_CausalMediation(mat4cma, pvalue=pvalue, ifPlot = TRUE, ifSavePlot = TRUE, ifCombinePlot = TRUE,
                    ifReturnPlots = FALSE, ifReturnFitModels = FALSE, ifPlotOnlySameDirect = FALSE,
                    ifSaveFitModels = TRUE, datalabel = "DC_CD16Mono_CD4Tcellfunctions", saveDir = cma_path) #mod = "spc",


