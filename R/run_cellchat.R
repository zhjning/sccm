# packs_to_check = c("CellChat","Seurat","R.filesets","future")
# load_packages(packs_to_check)
#
# `%+%` <- paste0

#' @title loadCellChatDB
#' @description load CellChat reference database for a specific species
#'
#' @param species current legal species are human, mouse, zebrafish.
#' @return a cell database or NULL
#' @export
loadCellChatDB <- function(species){

  packs_to_check = c("CellChat","Seurat","R.filesets","future")
  load_packages(packs_to_check)
  `%+%` <- paste0

  if (species == "human"){
    return(CellChatDB.human)
  } else if (species == "mouse"){
    return(CellChatDB.mouse)
  } else if (species == "zebrafish"){
    return(CellChatDB.zebrafish)
  } else {
    message(species %+% " is not available in CellChat reference database.")
  }
}
# CellChatDB <- loadCellChatDB("mouse")
# interaction_input <- CellChatDB$interaction



#' @title run_cellchat
#' @description run cellchat
#'
#' @param so seurat object
#' @param feature.name the cluster name used for resampling. In our examples, ct.major was used.
#' @param batch.name the batch name used for resampling. The variable used should be included in the meta.table from the analyzed seurat object. Default is "orig.ident".
#' @param output the path of the directory for cellchat results. Default is NULL.
#' @param minCNUM.tot the minCNUM.tot used for the bootstrap of seurat object. Default is 10.
#' @param sampelProp the proportion value used for the bootstrap of seurat object. Default is 0.5.
#' @param minCNUM.new the minCNUM.new used for the bootstrap of seurat object. Default is 10.
#' @param runs the runs used for the bootstrap of seurat object. Default is 10.
#' @param suggestedRunTimes the suggestedRunTimes used for the bootstrap of seurat object. Default is FALSE.
#' @param species the species used for reference CellChat database offered by CellChat.
#' @param storeResampledSO a boolean variable. If TRUE, not save the resampled seurat sub-objects used for cell chat analysis. Default is FALSE.
#' @param customedSymbolList default is NULL. Pass a named vector for single cell studies with customed defined gene names. The names of the vector are customed defined names, the corresponding values are official gene symbols.
#' @return the path of the directory for cellchat results.
#' @export
run_cellchat = function(so, feature.name, batch.name = "orig.ident", output=NULL,
                        minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10,
                        suggestedRunTimes = F, species = "mouse", 
                        storeResampledSO = FALSE, customedSymbolList = NULL){

  packs_to_check = c("CellChat","Seurat","R.filesets","future")
  load_packages(packs_to_check)
  `%+%` <- paste0

  if (!all(c(feature.name,batch.name) %in% colnames(so@meta.data))){
    message(paste_with_sep(setdiff(c(feature.name, batch.name),colnames(so@meta.data))) %+% " not found in the seurat object." )
    return(1)
  }
  if (is.null(output) | (!dir.exists(output))){
    oDir = file.path(getwd(), "cellchat.out")
    print("Output directory is not found. Saving results to the following directory:\n" %+% oDir)
  } else {
    oDir = output
  }
  createDir(oDir)
  # bug fixed 241218
  for (i in 1:uniqlen(so@meta.data[,batch.name])){
    if (as.numeric(packageVersion("Seurat")[1,1]) >= 5){
      ident_id <- unique(so@meta.data[,batch.name])[i]
      data_layers <- grep("^data", Layers(so), value = T)[i]
      so.sub <- subset(so, cells=Cells(so, assay = "RNA", layer = data_layers))
      so.sub.bts = bootstrapSO(so.sub, feature.name = feature.name, minCNUM.tot = minCNUM.tot, sampleProp = sampleProp, minCNUM.new = minCNUM.new, runs = runs, suggestedRunTimes = suggestedRunTimes,  returnNameOnly = T)
      if (storeResampledSO){
        saveRDS(so.sub.bts, file.path(oDir,ident_id %+% ".bts" %+% minCNUM.tot %+% ".prop" %+% sampleProp %+% "minCell" %+% minCNUM.new %+% ".rds"))
      }
      data.input <- GetAssayData(so, assay = "RNA", layer = data_layers)
    } else {
      so.sub <- subset(so, cells=Cells(so)[so@meta.data[,batch.name] == ident_id])
      so.sub.bts = bootstrapSO(so.sub, feature.name = feature.name, minCNUM.tot = minCNUM.tot, sampleProp = sampleProp, minCNUM.new = minCNUM.new, runs = runs, suggestedRunTimes = suggestedRunTimes,  returnNameOnly = T)
      if (storeResampledSO){
        saveRDS(so.sub.bts, file.path(oDir,ident_id %+% ".bts" %+% minCNUM.tot %+% ".prop" %+% sampleProp %+% "minCell" %+% minCNUM.new %+% ".rds"))
      }
      data.input <- GetAssayData(so, assay = "RNA", slot = "data")
    }
    if(!is.null(customedSymbolList)){
      idmapped = customedSymbolList[rownames(data.input)]
      idmapped[is.na(idmapped)] = rownames(data.input)[is.na(idmapped)]
      rownames(data.input) = idmapped        
    }
    meta <- data.frame(ct = so@meta.data[,feature.name], row.names = colnames(so), sample = so@meta.data[,batch.name])
    real_irun = 0
    for (irun in 1:length(so.sub.bts)){
      iso = so.sub.bts[[irun]]
      if ((meta[Cells(so),"ct"] %>% table %>% length) > 2){
      # if ((iso@meta.data[,feature.name] %>% table %>% length) > 2){
        real_irun = real_irun + 1
        iso.cellchat = createCellChat(object = data.input[,iso], meta = meta[iso,], group.by = "ct")
        # iso.cellchat = createCellChat(object = data.input[,Cells(iso)], meta = meta[Cells(iso),], group.by = "ct")
        iso.cellchat@DB = loadCellChatDB(species)
        iso.cellchat <- subsetData(iso.cellchat)
        future::multisession()
        iso.cellchat <- identifyOverExpressedGenes(iso.cellchat)
        iso.cellchat <- identifyOverExpressedInteractions(iso.cellchat)
        iso.cellchat <- projectData(iso.cellchat, PPI.mouse)
        iso.cellchat <- computeCommunProb(iso.cellchat, raw.use=T)
        iso.cellchat <- filterCommunication(iso.cellchat, min.cells=minCNUM.new)
        saveRDS(iso.cellchat, file.path(oDir,ident_id %+% "_run" %+% real_irun %+% ".cellchat.rds"))
        rm(iso.cellchat); gc()
        future::sequential()
      } else {
        print("run #" %+% irun %+% " skipped because cell types less than 3.")
        print(meta[iso,"ct"] %>% table)
        # print(iso@meta.data[,feature.name] %>% table)
      }
    }
  }
  saveRDS(list(oDir = oDir,
               feature.name = feature.name,
               minCNUM.tot = minCNUM.tot,
               sampleProp = sampleProp,
               minCNUM.new = minCNUM.new,
               runs = runs,
               suggestedRunTimes = suggestedRunTimes),
          file.path(oDir,"run_cellchat.config.rds"))
  return(oDir)
}
