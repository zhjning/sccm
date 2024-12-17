
#' @title run_iTALK
#' @description run gene expression to gain LR interactions
#' 
#' @param so a seurat object
#' @param feature.name the feature used to control sampling cell proportion
#' @param batch.name the feature used to separate cells from different resources
#' @param output the directory path to store results
#' @param useStoredResampledSO a boolean value. Default is FALSE. Set to TRUE to skip sample resampling by offering stored resampled seurat objects.
#' @param storedResampledSODir the path to stored resampled seurat objects. When useStoredResampledSO is TRUE, load the stored resampled seurat objects.
#' @param minCNUM.tot default is 10. The minimum number of cells in each cell cluster per feature per batch. Unsatisfied clusters will be removed.
#' @param sampleProp default is 0.5. Value ranges between 0-1. To Estimate the proportion of resampled cell number in each cell cluster.
#' @param minCNUM.new default is 10. The minimum number of cells in each resampled cell clusters per feature per batch. Unsatisfied clusters will be removed.
#' @param runs default is 10. The replicate number of resampled times.
#' @param storeResampledSO default is FALSE. a boolean value. If TRUE, store resampled seurat objects to output.
#' @param cutomedSymbolList default is NULL. Pass a named vector for single cell studies with customed defined gene names. The names of the vector are customed defined names, the corresponding values are official gene symbols.
#' @export
run_iTALK = function(so, feature.name = "ct.major", batch.name = "orig.ident",
          output = NULL, useStoredResampledSO = F, storedResampledSODir = NULL,
          minCNUM.tot = 10, sampleProp = 0.5, minCNUM.new = 10, runs = 10,
          storeResampledSO = FALSE, customedSymbolList = NULL){

  packs_to_check = c("iTALK","Seurat","R.filesets","future","Matrix")
  load_packages(packs_to_check)
  `%+%` <- paste0
  suggestedRunTimes = F

  comm_list<-c('growth factor','other','cytokine','checkpoint')

  if (is.null(output) | (!dir.exists(dirname(output)))){
    oDir = file.path(getwd(), "iTALK.out")
    print("Output directory is not found. Saving results to the following directory:\n" %+% oDir)
  } else {
    oDir = output
    print("Saving results to the following directory:\n" %+% oDir)
  }
  createDir(oDir)

  if (useStoredResampledSO){
    if (!is.null(storedResampledSODir)){
      if (dir.exists(storedResampledSODir)){
        so.sub.btsList = grep("rds", list.files(storedResampledSODir), value=T) %>%
          grep("bts", ., value = T) %>%
          grep("prop", ., value = T) %>%
          grep("minCell", ., value = T)
        for (ident_id in unique(so@meta.data[,batch.name])){
          so.sub.btsPath = grep(ident_id, so.sub.btsList, value=T)
          if (length(so.sub.btsPath) > 0){
            if (file.exists(file.path(storedResampledSODir, so.sub.btsPath))){
              so.sub.bts = loadRDS(file.path(storedResampledSODir, so.sub.btsPath))
              ## calculate
              data <-  GetAssayData(so, assay = "RNA", slot = "counts") # turn to cell x genes for rawParse
              if(!is.null(customedSymbolList)){
                idmapped = customedSymbolList[rownames(data.input)]
                idmapped[is.na(idmapped)] = rownames(data.input)[is.na(idmapped)]
                rownames(data.input) = idmapped        
              }
              meta <- data.frame(ct = so@meta.data[,feature.name], row.names = colnames(so), sample = so@meta.data[,batch.name])

              real_irun = 0
              for (irun in 1:length(so.sub.bts)){
                iso = so.sub.bts[[irun]]
                if ((meta[iso,feature.name] %>% table %>% length) > 2){
                # if ((iso@meta.data[,feature.name] %>% table %>% length) > 2){
                  real_irun = real_irun + 1
                  iso.hvg <- Matrix::rowMeans(data[,iso]) %>% sort(decreasing = T)
                  # iso.hvg <- Matrix::rowMeans(data[,Cells(iso)]) %>% sort(decreasing = T)
                  iso.hvg <- names(iso.hvg)[intersect(c(1:2500), which(iso.hvg > 0.25))]
                  iso.hvg.num <- iso.hvg %>% length
                  iso.data <- data[iso.hvg, iso] %>% t %>% data.frame
                  iso.data$cell_type <- meta[rownames(iso.data), "ct"]
                  iso.highly_exprs_genes <- iTALK::rawParse(iso.data, top_genes = iso.hvg.num, stats = "mean")
                  iso.iTALK = NULL
                  for(comm_type in comm_list){
                    iso.iTALK.tmp <- FindLR(iso.highly_exprs_genes,datatype='mean count',comm_type=comm_type)
                    iso.iTALK.tmp <- iso.iTALK.tmp[order(iso.iTALK.tmp$cell_from_mean_exprs*iso.iTALK.tmp$cell_to_mean_exprs,decreasing=T),]
                    iso.iTALK <- rbind(iso.iTALK, iso.iTALK.tmp)
                  }
                  iso.iTALK %>% distinct -> iso.iTALK
                  iso.iTALK = iso.iTALK %>%
                    dplyr::mutate(signal=cell_from_mean_exprs*cell_to_mean_exprs) %>%
                    dplyr::group_by(cell_from, cell_to) %>%
                    dplyr::summarise(score = mean(signal)) %>%
                    dplyr::arrange(-score)
                  saveRDS(iso.iTALK, file.path(oDir,ident_id %+% "_run" %+% real_irun %+% ".iTALK.rds"))
                  rm(iso.iTALK, iso, iso.hvg, iso.hvg.num, iso.data, iso.highly_exprs_genes); gc()
                  future::sequential()
                } else {
                  print("run #" %+% irun %+% " skipped because cell types less than 3.")
                  print(meta[iso,feature.name] %>% table)
                  # print(iso@meta.data[,feature.name] %>% table)
                }
              }
              saveRDS(list(oDir = oDir,
                           feature.name = feature.name,
                           minCNUM.tot = minCNUM.tot,
                           sampleProp = sampleProp,
                           minCNUM.new = minCNUM.new,
                           runs = runs,
                           suggestedRunTimes = suggestedRunTimes),
                      file.path(oDir,"run_iTALK.config.rds"))
            }
          }
        }
      }
    }
  } else {
    if (!all(c(feature.name,batch.name) %in% colnames(so@meta.data))){
      message(paste_with_sep(setdiff(c(feature.name, batch.name),colnames(so@meta.data))) %+% " not found in the seurat object." )
      return(1)
    }
    for (ident_id in unique(so@meta.data[,batch.name])){
      so.sub <- subset(so, cells=Cells(so)[so@meta.data[,batch.name] == ident_id])
      so.sub.bts = bootstrapSO(so.sub, feature.name = feature.name, minCNUM.tot = minCNUM.tot, sampleProp = sampleProp, minCNUM.new = minCNUM.new, runs = runs, suggestedRunTimes = suggestedRunTimes,  returnNameOnly = T)
      if (storeResampledSO){
        saveRDS(so.sub.bts, file.path(oDir,ident_id %+% ".bts" %+% minCNUM.tot %+% ".prop" %+% sampleProp %+% "minCell" %+% minCNUM.new %+% ".rds"))
      }
      ## calculate
      data <-  GetAssayData(so, assay = "RNA", slot = "counts") # turn to cell x genes for rawParse
      if(!is.null(customedSymbolList)){
        idmapped = customedSymbolList[rownames(data.input)]
        idmapped[is.na(idmapped)] = rownames(data.input)[is.na(idmapped)]
        rownames(data.input) = idmapped        
      }
      meta <- data.frame(ct = so@meta.data[,feature.name], row.names = colnames(so), sample = so@meta.data[,batch.name])

      real_irun = 0
      for (irun in 1:length(so.sub.bts)){
        iso = so.sub.bts[[irun]]
        if ((meta[iso,"ct"] %>% table %>% length) > 2){
        # if ((iso@meta.data[,feature.name] %>% table %>% length) > 2){
          real_irun = real_irun + 1
          iso.hvg <- Matrix::rowMeans(data[,iso]) %>% sort(decreasing = T)
          #iso.hvg <- Matrix::rowMeans(data[,Cells(iso)]) %>% sort(decreasing = T)
          iso.hvg <- names(iso.hvg)[intersect(c(1:2500), which(iso.hvg > 0.25))]
          iso.hvg.num <- iso.hvg %>% length
          iso.data <- data[iso.hvg, iso] %>% t %>% data.frame
          iso.data$cell_type <- meta[rownames(iso.data), "ct"]
          iso.highly_exprs_genes <- iTALK::rawParse(iso.data, top_genes = iso.hvg.num, stats = "mean")
          iso.iTALK = NULL
          for(comm_type in comm_list){
            iso.iTALK.tmp <- FindLR(iso.highly_exprs_genes,datatype='mean count',comm_type=comm_type)
            iso.iTALK.tmp <- iso.iTALK.tmp[order(iso.iTALK.tmp$cell_from_mean_exprs*iso.iTALK.tmp$cell_to_mean_exprs,decreasing=T),]
            iso.iTALK <- rbind(iso.iTALK, iso.iTALK.tmp)
          }
          iso.iTALK %>% distinct -> iso.iTALK
          iso.iTALK = iso.iTALK %>%
            dplyr::mutate(signal=cell_from_mean_exprs*cell_to_mean_exprs) %>%
            dplyr::group_by(cell_from, cell_to) %>%
            dplyr::summarise(score = mean(signal)) %>% # use mean instead of sum
            dplyr::arrange(-score)
          saveRDS(list(scoremat = iso.iTALK, cellnames = iso), file.path(oDir,ident_id %+% "_run" %+% real_irun %+% ".iTALK.rds"))
          rm(iso.iTALK, iso, iso.hvg, iso.hvg.num, iso.data, iso.highly_exprs_genes); gc()
          future::sequential()
        } else {
          print("run #" %+% irun %+% " skipped because cell types less than 3.")
          print(meta[iso,feature.name] %>% table)
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
            file.path(oDir,"run_iTALK.config.rds"))
  }
}


