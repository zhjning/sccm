
#' @title generate_triple_ct_model
#' @description generate triple cell type models
#' 
#' cells to fix: named list, name for cell position, value for cell name, use Null or "" as the name if don't need to filter the fixed position;
#'
#' @param alist_of_ct a character list of the name of cells
#' @param cells_to_fix a named list or character to filter the the result. position 1,2,3 represent sender, mediator, receiver, respectively. Names of the vector corresponding to the name of the cell types used in alist_of_ct.
#' @param keep_or_remove two values can be identified. "keep" to keep the combinations satisfying cells_to_fix's rules. "remove" to remove the combinations satifying cells_to_fix's rules.
#' @return a combination of the sender, mediator, and receiver to use for the causal mediation analysis.
#' @export
generate_triple_ct_model = function(alist_of_ct, cells_to_fix, keep_or_remove="keep"){

  packs_to_check = c("combinat","lavaan","stringr","cowplot")
  load_packages(packs_to_check)

  if (length(alist_of_ct) < 3){print("The input must contain at lease three types of cells."); return(1);}
  combinations_with_no_order = combn(alist_of_ct, 3)
  if (length(alist_of_ct) == 3){combinations_with_no_order = data.frame(combinations_with_no_order)}
  combinations_with_order = cbind(combinations_with_no_order,
                                  combinations_with_no_order[c(1,3,2),],
                                  combinations_with_no_order[c(3,1,2),],
                                  combinations_with_no_order[c(3,2,1),],
                                  combinations_with_no_order[c(2,3,1),],
                                  combinations_with_no_order[c(2,1,3),])
  global_filters = which(is.null(names(cells_to_fix)) | names(cells_to_fix)=="")
  sender_filters = which(names(cells_to_fix) == 1)
  mediator_filters = which(names(cells_to_fix) == 2)
  receiver_filters = which(names(cells_to_fix) == 3)

  global_flag = 0; sender_flag = 0; mediator_flag = 0; receiver_flag = 0;

  if (keep_or_remove == "keep"){
    if (length(global_filters) > 0){
      global_combinations = combinations_with_order[,combinations_with_order[1,] %in% cells_to_fix[global_filters] | combinations_with_order[2,] %in% cells_to_fix[global_filters] | combinations_with_order[3,] %in% cells_to_fix[global_filters]]
      global_flag = 1
    }
    if (length(sender_filters) > 0){
      sender_combinations = combinations_with_order[,combinations_with_order[1,] %in% cells_to_fix[sender_filters]]
      sender_flag = 1
    }
    if (length(mediator_filters) > 0){
      mediator_combinations = combinations_with_order[,combinations_with_order[2,] %in% cells_to_fix[mediator_filters]]
      mediator_flag = 1
    }
    if (length(receiver_filters) > 0){
      receiver_combinations = combinations_with_order[,combinations_with_order[3,] %in% cells_to_fix[receiver_filters]]
      receiver_flag = 1
    }

    new_combinations = data.frame()
    if (global_flag == 1){
      if (ncol(new_combinations) == 0){
        new_combinations = global_combinations
      } else {
        new_combinations = cbind(new_combinations, global_combinations)
      }
    }
    if (sender_flag == 1){
      if (ncol(new_combinations) == 0){
        new_combinations = sender_combinations
      } else {
        new_combinations = cbind(new_combinations, sender_combinations)
      }
    }
    if (mediator_flag == 1){
      if (ncol(new_combinations) == 0){
        new_combinations = mediator_combinations
      } else {
        new_combinations = cbind(new_combinations, mediator_combinations)
      }
    }
    if (receiver_flag == 1){
      if (ncol(new_combinations) == 0){
        new_combinations = receiver_combinations
      } else {
        new_combinations = cbind(new_combinations, receiver_combinations)
      }
    }
    new_combinations = apply(new_combinations,2,function(x) paste0(x, collapse=",")) %>% unlist %>% unique
    return(new_combinations)
  } else if (keep_or_remove=="remove"){
    new_combinations = combinations_with_order
    if (length(global_filters) > 0){
      new_combinations = new_combinations[,! (new_combinations[1,] %in% cells_to_fix[global_filters] | new_combinations[2,] %in% cells_to_fix[global_filters] | new_combinations[3,] %in% cells_to_fix[global_filters])]
    }
    if (length(sender_filters) > 0){
      new_combinations = new_combinations[,! (new_combinations[1,] %in% cells_to_fix[sender_filters])]
    }
    if (length(mediator_filters) > 0){
      new_combinations = new_combinations[,! (new_combinations[2,] %in% cells_to_fix[mediator_filters])]
    }
    if (length(receiver_filters) > 0){
      new_combinations = new_combinations[,! (new_combinations[3,] %in% cells_to_fix[receiver_filters])]
    }
    new_combinations = apply(new_combinations,2,function(x) paste0(x, collapse=",")) %>% unlist %>% unique
    return(new_combinations)
  }
}



#'@title normalized_interaction_strength
#'@param int.scores the vector for interaction scores.
#'@return the normalized interaction score. New range: 0-1.
normalized_interaction_strength = function(int.scores){
  int.scores = as.numeric(int.scores)
  int.scores.new = sapply(int.scores, function(x) (x-min(int.scores))/diff(range(int.scores))) %>% range
  return(int.scores.new)
}



#' @title generate_matrix_for_CMA
#' @description Generate the matrix used for Causal mediation analysis
#'
#' @param so a seurat object
#' @param feature.name the variable name of cell types used in the seurat object.
#' @param int.type the method used for measure outbound interaction strength. 'cellchat.raw' for interaction number estimated by cellchatDB. 'cellchat.p' for cellchat estimated interaction possibility scores. 'iTALK' for iTALK scores.
#' @param dataDir the directory of cell chat results
#' @param ctlist_to_check a character of cell combinations, each elements contains 3 cell types connected by comma.
#' @param funcToUse a list of grouped genes with list names show function names. Obtained from get_functions_to_view.
#' @param species the species for the reference database of cellchat database.
#' @param batch the name of the batch feature for the seurat object. Default is "orig.ident".
#' @param ifSave a boolean variable. save result matrix to outputFile. Default is FALSE.
#' @param ifReturn a boolean variable. Only works when ifSave is TRUE. if ifReturn is TRUE, return the maxtrix when save the cma matrix to outputFile. Default is TRUE.
#' @param outputFile default is NULL. If ifSave is TRUE and the directory assigned for outputFile doesn't exist, use the current working directory and save the results to mat4cma.rds.
#' @param customedSymbolList default is NULL. Pass a named vector for single cell studies with customed defined gene names. The names of the vector are customed defined names, the corresponding values are official gene symbols.
#' @return a data.frame which contains the interaction strength calculated by CellChat
#' @export
generate_matrix_for_CMA = function(so, feature.name, int.type, dataDir, ctlist_to_check, funcToUse, species = "mouse", batch = "orig.ident", ifSave = FALSE, ifReturn = TRUE, outputFile = NULL, customedSymbolList = NULL){

  packs_to_check = c("combinat","lavaan","stringr","cowplot")
  load_packages(packs_to_check)

  funcNames <- names(funcToUse)
  mat4cma = data.frame() # initiate the matrix for causal mediation analysis
  if (grep("cellchat",int.type,value=T) %>% length > 0){

    CellChatDB <- loadCellChatDB(species=species)
    interaction_input <- CellChatDB$interaction

    for (ident_id in unique(so@meta.data[,batch])){
      if ((grep(ident_id %+% "_run", list.files(dataDir), value = T) %>% length) > 0){
        maxrun = grep(ident_id %+% "_run", list.files(dataDir), value = T) %>%
          stringr::str_replace(ident_id %+% "_run", "") %>%
          stringr::str_replace(".cellchat.rds", "") %>%
          as.numeric %>% max
        for (irun in 1:maxrun){
          iso.cellchat <- loadRDS(file.path(dataDir,ident_id %+% "_run" %+% irun %+% ".cellchat.rds"))

          if (int.type=="cellchat.raw"){
            tot.ligands <- interaction_input$ligand %>% unique
            tot.receptors <- interaction_input$receptor %>% unique # not used

            iso <- subset(so, cells=colnames(iso.cellchat@data))
            Seurat::DefaultAssay(iso) <- "RNA"
            if (! is.null(customedSymbolList)){
              iso <- Seurat::AddModuleScore(iso,
                                            name = c("tot.ligands", "tot.receptors"),
                                            features = list(tot.ligands = intersect(names(customedSymbolList)[customedSymbolList %in% tot.ligands], rownames(iso)), 
                                                            tot.receptors = intersect(names(customedSymbolList)[customedSymbolList %in% tot.receptors], rownames(iso))),
                                            seed = 42)
            } else {
              iso <- Seurat::AddModuleScore(iso, 
                                            name = c("tot.ligands", "tot.receptors"),
                                            features = list(tot.ligands = intersect(tot.ligands, rownames(iso)),
                                                            tot.receptors = intersect(tot.receptors, rownames(iso))),
                                            seed = 42)
            }
            colnames(iso@meta.data)[c((ncol(iso@meta.data)-1):ncol(iso@meta.data))] <- colnames(iso@meta.data)[c((ncol(iso@meta.data)-1):ncol(iso@meta.data))] %>%
              stringr::str_replace("[0-9]","")
            Seurat::DefaultAssay(iso) <- "RNA"
            if (! is.null(customedSymbolList)){
              ifuncToUse <- lapply(funcToUse, function(x) intersect(rownames(iso), names(customedSymbolList)[customedSymbolList %in% x]))
            } else {
              ifuncToUse <- lapply(funcToUse, function(x) intersect(rownames(iso), x))
            }
            iso <- Seurat::AddModuleScore(iso, features = ifuncToUse, name = funcNames, seed = 42)
            colnames(iso@meta.data)[ncol(iso@meta.data)] <- colnames(iso@meta.data)[ncol(iso@meta.data)] %>%
              stringr::str_sub(.,start = 1,end = -2)

            for (cts in ctlist_to_check){
              message("Running for data:" %+% ident_id %+% " ," %+% "run #" %+% irun %+% ", cells:" %+% cts)
              cts <- cts %>% strsplit(",") %>% unlist
              Sender <- cts[1]; Mediator <- cts[2]; Receiver <- cts[3];
              # typo fixed here 241218
              if (length(intersect(rownames(iso.cellchat@net$prob), c(Sender,Mediator,Receiver))) == 3){
                sender.signal <- median(iso@meta.data[iso@meta.data[,feature.name] == Sender,"tot.ligands"]) # X
                mediator.signal <- median(iso@meta.data[iso@meta.data[,feature.name] == Mediator,"tot.ligands"]) # M
                receiverFuncOfInterest.signal <- median(iso@meta.data[iso@meta.data[,feature.name] == Receiver,funcNames]) # Y
                mat4cma <- rbind(mat4cma,
                                 c(ident_id, irun, paste0(cts, collapse=","),
                                   sender.signal,
                                   mediator.signal,
                                   receiverFuncOfInterest.signal))
                colnames(mat4cma) = c("Sample_id","Rep_id","CA_rls",
                                      "X","M","Y") # CA_rls stands for causual mediation relationship
              }
            }
          } else if (int.type=="cellchat.p"){
            iso.cellchat <- aggregateNet(iso.cellchat,return.object = T)

            iso <- subset(so, cells=colnames(iso.cellchat@data))
            Seurat::DefaultAssay(iso) <- "RNA"
            if (! is.null(customedSymbolList)){
              ifuncToUse <- lapply(funcToUse, function(x) intersect(rownames(iso), names(customedSymbolList)[customedSymbolList %in% x]))
            } else {
              ifuncToUse <- lapply(funcToUse, function(x) intersect(rownames(iso), x))
            }
            iso <- Seurat::AddModuleScore(iso, features = ifuncToUse, name = funcNames, seed = 42)
            colnames(iso@meta.data)[ncol(iso@meta.data)] <- colnames(iso@meta.data)[ncol(iso@meta.data)] %>%
              stringr::str_sub(.,start = 1,end = -2)

            for (cts in ctlist_to_check){
              message("Running for data:" %+% ident_id %+% " ," %+% "run #" %+% irun %+% ", cells:" %+% cts)
              cts <- cts %>% strsplit(",") %>% unlist
              Sender <- cts[1]; Mediator <- cts[2]; Receiver <- cts[3];
              if (length(intersect(rownames(iso.cellchat@net$weight), c(Sender,Mediator,Receiver))) == 3){
                # sender2mediator.weight <- iso.cellchat@net$weight[Sender,Mediator] # ,W1
                sender2receiver.weight <- iso.cellchat@net$weight[Sender,Receiver] # X,W2
                mediator2receiver.weight <- iso.cellchat@net$weight[Mediator,Receiver] # M,W3
                receiverFuncOfInterest.signal <- median(iso@meta.data[iso@meta.data[,feature.name] == Receiver,funcNames]) # Y
                mat4cma <- rbind(mat4cma,
                                 c(ident_id, irun, paste0(cts, collapse=","),
                                   sender2receiver.weight,
                                   mediator2receiver.weight,
                                   receiverFuncOfInterest.signal))
                colnames(mat4cma) = c("Sample_id","Rep_id","CA_rls",
                                      "X","M","Y") # CA_rls stands for causual mediation relationship
              }
            }
          } else {
            message("int.type, " %+% int.type %+% " not recognizable")
            return(1)
          }
        }
      }
    }
  } else if (grep("iTALK",int.type,value=T) %>% length > 0){

    for (ident_id in unique(so@meta.data[,batch])){
      if (grep(ident_id %+% "_run", list.files(dataDir), value = T) %>% length > 0){
        maxrun = grep(ident_id %+% "_run", list.files(dataDir), value = T) %>%
          stringr::str_replace(ident_id %+% "_run", "") %>%
          stringr::str_replace(".iTALK.rds", "") %>%
          as.numeric %>% max
        for (irun in 1:maxrun){
          iso.iTALK <- loadRDS(file.path(dataDir,ident_id %+% "_run" %+% irun %+% ".iTALK.rds"))

          iso <- subset(so, cells=iso.iTALK$cellnames)
          Seurat::DefaultAssay(iso) <- "RNA"

          if (! is.null(customedSymbolList)){
            ifuncToUse <- lapply(funcToUse, function(x) intersect(rownames(iso), names(customedSymbolList)[customedSymbolList %in% x]))
          } else {
            ifuncToUse <- lapply(funcToUse, function(x) intersect(rownames(iso), x))
          }
          iso <- Seurat::AddModuleScore(iso, features = ifuncToUse, name = funcNames, seed = 42)
          colnames(iso@meta.data)[ncol(iso@meta.data)] <- colnames(iso@meta.data)[ncol(iso@meta.data)] %>%
            stringr::str_sub(.,start = 1,end = -2)

          for (cts in ctlist_to_check){
            message("Running for data:" %+% ident_id %+% " ," %+% "run #" %+% irun %+% ", cells:" %+% cts)
            cts <- cts %>% strsplit(",") %>% unlist
            Sender <- cts[1]; Mediator <- cts[2]; Receiver <- cts[3];
            if (length(intersect(unique(iso.iTALK$scoremat$cell_from,iso.iTALK$scoremat$cell_to), c(Sender,Mediator,Receiver))) == 3){
              # sender2mediator.weight <- iso.cellchat@net$weight[Sender,Mediator] # ,W1
              sender2receiver.weight <- with(iso.iTALK$scoremat, score[cell_from == Sender & cell_to == Receiver]) # X,W2
              mediator2receiver.weight <- with(iso.iTALK$scoremat, score[cell_from == Mediator & cell_to == Receiver]) # M,W3
              receiverFuncOfInterest.signal <- median(iso@meta.data[iso@meta.data[,feature.name] == Receiver,funcNames]) # Y
              mat4cma <- rbind(mat4cma,
                               c(ident_id, irun, paste0(cts, collapse=","),
                                 sender2receiver.weight,
                                 mediator2receiver.weight,
                                 receiverFuncOfInterest.signal))
              colnames(mat4cma) = c("Sample_id","Rep_id","CA_rls",
                                    "X","M","Y") # CA_rls stands for causual mediation relationship
            }
          }
        }
      }
    }
  } else {
    message("The interaction type is illegal, Please use 'iTALK','cellchat.raw','cellchat.p' instead")
    return(1)
  }

  if (! ifSave){
    return(mat4cma)
  } else {
    if (!dir.exists(dirname(outputFile))){
      outputFile = file.path(getwd(),"mat4cma.rds")
    }
    saveRDS(mat4cma, outputFile)
    print("Matrix for causal mediation analysis has saved to " %+% outputFile)
    if (ifReturn){
      return(mat4cma)
    } else {return(NULL)}
  }
}



#' @title run_CausalMeidation
#' @description run causual mediation analysis
#'
#' causal mediation based on specific ligand-receptor interactions
#' deparacated parameter: mod the mode used to measure the signals. spc for specific signals, which only consider the possible interactions between the interested cells. gnl for general signals, which consider all signals of ligands and receptors from the cells.
#'
#' @param mat4mca a matrix generated by generate_matrix_for_CMA using the seurat object
#' @param useColor a boolean value.
#' @param colList a vector of color values. Default is NULL.
#' @param pvalue the threshold value to determine the significance.
#' @param ifPlot a boolean value. Whether to generate the plots for those cells displayed a partial/fully mediation relationship. We only generate plots for those show a consistent direction of mediation effects between Sender-Mediator pair and Mediator-Receiver pair. Default is TRUE.
#' @param ifSavePlot a boolean value. Masked by ifPlot. Only works when ifPlot is True.
#' @param ifCombinePlot a boolean value. Masked by ifSavePlot. Only works when ifSavePlot is True. If many causual mediation relationships detected, we recommend to display the relationships in a combined way. Each combined figure can display up to nine relationships.
#' @param ifReturnPlots a boolean value. Masked by ifPlot. Only works when ifPlot is True. Whether to print the plots in R IDE. If not, save them to saveDir with the prefix, . Otherwise,
#' @param ifReturnFitModels a boolean value.
#' @param ifSaveFitModels a boolean value.
#' @param ifPlotOnlySameDirect a boolean value. Default is TRUE. If true, only draw relationships that have same regulation directions in 'Sender->Mediator' and 'Mediator->Receiver'.
#' @param datalabel a character value. Default is "". To label saved fitted model with a customed name.
#' @param saveDir the path for the directory to save plots and fitted models.
#' @return return NULL or a list of results with fitted models and/or plots of accepted causal mediation relationships.
#' @export
run_CausalMediation = function(mat4mca,
                               useColor=TRUE,
                               colList=NULL,
                               pvalue=1e-3,
                               #mod = "spc",
                               ifPlot = TRUE,
                               ifSavePlot = TRUE,
                               ifCombinePlot = TRUE,
                               ifReturnPlots = FALSE,
                               ifReturnFitModels = FALSE,
                               ifSaveFitModels = TRUE,
                               ifPlotOnlySameDirect = TRUE,
                               datalabel = "",
                               saveDir = NULL){

  packs_to_check = c("combinat","lavaan","stringr","cowplot","extrafont")
  load_packages(packs_to_check)
  
  ## treated annotated cells like CD4+ T cell to CD4_T_cell
  mat4cma$CA_rls <- mat4cma$CA_rls %>% stringr::str_replace_all(" ","_") %>% stringr::str_replace_all("\\+","") %>% stringr::str_replace_all("\\*","")

  tot_ctList = mat4cma$CA_rls %>% strsplit(",") %>% unlist %>% unique

  if (useColor){
    if (is.null(colList) | (length(colList) < length(tot_ctList))){
      message("The color list to use is not provided or the number of provided colors is not enough.\nUse randomly generated color list instead.")
      colList = get_random_colorlist(tot_ctList)
    } else {
      colList = colList[1:length(tot_ctList)]
    }
  }

  pval4sig = pvalue
  # mod = "spc" # mod = "gnl" # # two types of modes, specific (spc) or general (gnl)

  if (is.null(saveDir) | (!dir.exists(saveDir))){
    saveDir = getwd()
  }

  fitList = list()
  fitList_same_sign = list()
  pList = list()

  for (cts in unique(mat4cma$CA_rls)){

    data4cma = mat4cma %>% filter(CA_rls == cts)
    cts = cts %>% strsplit(",") %>% unlist
    sender = cts[1]
    mediator = cts[2]
    receiver = cts[3]

    colnames(data4cma)[colnames(data4cma)=="X"] = sender
    colnames(data4cma)[colnames(data4cma)=="M"] = mediator
    colnames(data4cma)[colnames(data4cma)=="Y"] = receiver

    # if (mod == "spc"){
    #   colnames(data4cma)[colnames(data4cma)=="W1"] = sender
    #   colnames(data4cma)[colnames(data4cma)=="W3"] = mediator
    #   colnames(data4cma)[colnames(data4cma)=="W2"] = receiver
    # } else if (mod == "gnl"){
    #   colnames(data4cma)[colnames(data4cma)=="X"] = sender
    #   colnames(data4cma)[colnames(data4cma)=="M"] = mediator
    #   colnames(data4cma)[colnames(data4cma)=="Y"] = receiver
    # }

    if ((data4cma[,sender] %>% unique %>% length) > 1 & (data4cma[,mediator] %>% unique %>% length) > 1 ){
      model4cma <- '
      ## indirect effect
      ' %+% receiver %+% ' ~ b*' %+% mediator %+% '
      ## mediator
      ' %+% mediator %+% ' ~ a*' %+% sender %+% '
      ## direct effect
      ' %+% receiver %+% ' ~ c*' %+% sender %+% '
      ## indicrect effect (a*b)
         ab := a*b
      ## total effect
         total := c + a*b
      '
      fit<-sem(model=model4cma, data=data4cma, meanstructure=TRUE, std.lv=TRUE, estimator="MLM")
      # bug fixed according to error message 'operator is invalid for atomic vectors', 241218
      fit_summary<- parameterEstimates(fit)#summary(fit,fit.measures=TRUE);fit_summary[["pe"]]
      if (length(which(fit_summary[c(1,2), "pvalue"] < pval4sig)) == 2){
        fitList[[paste_with_sep(cts)]] = fit
        tmp_signs = sign(fit_summary[c(1,2), "est"])
        if (tmp_signs[1]==tmp_signs[2]){
          fitList_same_sign[[paste_with_sep(cts)]] = fit
        }
        if (ifPlot){
          if (ifPlotOnlySameDirect){
            if (tmp_signs[1]==tmp_signs[2]){
              p <- plot_cell_mediation_plot(fit, cols = named_vec(c("X","M1","Y"),
                                                                  colList[cts] %>% as.character),
                                            title_text = "Cell types:\n" %+% paste_with_sep(cts),
                                            base_size = 2.5, base_family = "mono",
                                            title_size = 6)
              pList[[paste_with_sep(cts)]] = p
            }
          } else {
            p <- plot_cell_mediation_plot(fit, cols = named_vec(c("X","M1","Y"),
                                                                colList[cts] %>% as.character),
                                          title_text = "Cell types:\n" %+% paste_with_sep(cts),
                                          base_size = 2.5, base_family = "mono",
                                          title_size = 6)
            pList[[paste_with_sep(cts)]] = p
          }
        }
      } else {print("No satisfied fit model remained.")}
    }
  }
  toReturn = list()
  if (ifSaveFitModels){
    saveRDS(list(fitList = fitList, fitList_same_sign = fitList_same_sign),
            file = file.path(saveDir, "sem_fitted_models." %+% datalabel %+% ".rds"))
  }
  if (ifReturnFitModels){
    toReturn[["fittedModels"]] = list(fitList = fitList, fitList_same_sign = fitList_same_sign)
  }
  if (ifPlot & ifSavePlot){
    if (ifCombinePlot){
      if (length(pList)>9){
        partn = ceiling(length(pList) / 9)
        for (pn in 1:partn){
          if (pn != partn) {
            p_combined = plot_grid(plotlist = pList[(1+(pn-1)*9):(pn*9)])
            savePlot(pObj = p_combined, pDir = saveDir, pType = "med_net", pName = "combined_part" %+% pn, pFmt = "pdf")
          } else {
            pnum = length(pList)-(1+(pn-1)*9)+1
            prow = ceiling(pnum/3)
            p_combined = plot_grid(plotlist = pList[(1+(pn-1)*9):length(pList)], ncol = 3)
            savePlot(pObj = p_combined, pDir = saveDir, pType = "med_net", pName = "combined_part" %+% pn, pFmt = "pdf", pHt = 7/3*prow)
          }
        }
      } else if (length(pList) > 0) {
        p_combined = plot_grid(plotlist = pList)
        savePlot(pObj = p_combined, pDir = saveDir, pType = "med_net", pName = "combined_all", pFmt = "pdf")
      }
    } else {
      if (length(pList) > 0){
        for (i in 1:length(pList)){
          savePlot(pObj = pList[[i]], pDir = saveDir, pType = "med_net", pName = names(pList)[i], pFmt = "pdf", pHt = 2.5, pWd = 2.5)
        }
      }
    }
  } else if (ifPlot){
    if (ifReturnPlots & length(pList) > 0){
      toReturn[["plots"]] = pList
    }
  }

  if (length(toReturn) > 0){
    print("Causual mediation analysis has finished.")
    if (ifReturnPlots){print("Plots stored in plots.")}
    if (ifReturnFitModels){print("Fitted models passed the significance test were stored in fittedModels.")}
    return(toReturn)
  } else {
    print("Causual mediation analysis has finished. Results were stored in " %+% saveDir)
  }
}






