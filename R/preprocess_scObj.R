
#' @title integrateSO
#' @description integrated seurat object
#'
#' This is a customed one-step function wrapped the integration pipeline based on functions imported from Seurat and harmony.
#'
#' @param solist a list of seurat objects to integrate
#' @return a combined seurat object
#'
integrateSO <- function(solist){

  packs_to_check = c("Seurat", "R.filesets", "ggplot2", "ggsci", "jzhang.utils", "dplyr", "harmony")
  load_packages(packs_to_check)

  so_combined = solist
  return(so_combined)
}



#' @title bootstrapSO
#' @description bootstrap seurat object
#'
#' @param iso A seurat object
#' @param feature.name The name of the feature, which come from the metatable of seurat object, to classfiy the cells
#' @param minCNUM.tot The minimum cell number for a cluster, cluster with cells less the minCNUM.tot will be removed from the analysis. Default 10.
#' @param sampleProp The porportion of cells used in resampling for each subcluster. Default 0.1.
#' @param minCNUM.new Default is 10. The expected minimum cells for resmple. If estimated cells used for resample, which is the cluster_cell_size*sampleProp, is smaller than the minimum cell number for resampled cells in each subcluster, use minimum cell number as the resample size. This value should be not greater than minCNUM.tot to avoid error.
#' @param runs Total runs for the bootstrap. Default is 5.
#' @param suggestedRunTimes default is False. If is True, the bootstrap won't run. The function will estimate the maximum runs for non-redundancy sampling results.
#' @param returnNameOnly default is True. If is True, only return the cellnames rather than the resampled objects.
#' @return a seurat object list or the cell names of a seurat object list or the estimated maximum run times.
#' @export
bootstrapSO <- function(iso, feature.name = "seurat_clusters", minCNUM.tot = 10, sampleProp = 0.1, minCNUM.new = 10, runs = 5, suggestedRunTimes = F, returnNameOnly = T){
  ## filter groups contains cells with cell number < minCNUM

  packs_to_check = c("Seurat", "R.filesets", "ggplot2", "ggsci", "jzhang.utils", "dplyr", "harmony")
  load_packages(packs_to_check)

  failed_ctlist = table(iso@meta.data[,feature.name]) %>% .[. < minCNUM.tot] %>% names
  if (length(failed_ctlist) < uniqlen(iso@meta.data[,feature.name])){
    if (length(failed_ctlist)!=0){
      iso = subset(iso, cells = Cells(iso)[! iso@meta.data[,feature.name] %in% failed_ctlist])
    }
    expected.cells.size = table(iso@meta.data[,feature.name])*sampleProp
    ## output suggestedRunTimes
    if (suggestedRunTimes){
      realProp = minCNUM.new/expected.cells.size
      maxTimes = prod(floor(expected.cells.size/minCNUM.new) %>% .[. > 0])
      print("The maximum times estimated for bootstrap: " %+% maxTimes)
      print("The real proportion for each clusters: " %+% paste0(paste(names(realProp),realProp,sep=":"), collapse=", "))
      return(maxTimes)
    } else {
      set.seed(round(runif(1,0,16666661)))
      cellnames = list()
      for (irun in 1:runs){
        cellnames[[irun]] = 0
      }
      for (ct in names(expected.cells.size)){
        if (expected.cells.size[ct] > minCNUM.new){
          for (irun in 1:runs){
            cellnames[[irun]] = c(cellnames[[irun]], sample(Cells(iso)[iso@meta.data[,feature.name]==ct], expected.cells.size[ct], replace = F))
          }
        } else if (minCNUM.new < expected.cells.size[ct]/sampleProp){
          for (irun in 1:runs){
            cellnames[[irun]] = c(cellnames[[irun]], sample(Cells(iso)[iso@meta.data[,feature.name]==ct], minCNUM.new, replace = F))
          }
        } else {
          for (irun in 1:runs){
            cellnames[[irun]] = c(cellnames[[irun]], Cells(iso)[iso@meta.data[,feature.name]==ct])
          }
        }
      }
      solist = list()
      for (irun in 1:runs){
        iso.resampled = subset(iso, cells = cellnames[[irun]][-1])
        solist[[irun]] = iso.resampled
      }
      if (returnNameOnly){
        solist.cellnames = lapply(solist, Seurat::Cells)
        return(solist.cellnames)
      } else {return(solist)}
    }
  } else {
    print(failed_ctlist %>% names %>% paste0(.,collapse=",") %+% " all failed")
    print("This Seurat Object not suitable for bootstrap, or parameters require adjustments.")
  }
}
