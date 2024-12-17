
#' @title load_file
#' @description load table files
#' 
#' @param filepath file path of the data table
#' @param delimiter Default is NULL.
#' @return data table
#' @export
load_file = function(filepath, delimiter=NULL, ...){

  packs_to_check = c("Seurat","harmony","umap", "R.filesets")
  load_packages(packs_to_check)

  if (! file.exists(filepath)){
    message(filepath %+% " not exists.")
    break
  } else {
    ext = strsplit(filepath, "\\.") %>% unlist %>% tail(., n=1)
    if (ext %in% c("txt","tsv","csv")){
      if (is.null(delimiter)){
        if (ext == "csv"){
          delimiter = ","
        } else {delimiter = "\t"}
      }
      dat = read.delim(filepath, sep = delimiter, ...)
    } else if (ext == "rds"){
      dat = loadRDS(filepath)
    }
    return(dat)
  }
}



#' @title load_data
#' @description import the single cell data and transfer data to seurat object.
#' 
#' @param dataDir the directory path for single cell data
#' @param dataType the format of the raw data. Accepted values are 10x, h5, rds.
#' @param metaPath the filepath of the meta file. Default is NULL.
#' @return a seurat object
#' @export
load_data = function(dataDir, dataType, metaPath = NULL){

  packs_to_check = c("Seurat","harmony","umap", "R.filesets")
  load_packages(packs_to_check)


  if (!is.null(metaPath)){
    meta = load_file(metaPath)
  }
  if (exists("meta", 2, inherits = FALSE)){
    if (!(ncol(so)==nrow(meta))){
      meta = NULL
    }
  } else {meta = NULL}

  if (dataType == "10x"){
    so = Seurat::Read10X(data.dir = dataDir)
    so = Seurat::CreateSeuratObject(so)
    so@meta.data = cbind(so@meta.data, meta)
  } else if (dataType == "h5"){
    so = Seurat::Read10X_h5(filename = grep("h5",list.file(dataDir), value=T))
    if (!is.null(meta)){so@meta.data = cbind(so@meta.data, meta)}
  } else if (dataType == "rds"){
    so = loadRDS(filename = grep("rds",list.file(dataDir), value=T))
    if (!is.null(meta)){so@meta.data = cbind(so@meta.data, meta)}
  } else {
    print(dataType %+% " is not identified. Please use 10x, h5, or rds")
  }
  return(so)
}



#' @title load_so
#' @description load newest seurat object
#' 
#' @param so_dir The directory of stored seurat object in rds format
#' @return seurat object
#' @export
load_so = function(so_dir){

  packs_to_check = c("Seurat","harmony","umap", "R.filesets")
  load_packages(packs_to_check)

  if (!dir.exists(so_dir)){
    print(so_dir %+% " does not exist!")
    break
  }
  soFiles = list.files(so_dir) %>% grep("rds", ., value=T)
  meta = grep("so.meta", soFiles, value = T)
  so = grep("so.[0-9]+.rds", soFiles, value  =T)
  if (length(so) == 1){
    so = loadRDS(file.path(so_dir, so))
  } else if (length(so) > 0){
    newest_date = so %>% stringr::str_replace(., "so\\.", "") %>% stringr::str_replace(., "\\.rds", "") %>% as.numeric %>% max
    so = loadRDS(file.path(so_dir, "so." %+% newest_date %+% ".rds"))
  } else {
    "Seurat object not find."
    break
  }
  if (length(meta) == 1){
    so@meta.data = loadRDS(file.path(so_dir, meta))
  } else if (length(meta) > 0){
    newest_date = meta %>% stringr::str_replace(., "so\\.meta\\.", "") %>% stringr::str_replace(., "\\.rds", "") %>% as.numeric %>% max
    so@meta.data = loadRDS(file.path(so_dir, "so.meta." %+% newest_date %+% ".rds"))
  }
  return(so)
}



#' @title update_so
#' @description save or update seurat object and its meta.data
#' 
#' @param so seurat object
#' @param so_dir the seurat directory
#' @param meta.only a boolean variable. Default is TRUE. If TRUE, only update metatable, seurat object will not be updated.
#' @export
update_so = function(so, so_dir, meta.only = T){

  packs_to_check = c("Seurat","harmony","umap", "R.filesets")
  load_packages(packs_to_check)

  if (!dir.exists(so_dir)){
    # print(so_dir %+% " does not exist!")
    createDir(so_dir)
  }
  soFiles = list.files(so_dir) %>% grep("rds", ., value=T)
  metaPath = grep("so.meta", soFiles, value = T)
  soPath = grep("so.[0-9]+.rds", soFiles, value  =T)
  dstamp = getdstamp()
  if (meta.only){
    if (length(metaPath) > 0){
      newest_date = metaPath %>% stringr::str_replace(., "so\\.meta\\.", "") %>% stringr::str_replace(., "\\.rds", "") %>% as.numeric %>% max
      if (newest_date == dstamp){
        print("Old meta file will be overwritten.")
      }
    }
    saveRDS(so@meta.data, file.path(so_dir, "so.meta." %+% dstamp %+% ".rds"))
  }
  else {
    if (length(metaPath) > 0){
      newest_date = metaPath %>% stringr::str_replace(., "so\\.meta\\.", "") %>% stringr::str_replace(., "\\.rds", "") %>% as.numeric %>% max
      if (newest_date == dstamp){
        print("Old meta file will be overwritten.")
      }
    }
    saveRDS(so@meta.data, file.path(so_dir, "so.meta." %+% dstamp %+% ".rds"))

    if (length(soPath) > 0){
      newest_date = soPath %>% stringr::str_replace(., "so\\.", "") %>% stringr::str_replace(., "\\.rds", "") %>% as.numeric %>% max
      if (newest_date == dstamp){
        print("Old seurat object will be overwritten.")
      }
    }
    saveRDS(so, file.path(so_dir, "so." %+% dstamp %+% ".rds"))

  }
}



#' @title separate_so
#' @description separate seurat object by batch feature
#' 
#' @return a list of seurat subobjects
#' @export
separate_so = function(so, batch = "orig.ident"){

  packs_to_check = c("Seurat","harmony","umap", "R.filesets", "dplyr")
  load_packages(packs_to_check)

  if (batch %in% colnames(so@meta.data)){
    solist = list()
    sublist = so@meta.data[,batch] %>% unique %>% as.character
    for (sb in sublist){
      solist[[sb]] = subset(so, cells = Cells(so)[so@meta.data[,batch] %>% as.character == sb])
    }
    return(solist)
  } else {
    print("The given batch, " %+% batch %+% " is not identified in the seurat object.")
  }
}



#' @title run_seurat
#' @description a customed pipeline for the preprocessing of seurat object.
#' 
#' @param so a seurat object.
#' @param label the annotation of seurat obejct.
#' @param batch the variable name of batch.
#' @param features_to_view default is NULL. a list of features for the visualization of feature distribution in the UMAP.
#' @param save a boolean value. Whether to save treated seurat object to so_dir. Default is FALSE.
#' @param so_dir default is NULL. The directory to store treated seurat object.
#' @param prior_metadata default is NULL. a customed meta table.
#' @param treat.cellname a boolean value. If TRUE, add label to seurat cell names.
#' @param mtlabel default is "MT-". You can define a customed mitochondrial gene list by yourself.
#' @param min_avg default is 0.01. Minimum expressed average gene.
#' @param assay default is "RNA".
#' @param scale_factor default is 1e+06 for CPM.
#' @param min_bin default is 50. Minimum bin value.
#' @param min_cutoff default is c(0.01, Inf).
#' @param dimnum default is "auto".
run_seurat = function (so, label = "", batch = NULL, features_to_view = NULL,
                       save = T, so_dir = NULL, prior_metadata = NULL, addNamedColorList = NULL,
                       treat.cellname = F, mtlabel = "MT-", min_avg = 0.01, assay = "RNA",
                       scale_factor = 1e+06, min_bin = 50, min_cutoff = c(0.01, Inf),
                       calcClustMarker = T,
                       dimnum = "auto", dim_sd = 1.25, dim_pval = 0.01, clustRes = 0.8) {

  packs_to_check = c("Seurat","harmony","umap", "R.filesets")
  load_packages(packs_to_check)

  if (dimnum == "auto"){
    ## dimnum range from 5-80
    dimnum = max(min(floor(ncol(so)/200),80),5)
  }
  print("Default dims to use: " %+% dimnum)

  # bug fixed 20241217
  if (length(mtlabel) > 0) { 
    mtlist = rownames(so)[grep(mtlabel, toupper(rownames(so)))]
  } else {
    mtlist = mtlabel
  }

  opath4seurat = "."
  if (is.null(so_dir)) {
    so_dir = file.path(opath4seurat, label)
  } else {
    so_dir = file.path(so_dir, label)
  }
  createDir(so_dir)
  if (treat.cellname) {
    so = Seurat::RenameCells(so, add.cell.id = label)
  }
  all.genes = rownames(so)
  # bug fixed 20241217
  if (as.numeric(packageVersion('Seurat')[1,1]) >= 5) {
    hiexpr.genes = apply(Seurat::GetAssayData(so, assay = "RNA", layer = "counts"), 1, mean)
  } else {
    hiexpr.genes = apply(Seurat::GetAssayData(so, assay = "RNA", slot = "counts"), 1, mean)
  }
  hiexpr.genes = apply(so@assays$RNA@counts, 1, mean)
  so <- Seurat::PercentageFeatureSet(so, features = intersect(mtlist,
                                                              rownames(so)), col.name = "percent.mt")
  so <- Seurat::NormalizeData(so, normalization.method = "LogNormalize",
                              scale.factor = scale_factor, assay = assay)
  # options(future.globals.maxSize = 891289600)
  # plan(strategy = "multicore", workers = length(unique(so$orig.ident)))
  so <- Seurat::ScaleData(so, vars.to.regress = batch, assay = assay,
                          features = names(hiexpr.genes)[hiexpr.genes > min_avg])
  so <- Seurat::FindVariableFeatures(so, selection.method = "mean.var.plot",
                                     num.bin = min_bin, mean.cutoff = min_cutoff, assay = assay)
  gc()
  so <- Seurat::RunPCA(so, assay = assay, npcs = dimnum, verbose = FALSE,
                       reduction.name = "pca", reduction.key = assay)
  p_elbowplot_rnapc = Seurat::ElbowPlot(so, ndims = dimnum,
                                        reduction = "pca")
  dimnum = with(p_elbowplot_rnapc$data, dims[stdev >= dim_sd])
  pdf(file = so_dir %+% "/p_elbow_rnapc.pdf", width = 4, height = 3)
  print(p_elbowplot_rnapc)
  dev.off()
  gc()
  so <- Seurat::JackStraw(so, num.replicate = 5, reduction = "pca",
                          dims = max(dimnum))
  so <- Seurat::ScoreJackStraw(so, reduction = "pca", dims = dimnum,
                               score.thresh = dim_pval)
  p_jackstrawplot_pca = Seurat::JackStrawPlot(so, dims = dimnum)
  pdf(file = so_dir %+% "/p_jackstraw_pca.pdf", width = 8.6, height = 5.4)
  print(p_jackstrawplot_pca)
  dev.off()
  selected_pca_dims = which(so@reductions$pca@jackstraw$overall.p.values[,2] < dim_pval)
  gc()
  if (save) {
    dstamp = getdstamp()
    saveRDS(so, file.path(so_dir, "so." %+% dstamp %+% ".rds"))
  }
  so <- Seurat::RunUMAP(so, dims = selected_pca_dims, verbose = FALSE, seed.use = 0)
  # bug fixed 20241217
  so <- Seurat::FindNeighbors(so, dims = selected_pca_dims, compute.SNN = TRUE, prune.SNN = 0)
  so <- Seurat::FindClusters(so, resolution = clustRes, verbose = F, group.singletons = T, n.iter = 10, random.seed = 0)
  colorlist = get_random_colorlist(so@meta.data[, "seurat_clusters"])
  colorlist = c(addNamedColorList, colorlist)
  p_dimplot <- Seurat::DimPlot(so, group.by = "seurat_clusters", pt.size = 2.5, reduction = "umap", label = TRUE, raster = TRUE) +
    ggplot2::scale_color_manual(values = colorlist) +
    Seurat::NoLegend() + Seurat::NoAxes()
  savePlot(pObj = p_dimplot, pDir = so_dir, pType = "umap", pName = "res0.8", pFmt = "pdf",
         pHt = 5.2, pWd = 5)
  p_umi <- Seurat::FeaturePlot(so, feature = "nCount_" %+% assay, pt.size = 2.5, order = TRUE, raster = TRUE) + Seurat::NoAxes()
  savePlot(pObj = p_umi, pDir = so_dir, pType = "umap", pName = "seqDepth", pFmt = "pdf", pWd = 5.2,
         pHt = 5)
  gc()
  if (save) {
    saveRDS(so@meta.data, file.path(so_dir, "so.meta." %+% getdstamp() %+% ".rds"))
  }
  if (!is.null(prior_metadata)) {
    for (ifeature in setdiff(grep("Main", colnames(prior_metadata), value = T),
                             c("nCount_" %+% assay, "nFeature_" %+% assay))) {
      colorlist = get_random_colorlist(prior_metadata[, ifeature])
      colorlist = c(addNamedColorList, colorlist)
      if (uniqlen(prior_metadata[, ifeature]) > 1 & uniqlen(prior_metadata[, ifeature]) < 80) {
        p = Seurat::DimPlot(so, group.by = ifeature, pt.size = 2.5, reduction = "umap", label = TRUE, raster = TRUE) +
          ggplot2::scale_color_manual(values = colorlist) +
          Seurat::NoAxes()
        savePlot(pObj = p, pType = "umap", pName = ifeature, pDir = so_dir, pFmt = "pdf", pWd = 6.5, pHt = 5)
      }
    }
  } else {
    if (!is.null(features_to_view)){
      for (ifeature in features_to_view) {
        print(ifeature)
        colorlist = get_random_colorlist(ifeature)
        colorlist = c(addNamedColorList, colorlist)
        if (uniqlen(so@meta.data[,ifeature]) > 1 & uniqlen(so@meta.data[, ifeature]) < 80) {
          p = Seurat::DimPlot(so, group.by = ifeature, pt.size = 2.5, reduction = "umap", label = TRUE, raster = TRUE) +
            ggplot2::scale_color_manual(values = colorlist) +
            Seurat::NoAxes()
          savePlot(pObj = p, pType = "umap", pName = ifeature, pDir = so_dir, pFmt = "pdf", pWd = 6.5, pHt = 5)
        }
      }
    }
  }
  if (calcClustMarker){
    tryCatch(
      {so <- Seurat::BuildClusterTree(so, assay = assay, reduction = "umap")},
      error = function(cond){message(cond);message("\nBuildCluterTree failed in " %+% label);return(NA)},
      warning = function(cond){message(cond);message("\nBuildCluterTree failed in " %+% label);return(NA)})
    all.markers <- Seurat::FindAllMarkers(so, assay = assay, random.seed = 19061391)
    marker_filepath = so_dir %+% "/t_markers.cluster.csv"
    write.csv(all.markers, marker_filepath, quote = F)    
  }
  gc()
  if (!save) {
    return(so)
  }
  else {
    return(NULL)
  }
}



#' @title run_combined_seurat
#' @description run combined seurat object
#'
#' @param solist_dir the path of directory to the directories of stored seurat objects
#' @param dimnum the maximum dimension to use
#' @param save_dir the output directory of integrated object.
#' @param filterRules default is NULL. Use the rules to filter seurat objects before the integration. Using a string like 'CellType %in% c("EC","T_cell","Cancer") & stringr::str_detect(.$CellTypeStatus, "CT_confirmed") & nCount_RNA > 500 & nFeature_RNA > 200' to filter the result.
#' @param acceptRules default is NULL. Use the rules to determine whether to include a filtered seurat object or not. Syntax only work with 'CellType,>,50' to determine.
#'
#' @export
#'
run_combined_seurat = function(solist_dir, save_dir = NULL, ifReturn = FALSE,
                               filterRules = NULL, acceptRules = NULL,
                               dimnum = 50, k.anchor = 5, k.filter = 20, k.weight = 10){

  packs_to_check = c("Seurat","harmony","umap", "R.filesets")
  load_packages(packs_to_check)

  if (is.null(save_dir)){
    save_dir = file.path(dirname(solist_dir), "so_integrated")
  } else if (! dir.exists(save_dir)) {
    createDir(save_dir)
  }

  if (! dir.exists(solist_dir)){
    print("Directory for seurat objects not found:\n" %+% solist_dir)
    break
  }
  solist_dir = normalizePath(solist_dir)

  samples = list.dirs(solist_dir, recursive = F)
  names(samples) = stringr::str_replace(samples, solist_dir %+% "/", "")
  solist = list()
  for (sp in 1:length(samples)){
    tmp_so = load_so(so_dir = samples[[sp]])
    if (! is.null(filterRules)){
      tmp_cellnames = tmp_so@meta.data %>% dplyr::filter(rlang::eval_tidy(rlang::parse_expr(filterRules))) %>% rownames
      if (length(tmp_cellnames) > 0){
        tmp_so = subset(tmp_so, cells = tmp_cellnames)
        if (! is.null(acceptRules)){
          tmp_rules = strsplit(acceptRules, ",") %>% unlist()
          tmp_0 = tmp_so@meta.data[,tmp_rules[1]] %>% table
          tmp_1 = eval(parse(text=paste_with_sep(c('tmp_0[tmp_0',tmp_rules[2],tmp_rules[3],']'),' ')))
          if (length(tmp_0)==length(tmp_1) & length(tmp_0) > 0){
            solist[[names(samples)[sp]]] = tmp_so
            update_so(so = tmp_so,
                      so_dir = file.path(save_dir,#dirname(samples[[sp]]),
                                                      "filtered_" %+% basename(samples[[sp]])),
                      meta.only = F)
          }
        } else {solist[[names(samples)[sp]]] = tmp_so}
      }
    } else {solist[[names(samples)[sp]]] = tmp_so}
  }

  full.genes = rownames(solist[[1]])
  cnums = sapply(solist, ncol) %>% sort(decreasing=TRUE)
  if (min(cnums) < k.filter){
    print("k.filter " %+% k.filter %+% " is greater than so object's cell number " %+% min(cnums))
    print("seurat objects with cell numbers less than k.filter will be removed.")
    cnums = cnums[cnums > k.filter]
    solist = solist[names(cnums)]
  }
  if (length(solist) > 1){
    for (i in 2:length(solist)){
      print("try to merge " %+% names(solist)[i])
      if (i == 2){
        merged = solist[[names(cnums)[1]]]
      }
      merged.anchor = FindIntegrationAnchors(list(merged,solist[[names(cnums)[i]]]),
                                             dims = 1:dimnum,
                                             k.anchor = k.anchor,
                                             k.filter = k.filter)
      merged =IntegrateData(merged.anchor, k.weight = k.weight)
      if (i %% 5 == 0){
        update_so(so = merged, so_dir = save_dir, meta.only = F)
      }
      rm(merged.anchor)
      gc()
    }
    rm(solist)
    gc()
    ## BT1240 abandomed

    var.genes.integrated = merged@assays$integrated@var.features
    DefaultAssay(merged) = "RNA"
    full.genes = rownames(merged)

    update_so(so = merged, so_dir = save_dir, meta.only = F)

    DefaultAssay(merged) = "RNA"
    merged <- Seurat::ScaleData(merged)
    features_for_pca = setdiff(var.genes.integrated,rownames(merged)[apply(GetAssayData(merged, slot = "counts", assay = "RNA"),1,sum) < ncol(merged)*0.0005])
    merged <- Seurat::RunPCA(merged, npcs = dimnum, features = features_for_pca)
    p_elbowplot_rnapc = ElbowPlot(merged, ndims = dimnum, reduction="pca")

    selected_pca_dims = with(p_elbowplot_rnapc$data, dims[stdev >= 1.25])
    print(selected_pca_dims)
    merged <- RunHarmony(merged, group.by.vars = c("orig.ident"), reduction = "pca", dims.use = selected_pca_dims)
    merged <- RunUMAP(merged, dims=selected_pca_dims, reduction = "harmony", verbose=F, seed.use=0)
    merged <- RunTSNE(merged, dims=selected_pca_dims, reduction = "harmony", check_duplicates = FALSE)
    merged <- FindNeighbors(merged, dims=selected_pca_dims, compute.SNN=T, k.param = 10, verbose=F, force.recalc=T)
    merged <- FindClusters(merged, resolution=1, verbose=F)

    update_so(so = merged, so_dir = save_dir, meta.only = F)
    print("Merged seurat object saved to " %+% save_dir)

    if (ifReturn){
      return(merged)
    }
  }
}



#' @title estimate_pretreated_so_by_refMarkers estimate pretreated seurat object by given reference markers
#'
#' @description write a function to test consistency between known annotated labels and seurat clusters
#' a reference cell markers for known cell type labels
#' For those cell types with several seurat clusters, to check:
#' 1. Whether those subclusters all express reference markers
#' 2. Is there any significant novel markers for the specific cluster which were overlapped between distinct samples.
#' Warning clusterName should contains elements the same aligned to those used in the cluster marker file
#'
#' @param so seurat object
#' @param refList a given reference cell type markers
#' @param clusterMarkerPath the path of calculated markers generated by seurat find_markers
#' @param tissue The tissue name, used to identify markers in same cell types but from different tissues, such as Epithelial.CRC and Epithlial.LC
#' @param annotateName The feature name for given cell type
#' @param clusterName The feature name for cluster
#' @return a cluster-marker information table, cluster as subclusters for given cell types with novel markers or cluster contains several given cell types
#' @export
estimate_pretreated_so_by_refMarkers = function(so, refList, clusterMarkerPath, tissue = NULL, annotateName = "CellType", clusterName = "seurat_clusters"){

  packs_to_check = c("Seurat","harmony","umap", "R.filesets")
  load_packages(packs_to_check)

  markers = load_file(clusterMarkerPath)
  seurat_cluster_per_CellType = split(so@meta.data[,clusterName], so@meta.data[,annotateName])
  final_signal = data.frame() # ct, cluster, avg.pos, avg.neg
  for (ct in names(seurat_cluster_per_CellType)){
    tmp.clusters.num = table(seurat_cluster_per_CellType[[ct]]) %>% .[. > 10]
    if (length(tmp.clusters.num) > 0){
      tmp.cellnames = Cells(so)[so@meta.data[,annotateName] == ct & so@meta.data[,clusterName] %in% names(tmp.clusters.num)]
      if (is.null(refList[[ct]])){
        if (is.null(refList[[ct %+% "." %+% tissue]])){
          print(ct %+% " markers are not found. Skip...")
        } else {
          tmp.markers = refList[[ct]]
        }
      } else {
        tmp.markers = refList[[ct]]
      }
      tmp.markers = intersect(tmp.markers %>% toupper, rownames(so) %>% toupper)
      tmp.mlen = length(tmp.markers)
      if (tmp.mlen > 0){
        if (tmp.clusters.num %>% length == 1){
          tmp.markers.expr = AverageExpression(subset(so, cells = tmp.cellnames), features = tmp.markers, assays = "RNA", slot = "counts")
          tmp.markers.expr.ctrl = AverageExpression(subset(so, cells = setdiff(Cells(so),tmp.cellnames)), features = tmp.markers, assays = "RNA", slot = "counts")
          tmp.res = data.frame(ct = rep(ct, tmp.mlen),
                               gene = tmp.markers, cluster = rep(names(tmp.clusters.num), tmp.mlen),
                               pos = tmp.markers.expr$RNA[,1],
                               neg = tmp.markers.expr.ctrl$RNA[,1])
          tmp.res$cnum = tmp.clusters.num[tmp.res$cluster]
          rownames(tmp.res) = NULL
          if (nrow(final_signal) == 0){
            final_signal = tmp.res
          } else {
            final_signal = rbind(final_signal, tmp.res)
          }
        } else {
          tmp.markers.expr = AverageExpression(subset(so, cells = tmp.cellnames), features = tmp.markers, group.by = clusterName, assays = "RNA", slot = "counts")
          tmp.markers.expr.ctrl = AverageExpression(subset(so, cells = setdiff(Cells(so),tmp.cellnames)), features = tmp.markers, assays = "RNA", slot = "counts")
          tmp.res = data.frame(ct = rep(ct, tmp.mlen*length(tmp.clusters.num)),
                               gene = rep(tmp.markers, length(tmp.clusters.num)),
                               cluster = rep(names(tmp.clusters.num), each = tmp.mlen),
                               pos = reshape2::melt(tmp.markers.expr$RNA)$value,
                               neg = rep(tmp.markers.expr.ctrl$RNA[,1], length(tmp.clusters.num)))
          tmp.res$cnum = tmp.clusters.num[tmp.res$cluster]
          rownames(tmp.res) = NULL
          if (nrow(final_signal) == 0){
            final_signal = tmp.res
          } else {
            final_signal = rbind(final_signal, tmp.res)
          }
        }
      }
    }
  }

  ## check whether reference marker genes can be used to identify cell type in subclusters
  suspecious_clusters = detect_suspecious_clusters(final_signal)

  if (length(suspecious_clusters$confirmed) > 0){
    ## for cell types with only one confirmed subclusters
    tmp.cts = suspecious_clusters$confirmed %>% str_split(",", simplify = T) %>% .[,1] %>% table %>% .[. == 1] %>% names
    if (length(tmp.cts) > 0){
      tmp.clusters = suspecious_clusters$confirmed %>% str_split(",", simplify=T) %>% as.data.frame %>% filter(V1 %in% tmp.cts) %>% pull(V2)
      tmp.markers = data.frame(cluster = tmp.clusters, newMark = "", ct = tmp.cts, status = "CT_confirmed_refMark")
      if (nrow(tmp.markers) > 0){
        marker_info = tmp.markers
      }
    }
    ## for cell types with several confirmed subclusters
    tmp.cts = suspecious_clusters$confirmed %>% str_split(",", simplify = T) %>% .[,1] %>% table %>% .[. > 1] %>% names
    if (length(tmp.cts) > 0){
      tmp.clusters = suspecious_clusters$confirmed %>% str_split(",", simplify=T) %>% as.data.frame %>% filter(V1 %in% tmp.cts) %>% pull(V2)
      tmp.cts = sapply(tmp.clusters, function(x) grep("," %+% x %+% "$", suspecious_clusters$confirmed, value=T) %>% stringr::str_replace("," %+% x, "") %>% paste_with_sep)
      tmp.markers = markers %>%
        dplyr::mutate(pct.diff = pct.1-pct.2) %>%
        dplyr::filter(cluster %in% tmp.clusters, pct.diff > 0.6, pct.2 < 0.2) %>%
        dplyr::arrange(cluster, -pct.diff) %>%
        group_by(cluster) %>% top_n(n=3) %>%
        summarise(paste_with_sep(gene)) %>%
        mutate(ct = tmp.cts[as.character(cluster)], status = "CT_confirmed_newSubMark")
      colnames(tmp.markers)[2] = "newMark"
      if (nrow(tmp.markers) > 0){
        if (exists("marker_info")){
          marker_info = rbind(marker_info, tmp.markers)
        } else {marker_info = tmp.markers}
      }
    }
  }

  ## for clusters without confirmed types
  if (length(suspecious_clusters$unconfirmed) > 0){
    tmp.clusters = suspecious_clusters$unconfirmed %>% str_split(",", simplify=T) %>% .[,2]
    if (length(tmp.clusters) > 0){
      tmp.cts = sapply(tmp.clusters, function(x) grep("," %+% x %+% "$", suspecious_clusters$unconfirmed, value=T) %>% stringr::str_replace("," %+% x, "") %>% paste_with_sep)
      tmp.markers = markers %>%
        dplyr::mutate(pct.diff = pct.1-pct.2) %>%
        dplyr::filter(cluster %in% tmp.clusters, pct.diff > 0.6, pct.2 < 0.2) %>%
        dplyr::arrange(cluster, -pct.diff) %>%
        group_by(cluster) %>% top_n(n=3) %>%
        summarise(paste_with_sep(gene)) %>%
        mutate(ct = tmp.cts[as.character(cluster)], status = "CT_unconfirmed_newSubMark")
      colnames(tmp.markers)[2] = "newMark"
      if (nrow(tmp.markers) > 0){
        if (exists("marker_info")){
          marker_info = rbind(marker_info, tmp.markers)
        } else {marker_info = tmp.markers}
      }
    }
  }


  ## for clusters with mixed raw annotation
  tmp.markers = so@meta.data[so@meta.data[,clusterName] %in% setdiff(so@meta.data[,clusterName] %>% unique, marker_info$cluster),c(annotateName, clusterName)] %>%
    distinct %>% group_by(eval(parse(text=clusterName))) %>% summarise(paste_with_sep(eval(parse(text=annotateName))))
  if (length(tmp.markers) > 0){
    colnames(tmp.markers) = c("cluster", "ct")
    tmp.cts = named_vec(tmp.markers$cluster, tmp.markers$ct)
    tmp.markers = markers %>%
      dplyr::mutate(pct.diff = pct.1-pct.2) %>%
      dplyr::filter(cluster %in% names(tmp.cts), pct.diff > 0.6, pct.2 < 0.2) %>%
      dplyr::arrange(cluster, -pct.diff) %>%
      group_by(cluster) %>% top_n(n=3) %>%
      summarise(paste_with_sep(gene)) %>%
      mutate(ct = tmp.cts[as.character(cluster)], status = "CT_unconfirmedMixed_newSubMark")
    colnames(tmp.markers)[2] = "newMark"
    if(length(tmp.cts) > nrow(tmp.markers)){
      tmp.markers = rbind(tmp.markers,
                          data.frame(cluster = setdiff(names(tmp.cts), tmp.markers$cluster),
                                     newMark = "",
                                     ct = tmp.cts[setdiff(names(tmp.cts), tmp.markers$cluster)],
                                     status = "CT_unconfirmedMixed_NoMark"))
    }
    if (nrow(tmp.markers) > 0){
      if (exists("marker_info")){
        marker_info = rbind(marker_info, tmp.markers)
      } else {marker_info = tmp.markers}
    }
  }
 if (nrow(marker_info) > 0){
   return(marker_info)
 } else {print("Generating table for marker info failed.")}
}



#' @title detect_suspecious_clusters
#' @description Detecting suspecious clusters
#' 
#' @param final_signal_table a table per row shows: ct, gene, cluster, pos (gene expression in target cluster), neg (gene expression in control cells), cnum (cell number detected in the cluster)
#' @return a vector of elements shows cell type and cluster glued by comma
#' @export
detect_suspecious_clusters = function(final_signal_table){

  packs_to_check = c("Seurat","harmony","umap", "R.filesets")
  load_packages(packs_to_check)

  clusters_to_test = final_signal_table %>% dplyr::select(ct, cluster) %>% distinct %>% apply(., 1, paste_with_sep)
  clusters_confirmed = final_signal_table %>% mutate(diff = pos - neg, diff.prop=abs(pos - neg)/max(pos, neg)) %>% dplyr::filter(diff >= 1 | diff.prop >= 0.1) %>% dplyr::select(ct, cluster) %>% distinct %>% apply(., 1, paste_with_sep)
  clusters_unconfirmed = setdiff(clusters_to_test, clusters_confirmed)
  return(list(confirmed = clusters_confirmed, unconfirmed = clusters_unconfirmed))
}
