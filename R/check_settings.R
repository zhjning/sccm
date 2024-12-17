## Check Package Availability of Dependent Packages ----

# Define the list of packages to check
packs_to_check = c("lavaan", "plyr", "dplyr", "extrafont", "R.filesets", "ggplot2", "ggsci", "Seurat", "umap", "combinat","stringr","msigdbr","harmony","usethis")

# Get the list of currently installed packages
installed_packs = installed.packages() 

# Define custom infix operators for string concatenation
`%+%` <- paste0 # Concatenate strings without any separator
`%++%` <- function (...) {paste(sep = ".", ...)} # Concatenate strings with a period as the separator

# Install 'devtools' and 'BiocManager' if not already installed
if (! "devtools" %in% installed_packs){
  install.packages("devtools")
}
if (! "BiocManager" %in% installed_packs){
  install.packages("BiocManager")
}

# Loop through each package and install if not already installed
for (pack in packs_to_check){
  if (!pack %in% installed_packs){
    tryCatch({install.packages(pack)},
              error = function(cond){message(cond);return(NA)},
              warning = function(cond){message(cond);return(NA)},
              finally = {message("install.packages tried")})
    
    if (!pack %in% installed_packs){
    tryCatch({BiocManager::install(pack)},
              error = function(cond2){message(cond2);return(NA)},
              warning = function(cond2){message(cond2);message(pack %+% " is not available.\n" %+% "Please search github and install by devtools::install_github");return(NA)},
              finally = {message("BiocManager::install tried")})
    }
  } else {
    message(pack %+% " already installed.")
  }
}


## Load Utilities ----
options(stringsAsFactors = F)


## Functions ----

#' @title paste_with_sep
#' @description Gluing elements together with the given delimiter.
#' 
#' @param eles Elements to link together
#' @param delimiter The delimiter to use (default is ",").
#' @return A string with elements linked by selected delimiter.
#' @export
paste_with_sep = function(eles, delimiter=","){
  return(paste0(eles,collapse=delimiter))
}



#' @title named_vec
#' @description Generating a named vector.
#' 
#' @param inames A vector of name list.
#' @param ivalues A vector of value list.
#' @return A named vector.
#' @export
named_vec = function (inames, ivalues){
  if (length(inames) != length(ivalues)) {
    warning("the lengths of the names and values are not equal! only returns the minimum length.")
    alen = min(length(inames), length(ivalues))
    inames = inames[1:alen]
    ivalues = ivalues[1:alen]
  }
  names(ivalues) = inames
  return(ivalues)
}



#' @title clear_scCellMediator_workspace()
#' @description identify and clear the work space of temporal files generated
#' 
#' @export
clear_scCellMediator_workspace = function(){
  rm(list = ls())
}



#' @title load_packages
#' @description Loading the packages required by scCellMediator
#' 
#' @param packs_to_check a vector of packages need to be checked
#' @export
load_packages = function(packs_to_check){
  if (all(packs_to_check %in% installed.packages() == TRUE)){
    loaded = check_loaded(packs_to_check)
    if (!is.null(loaded)){
      for (pack in loaded){
        library(pack, character.only = T, quietly = TRUE)
        message(pack %+% " loaded.")
      }
    }
  }
}



#' @title uniqlen
#' @param avec a vector
#' @return the number of unique elements in the vector
#' @export
uniqlen = function(avec){
  return(length(unique(avec)))
}



#' @title check_loaded
#' @description  Check whether required packages have been loaded. If not, load these packages.
#' 
#' @export
check_loaded = function(loaded_to_check = NULL){
  loaded_packages = (.packages())
  if (is.null(loaded_to_check)){
    loaded_to_check = packs_to_check
  }
  check_id = which(!loaded_to_check %in% loaded_packages)
  if (length(check_id) > 0){
    return(loaded_to_check[check_id])
    # print(paste_with_sep(loaded_to_check[check_id]) %+% " are not loaded. Please check, install and load them before run.")
  } else {
    return(NULL)
  }
}



#' @title get_random_colorlist
#' @description Get random generated color list
#' 
#' @param valList a list of value
#' @param named if return the name of input valList
#' @return the generated color list.
#' @export
get_random_colorlist = function(valList, named = T){
  check_loaded("grDevices")
  DISTINCT28 <- get_distinct_colorlist()
  colList <- (grDevices::colorRampPalette(sample(DISTINCT28,
                                                 size = min(uniqlen(DISTINCT28), uniqlen(valList)), replace = F)))(uniqlen(valList))
  if (named) {
    names(colList) <- unique(valList)
  }
  return(colList)
}



#' @title get_distinct_colorlist
#' @description Get 28 most distinct customed color values
#' 
#' @return a list of color list
get_distinct_colorlist = function(){
  DISTINCT28 = c("#e6194b","#3cb44b","#ffe119","#4363d8","#f58231","#911eb4","#46f0f0","#f032e6",
                 "#bcf60c","#fabebe","#008080","#e6beff","#9a6324","#fffac8","#800000","#aaffc3",
                 "#808000","#ffd8b1","#000075","#808080","#ff88ff","#000000","darkblue","darkorange",
                 "steelblue","cadetblue3","sienna1","lightyellow")
  return(DISTINCT28)
}




#' @title savePlot
#' @description save ggplot2 object
#'
#' @param pObj the ggplot2 object
#' @param pDir the path of directory to output.
#' @param pType types of plot, such as barplot, ggplot, circleplot, etc.
#' @param pName a recognizable name to the plot
#' @param pWd the width of the output plot. Default is 7.
#' @param pHt the height of the output plot. Default is 7.
#' @param pFmt the format of the figure. Default is .pdf.
#' @param dStamp a boolean value to decide whether add time stamp to the output filename.
#' @export
savePlot = function (pObj, pDir = ".", pType = "ggplot", pName = "", pWd = 7,
          pHt = 7, pFmt = "pdf", dStamp = T){
  if (dStamp) {
    ifilename = file.path(pDir, "p" %++% pType %++% pName %++%
                            getdstamp() %++% pFmt)
  }
  else {
    ifilename = file.path(pDir, "p" %++% pType %++% pName %++%
                            pFmt)
  }
  createDir(pDir)
  ggplot2::ggsave(filename = ifilename, plot = pObj, device = pFmt,
                  width = pWd, height = pHt)
}



#' @title createDir
#' @description If directory not exist, create the directory. Otherwise, ignore the command.
#' 
#' @param dirpath The directory path to check.
createDir = function (dirpath){
  if (!dir.exists(dirpath)) {
    dir.create(dirpath, recursive = T)
  }
}



#' @title get_functions_to_view
#' @description get the gene sets of the query functions with given names.
#' 
#' @param keywords keywords to inspect or full names to extract
#' @param mod two values accepted. There are two modes. 'inspect' to see available terms to use. 'catch' to obtain gene sets with given terms.
#' @param overlapKeywords a boolean value. Default is FALSE. If TRUE, only terms with all keywords occurred in the term name will be used.
#' @param species the species used for the analysis.
#' @param cat MSigDB category to search. If want to search all, set cat to NULL. Default is C5.
#' @return a list of given functions
#' @export
get_functions_to_view = function(keywords, mod = "inspect", overlapKeywords = FALSE, species = "mouse", cat = NULL){

  if (species == "mouse") {
    species = "Mus musculus"
  } else if (species == "human"){
    species = "Homo sapiens"
  }
  if (is.null(cat)){
    funcDB = msigdbr::msigdbr(species = species)
  } else {
    if (cat %in% c(paste0("C", 1:8), "H")){
      funcDB = msigdbr::msigdbr(species = species, category = cat)
    } else if (cat %in% c("CGN","CGP","CM","CP","CP:BIOCARTA","CP:KEGG","CP:PID","CP:REACTOME","CP:WIKIPATHWAYS","GO:BP","GO:CC","GO:MF","HPO","IMMUNESIGDB","MIR:MIR_Legacy","MIR:MIRDB","")){
      
    } else {
      
    }
  }

  keywords = toupper(keywords)

  funcToView = data.frame()
  for (kw in keywords){
    funcToView = rbind(funcToView, funcDB %>% filter(gs_name %in% (grep(kw, gs_name,value=T) %>% unique)))
  }
  funcToView = dplyr::distinct(funcToView)
  if (mod == "inspect"){
    funcToView = funcToView %>% pull(gs_name) %>% unique
    if (overlapKeywords){
      for (kw in keywords){
        funcToView = grep(kw,funcToView,value=T)
      }
    }
    if (length(funcToView) > 0){print(funcToView);return(funcToView)} else {print("No terms satisfied those keywords!")}
  } else if (mod == "catch"){
    funcToUse = intersect(funcDB$gs_name %>% unique, keywords)
    funcGeneSets = list()
    if (funcToUse > 0){
      for (func in funcToUse){
        funcGeneSets[[func]] = funcDB %>% filter(gs_name == func) %>% pull(gene_symbol) %>% unique
      }
      return(funcGeneSets)
    } else {
      print("No terms found in the database with given term names!")
    }
  }
}



#' @title getdstamp
#' @description get date timestamp
#' 
#' @param strlen Default length is 6. Show year, month, date without delimiter.
#' @return a numbered string
#' @export
getdstamp = function (strlen = 6) {
  if (strlen == 6) {
    dstamp = format(Sys.time(), "%Y%m%d") %>% 
      stringr::str_sub(start = 3)
  } else if (strlen == 8) {
    dstamp = format(Sys.time(), "%Y%m%d")
  } else {
    dstamp = Sys.time() %>% 
      stringr::str_replace_all("-", "_") %>% 
      stringr::str_replace_all(":", "_") %>% 
      stringr::str_replace_all(" ", "_")
  }
  return(dstamp)
}
