% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_cellchat.R
\name{run_cellchat}
\alias{run_cellchat}
\title{run_cellchat}
\usage{
run_cellchat(
  so,
  feature.name,
  batch.name = "orig.ident",
  output = NULL,
  minCNUM.tot = 10,
  sampleProp = 0.5,
  minCNUM.new = 10,
  runs = 10,
  suggestedRunTimes = F,
  species = "mouse",
  storeResampledSO = FALSE,
  customedSymbolList = NULL
)
}
\arguments{
\item{so}{seurat object}

\item{feature.name}{the cluster name used for resampling. In our examples, ct.major was used.}

\item{batch.name}{the batch name used for resampling. The variable used should be included in the meta.table from the analyzed seurat object. Default is "orig.ident".}

\item{output}{the path of the directory for cellchat results. Default is NULL.}

\item{minCNUM.tot}{the minCNUM.tot used for the bootstrap of seurat object. Default is 10.}

\item{minCNUM.new}{the minCNUM.new used for the bootstrap of seurat object. Default is 10.}

\item{runs}{the runs used for the bootstrap of seurat object. Default is 10.}

\item{suggestedRunTimes}{the suggestedRunTimes used for the bootstrap of seurat object. Default is FALSE.}

\item{species}{the species used for reference CellChat database offered by CellChat.}

\item{storeResampledSO}{a boolean variable. If TRUE, not save the resampled seurat sub-objects used for cell chat analysis. Default is FALSE.}

\item{customedSymbolList}{default is NULL. Pass a named vector for single cell studies with customed defined gene names. The names of the vector are customed defined names, the corresponding values are official gene symbols.}

\item{sampelProp}{the proportion value used for the bootstrap of seurat object. Default is 0.5.}
}
\value{
the path of the directory for cellchat results.
}
\description{
run cellchat
}
