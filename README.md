# sccm
`sccm` is a R computational framework that predicts triple-cell relationship regarding the biological functions of a certain cell type, developed and mained by [Jianing Zhang](https://www.undetermined.xx)

# Installation
To install the developmental version from Github:

```R
if (!require(remotes)) install.packages("remotes");
if (!require(R.filesets)) install.packages("R.filesets");
if (!require(stringr)) install.packages("stringr");
if (!require(dplyr)) install.packages("dplyr");
if (!require(purrr)) install.packages("purrr");
if (!require(tibble)) install.packages("tibble");
if (!require(Seurat)) remotes::install_github("satijalab/seurat", build_vignettes = TRUE);

remotes::install_github("zhjning/sccm", build_vignettes = TRUE)
```

`sccm` can accept the estimated ligand-receptor signal matrices from CellChat or iTALK. To prepare the CellChat and iTALK matrices displayed in the tutorial, CellChat or iTALK requires to be installed.

```R
# optional, recommend to install at least one of the cell-cell communication analysis tool.
if (!require(CellChat)) remotes::install_github("sqjin/CellChat");
if (!require(iTALK)) remotes::install_github("Coolgenome/iTALK");
```




