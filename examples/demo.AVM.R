
############################################
## known reference markers for cell types ##
############################################

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


##############################
## Load lung cancer dataset ##
##############################

dataDir = "data/E-MTAB-6149/LungCancer/LC_counts/"
dataType = "10x"
metafile = "data/E-MTAB-6149/LungCancer/2097-Lungcancer_metadata.csv"

treatedDir = "data/E-MTAB-6149/LungCancer/LC_treated/"
createDir(treatedDir)
