#biological context 
#code for installing annotation dbi

#code for installing org.Hs.eg.db 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gage")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gageData")

#covert gene symbol to gene id
gene <- c("ACTB", "ACTR3B" , "ADM",  "ANGPTL4", "ANLN", "AXL",  "BAG1", "BCL2", "BIRC5", "BLVRA", 
          "CAV1", "CCNB1","CCNE1","CD24", "CDC6", "CDC20","CDH3", "CENPF","CEP55","CLDN3",  
          "CLDN4","CLDN7","CRYAB","CXXC5","DDIT4","DDR1", "DSP",  "EGFR", "EMP3", "EPCAM",  
          "ERBB2","ERBB3","ERBB4","ESR1", "ESRP1","EVI2A","EXO1", "F11R", "FABP5","FBN1",
          "FGFR4","FLVCR2", "FOXA1", "FOXC1","GAL",  "GATA3","GNG11","GPR160",  "GRB7", "GRHL2",  
          "JUP",  "KIF2C","KRT5", "KRT8", "KRT14","KRT17","KRT19","LEPRE1",  "LHFP", "MAPT",
          "MDM2", "MELK", "MET",  "MIA",  "MKI67","MLPH", "MMP11","MPP1", "MRPL19",  "MYBL2",  
          "MYC",  "NAT1", "KNTC2","NDRG1","NT5E", "CDCA1","ORC6L","PGR",  "PHGDH","PIK3CA", 
          "PLOD1","PNP",  "PSMC4","PTEN", "PTTG1","PVRL3","RAB25","RAD17","RAD50","RB1", 
          "RPLP0","RRAGD","RRM2", "SF3A1","SFRP1","SH2B3","SLC16A3", "SLC39A6", "SPINT1",  "SPINT2", 
          "TMEM45B", "TMEM158", "TYMS", "UBE2C","UBE2T","UCHL1","VAMP8","VEGFA","VIM",  "ZEB1")

geneentrezid <- select(org.Hs.eg.db, keys =gene, 
       columns = c("ENTREZID"), keytype="SYMBOL")
geneentrezid2 <- select(org.Hs.eg.db, keys =gene, 
       columns = c("ENTREZID"), keytype="SYMBOL")
#58   LEPRE1      <NA>
 # 59     LHFP      <NA>
 # 73    KNTC2      <NA>
 # 76    CDCA1      <NA>
 # 77    ORC6L      <NA>
  #86    PVRL3      <NA>
#deleting rows 
resdf2 <- resdf
resdf2 <- resdf2[-c(58, 59,73,76,77,86),]


geneentrezid <- na.omit(geneentrezid)

rownames(geneentrezid)<- make.unique(geneentrezid[,2]) #nonadata2[,1]
#the code below didnt work because it was adding x to the numeric rowe name so make.unique prevented that
#rownames(resdf2)<- make.names( geneentrezid[,2], unique = TRUE) #nonadata2[,1]
rownames(resdf2)<- make.unique(geneentrezid[,2])
#merge dataframes by rownames
resdf2geneid <- (resdf2)
columns(org.Hs.eg.db)

foldchanges <- as.data.frame (resdf2$log2FoldChange, row.names = row.names(resdf2))
resdf2$symbol <- mapIds(org.Hs.eg.db, keys = row.names(resdf2),
                     column = "SYMBOL",
                     keytype = "ENTREZID",
                     multivalS = "first")
resdf2$entrez = mapIds(org.Hs.eg.db, keys = row.names(resdf2),
                    column = "ENSEMBL",
                    keytype = "ENTREZID",
                    multivalS = "first")
resdf2$name = mapIds(org.Hs.eg.db, keys = row.names(resdf2),
                  column = "GENENAME",
                  keytype = "ENTREZID",
                  multiVals = "first")
mapi

resdf2 <- resdf

install.packages("org.Hs.eg.db", repos="http://bioconductor.org/packages/3.2/bioc")

library(gage)
  library(gageData)
data("kegg.sets.hs")
data("sigmet.idx.hs")
keggres =gage(exprs = foldchanges, gsets = kegg.sets.hs, same.dir = TRUE)
data("go.subs.hs")
data("go.sets.hs")
library(tidyverse)
library(dplyr)
library(pathview)


gobpsets = go.sets.hs[go.subs.hs$BP]
gobpres =gage(foldchanges, gsets = gobpsets, same.dir = TRUE)

keggrespathways = data.frame(id = rownames(keggres$greater, keggress$greater)) %>%
                               tibble::as_tibble()%>%
                               filter(row_number() <= 20) %>%
                               .$id %>%
                               as.character()

keggrespathways
keggresids = substr(keggrespathways, start = 1, stop = 8)
keggresids
tmp = sapply(keggresids, function(pid) pathview(gene.data = foldchanges, pathway.id = pid, species = "hsa"))
#go.subs.hs$BP
#foldchanges2 <- resdf2$log2FoldChange
#names(foldchanges2)= resdf2$entrez
#head(foldchanges2)



#kegg.hs = kegg.gsets("human", id.type = "entrez")
#kegg.hs.sigmet = kegg.hs$kg.sets[kegg.hs$sigmet.idx]

#keggress = gage(foldchanges, gsets = kegg.hs.sigmet, same.dir = TRUE)
#lapply(keggress, head)

#greaters = keggress$greater
#lesser = keggress$less

#keggresspathways =data.frame(id = rownames(keggress$greater), keggress$greater)%>%
  #as.data.frame()%>%
  #rowSums(2) >=4 %>%
  #.$id%>%
  #as.character()
#keggresspathways
#row_num
#tb