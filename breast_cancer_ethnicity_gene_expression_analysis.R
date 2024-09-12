# Script to analyse gen expression data

#setwd("C:/Users/schol/OneDrive/Desktop/firsttrial2")

#load libraies
#install libraries 
#(install tidyverse from tools install. on intalling tidyvere dplyr comes with it)
#(to install GEOquery use this code: 
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("GEOquery")) at some point during the installation you will be require you to make an input "update all/some/none? [a/s/n]:" just input small letter a and press enter. then you'll get another prompt to make an input after "gert 2.0,1 2.1.0 True" at this point check the bottom icon bar on your system and pick icon that shows your download progress. you see a message saying something about compilation just select yes, and your download start again and then will be completed.
#if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("DESeq2")

#
install.packages("BiocManager")

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("apeglm")


library("dplyr")
library("tidyverse")
library("GEOquery")
ins
#DES
library("GenomeInfoDbData")
library("BiocManager")
library("apeglm")

################META Data#########################################################################
# read in data (my data is the same working directory)
#read.delim( file = "\Users\schol\Downloads\GSE229005_raw_data (1).txt", header = TRUE, sep = "\t")
#txt.data = read.delim(file.choose())
data <- read.delim("GSE229005_raw_data.txt")
#filter NA in data
newdata <- data[1:110,]
#get metadata that coldata GEO
metadata <- getGEO(GEO = "GSE229005", GSEMatrix = TRUE)
#data2 <- read.delim("GSE229005_series_matrix.txt")

#GSE229005_series_matrix.txt
newmetadata <- pData(phenoData(metadata[[1]]))
#selecting desired data
subsetnewmetadata <- select(newmetadata, c(1,10,12,22,36,38))
modimetadata <- newmetadata %>%
  select(1,10,22,36,38) %>%
  rename(character = characteristics_ch1)  %>%
  mutate(character = gsub("tissue: ", "", character))

###############COUNTDATA###############################################################################
#rename data column and gather to make into along format
datalong <- data %>%
  rename(gene = ID_REF)%>%
  gather(key = "samples", value = "FPKM",, -gene)

#join dataframes = datalong + modimetadata
datalong2 <- datalong %>%
  left_join(., modimetadata, by = c("samples" = "title"))
#checking foe na vaules
is.na(datalong2)
sum(is.na(datalong2))


#delete na value
nonadata<- na.omit(datalong2)

#renaming ROW
colnames(nonadata)[6] <- "country"
colnames(nonadata)[7] <- "tissue"

nonadata2<- select(data, !c(345,346,347,348,349,350))

rownames(nonadata2)<- make.names( nonadata2[,1], unique = TRUE) #nonadata2[,1]
#DELETEING UNWANTED ROW
nonadata2 <- nonadata2[2:344]


######################################################################################
# REARRANGED IMPORTANT CODE USED RUNNING MY ANALYSIS##################################
######################################################################################
data <- read.delim("GSE229005_raw_data.txt")
#filter NA in data
newdata <- data[1:110,]
#get metadata that coldata GEO
metadata <- getGEO(GEO = "GSE229005", GSEMatrix = TRUE)
#data2 <- read.delim("GSE229005_series_matrix.txt")

#wranging matrix data
#GSE229005_series_matrix.txt
newmetadata <- pData(phenoData(metadata[[1]]))
newmetadatarow <- newmetadata[ ,1:39] #this is just for duplication of the orginal newmetadat
#changing row names
rownames(newmetadatarow)<- newmetadatarow[,1] #making the samples the role name
#deleteting unwanted row
newmetadatarow2<- newmetadatarow[2:39] # NOT NECCESSARY
condition <- select(newmetadatarow2, c(37))
condition$tissue <- factor(condition$tissue, levels = c("Normal", "Basal", "Her2", "LumA", "LumB"))


########################################################################################################
#WRANGLING COUNTDATA###################################################################################
###############COUNTDATA###############################################################################
#renaming data column and gather to make into along format
datalong <- data %>%
  rename(gene = ID_REF)%>%
  gather(key = "samples", value = "FPKM",, -gene)

#join dataframes = datalong + modimetadata
datalong2 <- datalong %>%
  left_join(., modimetadata, by = c("samples" = "title"))

#checking foe na vaules
is.na(datalong2)
sum(is.na(datalong2))

#delete na value
nonadata<- na.omit(datalong2)

#renaming ROW
colnames(nonadata)[6] <- "country"
colnames(nonadata)[7] <- "tissue"

rownames(nonadata2)<- make.names( nonadata2[,1], unique = TRUE) #nonadata2[,1]

#DELETEING UNWANTED ROW
nonadata2 <- nonadata2[2:344]

#delecting colunms with na values in orignal data named data (then the number of columns matched the number of samples in the experiment)
#nonadata2<- !select(data, c(345,346,347,348,349,350)) #not sure if it the line under or this 
nonadata2<- select(data, !c(345,346,347,348,349,350))

#THIS ONE WORKED renaming rows with numbers as rownames which is refered to as unique values
rownames(nonadata2)<- make.names( nonadata2[,1], unique = TRUE)
#DELETEING UNWANTED column which has been turned to rownames
nonadata2 <- nonadata2[2:344]
#REMOVE NA VALUES
nonadata22 <- na.omit(nonadata2)
#checking if countdata and coldata matched
rownames(condition) == colnames(nonadata22)

#factoring the experiment
condition <- select(newmetadatarow2, c(37))
condition$tissue <- factor(condition$tissue, levels = c("Normal", "Basal", "Her2", "LumA", "LumB"))

#################begining of gene expression analysis using DESeq####################################
dds <- DESeqDataSetFromMatrix(nonadata22, condition, ~ tissue)
dds <- dds[rowSums(counts(dds)>2) >=4]
dds <-DESeq(dds)
resultsNames(dds)

res <- results(dds)
#checking the basal breast cancer subtype results###############################################
basalresults <- results(dds, name =  "tissue_Basal_vs_Normal" )
#resdf - this dataset in global environment side is the same as resdfbasal data set below
resdfbasal <- as.data.frame(basalresults)


summary(resdfbasal$log2FoldChange )

#save the resultsnames into a variable
result_name <- resultsNames(dds)
#to remove intercept from the results
resultname2 <- result_name[-1]
resultname2
res_shrunken <- lfcShrink(dds, coef =2 )
res_shrunkendf <- as.data.frame(res_shrunken)
#res_shrunken <- lfcShrink(dds, coef = result_name, res = res)


#reslfc <- lfcshrink(resdfbasal, coef = 2)

(resOrdered<- res[order(res$padj),])
summary(res_shrunken)
#just tried with res below and the result was different from that of res_shrunkedf
#summary(res)


flt_vs_fc <- as.data.frame(res$log2FoldChange)

head(flt_vs_fc)

plotMA(res_shrunken, ylim = c(-2, 2))

newsubset <-subset(SECTECTNEWMETADATROW2, country == "Nig")

newsubsetnig <-subset(SECTECTNEWMETADATROW2, country == "Nig")

rownames(newsubsetnig) & colnames(nonadata2)



##########################################################################################
#pathway analysis
#########################################################################################

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
