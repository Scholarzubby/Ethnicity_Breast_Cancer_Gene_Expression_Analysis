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

#didnt work for renaming my code
colnames(nonadata) <- c(gene, sample, fpkm, character, description, country, tissues)

#renaming row (the first function i used didnt work, so i used this)
colnames(nonadata)[6] <- "country"
colnames(nonadata)[7] <- "tissue"

nigsubdata <- filter(nonadata, country== "Nig")
aasubdata <- filter(nonadata, country== "AA")
causubdata <- filter(nonadata, country== "Cau")


#selecting nigerian basal-like breast cancer data
nigsubdatabasal <- filter(nigsubdata, tissue== "Basal")
aasubdatabasal <- filter(aasubdata, tissue== "Basal")
causubdatabasal <- filter(causubdata, tissue== "Basal")

#checking for lowly expressed genes
nolowreadnigs <- rowSums(count(nigsubdatabasal)) >= 10
head(nolowreadnigs)

nigsubdatabasalnot <- filter(nigsubdata, !tissue== "Basal")
nigsubdatabasalher2not <- filter(nigsubdatabasalnot, !tissue== "Her2")
nigsubdatabasalher2normalnot <- filter(nigsubdatabasalher2not, !tissue== "Normal")

#run deseQ #didint work
workingdata <- DESeq(nigsubdatabasalher2normalnot)






################################################

#colnames(nonadata)[6] <- "country"
#colnames(nonadata)[7] <- "tissue"

#name rownames for columnname all in countdata
rownames(nonadata2) <- nonadata2[,1]
row.names(nonadata2)


all(colnames(data) %in% row.names(nonadata2))
..
datatry <- read.table(file = "GSE229005_raw_data.txt", header = TRUE, sep = "", nrows = 343)
#find issues of the countdata and coldata not matching s because of NA
nonadata2<- !select(data, c(345,346,347,348,349,350))
nonadata2<- select(data, !c(345,346,347,348,349,350))

colnames(nonadata)<- nonadata2[,1]
rownames(nonadata) <- nonadata2[,2]
newmetadatarow <- newmetadata[ ,1:39]
#changing row name
rownames(newmetadatarow)<- newmetadatarow[,1]

#DIDNT WORK FOR CHANGING ROW NAME BECAUSE OR NUMERS WERE THE ROWNAMES
rownames(nonadata2)<- nonadata2[,1]
#WHILE THIS ONE WORKED FOR NUMBERS
rownames(nonadata2)<- make.names( nonadata2[,1], unique = TRUE) #nonadata2[,1]
#DELETEING UNWANTED ROW
nonadata2 <- nonadata2[2:344]

#colnames(newmetadatarow)<- newmetadatarow[,2]
#deleteting unwanted row
newmetadatarow2<- newmetadatarow[2:39]
all(colnames(nonadata2) %in% row.names(newmetadatarow2))
##nonadata2 <- nonadata2[,1:344]
###try <- nonadata2[1:344]

#selecting columns
SECTECTNONADATA2 <- select(nonadata2, c(6,7))
SECTECTNEWMETADATROW2 <- select(newmetadatarow2, c(35,37))

#renaming the rows
colnames(SECTECTNEWMETADATROW2)[1] <- "country"
colnames(SECTECTNEWMETADATROW2)[2] <- "tissue"

#filtering the tissue
SECTECTNEWMETADATROW2tissue <- filter(SECTECTNEWMETADATROW2, !tissue== "Basal")
SECTECTNEWMETADATROW2tissueher2 <- filter(SECTECTNEWMETADATROW2tissue,  !tissue== "Her2")
SECTECTNEWMETADATROW2tissueher2normal <- filter(SECTECTNEWMETADATROW2tissueher2, !tissue== "Normal")


#SECTECTNANODATA2tissue <- filter(nonadata2, !tissue== "Basal")
#SECTECTNANAODATA2tissueher2 <- filter(SECTECTNANODATA2tissue,  !tissue== "Her2")
#SECTECTNANAODATA2tissueher2normal <- filter(SECTECTNANAODATA2tissueher2, !tissue== "Normal")
TRY5 <- match(rownames(SECTECTNEWMETADATROW2tissueher2normal), colnames(nonadata2))



#factoring

maindeq <- DESeq(SECTECTNEWMETADATROW2)

#REMOVING UNMACTHING ROWS
MATCHING <- which()

#TRY2<- select_if(nonadata2,matches(row.names(SECTECTNEWM ETADATROW2tissueher2normal)))
TRY3 <- matches(SECTECTNEWMETADATROW2tissueher2normal, "nonadata2")
(all(colnames(nonadata2) = row.names(SECTECTNEWMETADATROW2tissueher2normal))
)
y <- data.frame( )

colnames(subsetnewmetadata)[6] <- "country"
colnames(subsetnewmetadata)[7] <- "tissue"


#factor leveling
subsetnewmetadatalevel<- factor(SECTECTNEWMETADATROW2$tissue, levels = c("Norma", "Basal", "Her2", "LumA", "LumB"))

subsetnewmetadatalevelfactor <- lapply(subsetnewmetadatalevel, as.factor)

#dds <- DESeqDataSetFromMatrix(countData = nonadata2, 
#colData = subsetnewmetadatalevelfactor,
#design = ~ tissue)

dim(nonadata2)

#make SUBSETNETMETADATALEVELFACTOR A DATDFRAME

subsetnewmetadatalevelfactordf<- as.data.frame(subsetnewmetadatalevelfactor)

dds <- DESeqDataSetFromMatrix(countData = nonadata2, 
                              colData = subsetnewmetadatalevelfactordf,
                              design = ~ tissue)



newmetadatarow2na<- na.omit(newmetadatarow2)
nonadata22 <- na.omit(nonadata2)
newmetadatarow2 
colnames(newmetadatarow2na)[37] <- "tissue"

dds <- DESeqDataSetFromMatrix(countData = nonadata22, 
                              colData = newmetadatarow2na,
                              design = ~ tissue)

nonadata22level<- factor(nonadata22$tissue, levels = c("Normal", "Basal", "Her2", "LumA", "LumB"))


#nonadata22levelfactor <- lapply(nonadata22level, as.factor)

dds <- DESeqDataSetFromMatrix(countData = nonadata22level, 
                              colData = newmetadatarow2na,
                              design = ~ tissue)



#checking the rownames amd colnames
rownames(newmetadatarow2) %in% colnames(nonadata22)



#####################################################################################
#make level with dataframe
#samepleinfo <- data.frame(newmetadatarow2$tissue, levels = c("Normal", "Basal", "Her2", "LumA", "LumB"))

# create level factor
newmetadatarow2na <- factor(newmetadatarow2na$tissue, levels = c("Normal", "Basal", "Her2", "LumA", "LumB"))

rownames(newmetadatarow2na) == colnames(nonadata22)
#selecting only contion
condition <- select(newmetadatarow2, c(37))
condition$tissue <- factor(condition$tissue, levels = c("Normal", "Basal", "Her2", "LumA", "LumB"))

rownames(condition) == colnames(nonadata22)

dds <- DESeqDataSetFromMatrix(nonadata22, condition, ~ tissue)


######################################################################################
IMPORTANT CODE USED FOR RUNNING MY ANALYSIS
##########################################################################################
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
newmetadatarow2<- newmetadatarow[2:39]

#wrangling countdata
#delecting colunms with na values in orignal data named data (then the number of columns matched the number of samples in the experiment)
#nonadata2<- !select(data, c(345,346,347,348,349,350)) #not sure if it the line under or this 
nonadata2<- select(data, !c(345,346,347,348,349,350))
#THIS ONE WORKED renaming rows with numbers which is refered to as unique values in R
rownames(nonadata2)<- make.names( nonadata2[,1], unique = TRUE)
#DELETEING UNWANTED column which has been turned to rownames
nonadata2 <- nonadata2[2:344]


condition <- select(newmetadatarow2, c(37))
condition$tissue <- factor(condition$tissue, levels = c("Normal", "Basal", "Her2", "LumA", "LumB"))

rownames(condition) == colnames(nonadata22)

dds <- DESeqDataSetFromMatrix(nonadata22, condition, ~ tissue)


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


#wrangling countdata
#delecting colunms with na values in orignal data named data (then the number of columns matched the number of samples in the experiment)
#nonadata2<- !select(data, c(345,346,347,348,349,350)) #not sure if it the line under or this 
nonadata2<- select(data, !c(345,346,347,348,349,350))
#THIS ONE WORKED renaming rows with numbers which is refered to as unique values in R
rownames(nonadata2)<- make.names( nonadata2[,1], unique = TRUE)
#DELETEING UNWANTED column which has been turned to rownames
nonadata2 <- nonadata2[2:344]
#REMOVE NA VALUES
nonadata22 <- na.omit(nonadata2)


rownames(condition) == colnames(nonadata22)

#factoring
condition <- select(newmetadatarow2, c(37))
condition$tissue <- factor(condition$tissue, levels = c("Normal", "Basal", "Her2", "LumA", "LumB"))


dds <- DESeqDataSetFromMatrix(nonadata22, condition, ~ tissue)
dds <- dds[rowSums(counts(dds)>2) >=4]
dds <-DESeq(dds)
resultsNames(dds)

res <- results(dds)

basalresults <- results(dds, name =  "tissue_Basal_vs_Normal" )
#resdf - this dataset in global environment side is the same as resdfbasal data set below
resdfbasal <- as.data.frame(basalresults)


summary(resdfbasal$log2FoldChange )

#save the resultsnames into a variable
result_name <- resultsNames(dds)
#to remove intercept from the results
resultname2 <- result_name[-1]
resultname2

#shrink fold change
#didnt work for shrinking the data
#res_shrunken <- lfcShrink(dds, coef = result_name, res = res)
#worked for shrinking the data
res_shrunken <- lfcShrink(dds, coef =2 )
#chane to data frame
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



