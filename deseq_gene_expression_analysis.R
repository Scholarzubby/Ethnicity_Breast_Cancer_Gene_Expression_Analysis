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
#IMPORTANT CODE USED FOR RUNNING MY ANALYSIS
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

#preparing to extract the sample ids form col to row in order to match the data and select nig data for analysis
#first seecting the col in coldata that has the sample ids 
nigcolnames <- select(nigsubdata, c(2))
#turn the cols to row names
colnames(nigcolnames) <- nigcolnames[1:20020,]
rownamescol <- nigcolnames[2,]
rownames(nigcolnames) <- make.names( nigcolnames[,1], unique = TRUE)
row.names()
rownamescol$row.names(nigcolnames)

#transposing data to make first column rownames
nigcolnamestrans <- setNames(data.frame(t(nigcolnames[, -1])), nigcolnames[ , 1])
niganalysisdata <- inner_join(nonadata2, nigcolnamestrans)
commoncolnig <- generics::intersect(nonadata2, nigcolnamestrans)
intersectcol <- intersect(colnames(nonadata2), colnames(nigcolnamestrans))
#mergenigcol <- merge(nonadata2, nigcolnamestrans, by.x = nonadata2, by.y = nigcolnamestrans)
mergenigcol<- as.data.frame(intersect(names(nonadata2), names(nigcolnamestrans)))
mergenigcol2<- as.data.frame(intersect(colnames(nonadata2), colnames(nigcolnamestrans)))
mergenigcol3 <-  setNames(data.frame(t(mergenigcol2[, -1])), mergenigcol2[ , 1])
mergenigcol3trans <-  setNames(data.frame(t(mergenigcol2[, -1])), mergenigcol2[ , 1])
niganalysisdata <- inner_join(nonadata2, mergenigcol3trans)
intersectcol2 <- intersect(colnames(nonadata2), colnames(mergenigcol3trans))
###generics::intersect(nonadata2, newsubsetnig)

#selecting columuns with nigeria data
nigfirstset <- nonadata22[ ,1:160]
nigsecondset <- nonadata22[ ,267:285]
nigthirdtset <- nonadata22[ ,329:331]

#merging te selected data
mergednigdata <- merge(nigfirstset, nigsecondset)
mergednigdata <- merge(mergednigdata, nigthirdtset)
#checking if the subseted nig meta and counts data are the same
rownames(newsubsetnig) == colnames(mergednigdata)

###########################################################################################
#this analysis is to compare the gene expression of different breast caner subtypes in Nigerian breast cancer patients
###########################################################################################

conditionnig <- select(newsubsetnig, c(2))
conditionnig$tissue <- factor(conditionnig$tissue, levels = c("Normal", "Basal", "Her2", "LumA", "LumB"))

mergednigdata2 <- mergednigdata[1:110,]
ddsnig <- DESeqDataSetFromMatrix(mergednigdata2, conditionnig, ~ tissue)

ddsnig <- ddsnig[rowSums(counts(dds)>2) >=4]
ddsnig <-DESeq(ddsnig)
resultsNames(ddsnig)

resnig <- results(ddsnig)

basalresultsnig <- results(ddsnig, name =  "tissue_Basal_vs_Normal" )
#resdf - this dataset in global environment side is the same as resdfbasal data set below
resdfbasalnig <- as.data.frame(basalresultsnig)


summary(resdfbasalnig$log2FoldChange )

#save the resultsnames into a variable
result_namenig <- resultsNames(ddsnig)
#to remove intercept from the results
resultnamenig2 <- result_namenig[-1]
resultnamenig2
res_shrunkennig <- lfcShrink(ddsnig, coef =2 )
res_shrunkendfnig <- as.data.frame(res_shrunkennig)
#res_shrunken <- lfcShrink(dds, coef = result_name, res = res)


#reslfc <- lfcshrink(resdfbasal, coef = 2)

(resnigOrdered<- resnig[order(resnig$padj),])
summary(res_shrunkennig)
#just tried with res below and the result was different from that of res_shrunkedf
#summary(res)


normal_vs_basal <- as.data.frame(resnig$log2FoldChange)

head(normal_vs_basal )

plotMA(res_shrunkennig, ylim = c(-2, 2))




normal_vs_basal <- as.data.frame(resnig$log2FoldChange)

head(normal_vs_basal )

plotMA(res_shrunkennig, ylim = c(-2, 2), colNonsig = "")

ggplot(resnig, aes(x = log10(baseMean), y = log2FoldChange, colour = 1)) + geom_point()

plotMA(res_shrunken, ylim = c(-2, 2))

plotMA(res, ylim = c(-2, 2), colNonSig = "gray32", colSig = "red3",)

plotMA(res_shrunken, ylim = c(-2, 2), colNonSig = "gray32", colSig = "red3",)

plotMA(res_shrunkendf, ylim = c(-2, 2), colNonSig = "gray32", colSig = "red3", colLine = "ff000080", log = "x", cex = 0.45, xlab = "mean expression", ylab = "log fold change")



###########################################################################################
#getting biological context through pathway analysis
###########################################################################################
#installing packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")

#installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
library(AnnotationDbi)


#############################################################################################
#normalize the dds data
############################################################################################
ddsnig <-DESeq(ddsnig)

normcount <- counts(dds, normalized = T)
heatmap1 <- as.data.frame(normcount)
heatmap2 <- as.data.frame(resOrdered)
heatmap3 <- as.data.frame(resOrdered)
#add yes or no column
heatmap2$sig <- ifelse(heatmap2$padj <= 0.05, "yes", 'no')
#ma plot
ggplot(heatmap2, aes(x = log10 (baseMean), y = log2FoldChange, colour = sig)) + geom_point()


#volcano plot
ggplot(heatmap2, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) + geom_point()

#pheatmap
signs <- subset(heatmap2, padj <= 0.05)
signs2 <- subset(heatmap3, padj <= 0.05)

#merge and filter
allsigns <- merge(heatmap1, , by = 0)

#trying to select columuns that match the data of country of interest
#this diffculty is to be used the filtered from metadata (coldata) with specific contries to extract the
#info in countdata wih is based on ID and not country.
#nigmatch<- mergenigcol3trans[, names(mergenigcol3trans) %in% names(nonadata2)]
#nigmatchwithvalue <- rows_update(nigmatch, nonadata2) -- error from this code made me realised i swicthed X and Y
#The error message was All columns in `y` must exist in `x`. the problem was that i made the bigger dataframe Y,
#while the x was now the smaller data thats the subsetted dataframe
#I then swicthed x and y in the code below and it worked
nigmatch2<- rows_update(nonadata2, nigmatch)

#filtering data to get only the Cau  patient data from coldata/metadta
newsubsetcau <-subset(SECTECTNEWMETADATROW2, country == "Cau")

#filtering data to get only the aa  patient data from coldata/metadta
##newsubsetAA <-subset(SECTECTNEWMETADATROW2, country == "AA")
#first seecting the col in coldata that has the sample ids 
caucolnames <- select(causubdata, c(2))
aacolnames <- select(aasubdata, c(2))

#making rownames colnames using transpose
caucolnamestrans <- setNames(data.frame(t(caucolnames[, -1])), caucolnames[ , 1])
aacolnamestrans <- setNames(data.frame(t(aacolnames[, -1])), aacolnames[ , 1])

#select only target patients using the same method used in the nigeria dataset to filtering and mergering the columnnames and extracting the values
caumatch2<- !rows_update(caucolnamestrans, nonadata2)
aamatch2<- rows_update(nonadata2, aacolnamestrans)
#selecting the id"s of nig data
nigcoldatamatch <- select(nonadata2, c("BC001.14", "BC002.13",,"BC003.14", "BC004.13",,"BC004.14", "BC004.14_206301581116_01_01", "BC006.12_206301581116_01_02", "BC011.13","BC016.06", "BC017.13", "BC018.13", "BC020.11_206301581116_01_08", "BC020.12", "BC023.13", "BC026.13", "BC030.13", "BC031.13", "BC033.13", "BC033.13_206422420217_01_08", "BC034.14_206422400217_01_07", "BC035.06", "BC036.13", "BC036.13_206422420217_01_09", "BC037.14", "BC041.12_206422590217_01_03", "BC041.13_206422420217_01_11", "BC044.12_31020642251_01_01", "BC044.13.1", "BC044.13.2", "BC044.14_206422550217_01_01" , "BC046.07", "BC046.12_31020642251_01_02" ,"BC046.14", "BC047.06", "BC047.10",, "BC047.14", "BC047.14_206422420217_01_01", "BC048.14", "BC048.14_206422420217_01_12", "BC051.08", "BC051.14", "BC051.14_206422550217_01_02", "BC054.13", "BC055.14","BC055.14_206422550217_01_04", "BC056.14_206422400217_01_11", "BC059.14",, "BC059.14..repeat.", "BC060.13", "BC062.14",,"BC063.13", "BC066.08","BC066.14", "BC075.13", "BC080.13", "BC081.13_206422550217_01_11", "BC082.13", "BC084.06","BC084.13",, "BC085.08","BC085.12_31020642252_01_05" , "BC088.13", "BC092.13", "BC095.13", "BC097.13", "BC098.13", "BC098.13_206422550217_01_12", "BC100.12_206422610217_01_04", "BC103.13", "BC116.13" ,"BC120.13", "BC122.13", "BC123.13",, "BC130.13", "BC133.12_31020642252_01_08", "BC133.13_206422550217_01_07", "BC134.13", "BC135.13",,"BC147.13", "BC150.13", "BC153.05", "BC154.13", "BC160.08", "BC161.13", "BC164.13", "BC169.11_31020642253_01_02", "BC178.11_31020642253_01_05", "BC190.08", "BC192.11_31020642253_01_07", "BC219.07", "BC221.06", "H100.11", "H1015.10", "H1048.12_206422610217_01_10", "H1076.13", "H108.99", "H1151.13", "H1153.13", "H1212.12_31020642253_01_09", "H1309.03", "H1318.10", "H1340.07", "H1345.06", "H135.14", "H142.14", "H151.07", "H1543.13", "H1551.12_31020642253_01_11", "H1587.14", "H1587.14_206422710217_H1626.14_03", "H1620.13", "H171.03", "H1832.12_31020642253_01_10", "H1836.09", "H1873.07", "H2016.08", "H2118.08", "H2120.04", "H2163.08", "H2165.08", "H2368.12_206422720217_01_09", "H2374.10", "H2444.07", "H2469.08", "H2491.04", "H2513.09", "h2517.08", "H2599.03", "H2741.06", "H29.12_206422720217_01_08", "H3077.08", "H3285.10", "H3286.10", "H3371.05", "H3409.09", "H499.12_206422720217_01_10", "H500.12_206422550217_01_05", "H600.13", "H601.13", "H694.13", "H76.13", "H773.03", "H806.14", "H806.14_206422710217_H1626.14_02", "H809.14", "H826.14", "H826.14_206422710217_H1626.14_06", "H840.13", "H841.13", "N010820", "N010822", "N010828", "N010829", "N010842", "N010851", "N010858", "N010862", "N010867", "N010868", "N010877", "TB154.FFPE", "TB170.FFPE", "TB170.FFPE.REPEAT", "TB180", "TB189", "TB217.FFPE", "TB245", "TB251", "TB251.FFPE.REPEAT", "TB263.FFPE", "TB280.FFPE", "TB283", "TB310.FFPE", "TB313", "X10_BX_Ng", "X13_BX_NG", "X14_BX_NG", "X15_BX_NG", "X16_BX_NG", "X3_bx_ng", "X5_bx_ng", "X6_BX_NG"
))
#for AA
newsubsetaa <-subset(SECTECTNEWMETADATROW2, country == "AA")
caucolnamestrans <- setNames(data.frame(t(caucolnames[, -1])), caucolnames[ , 1])
aacolnamestrans2 <- setNames(data.frame(t(newsubsetaa[, -1])), newsubsetaa[ , 1])

#copying out the colnames of the trans of the country data set
colnames(aacolnamestrans)



aacolnames <- select(aasubdata, c(2))
#turn the cols to row names
rownames(aacolnames) <- make.names( aacolnames[,1], unique = TRUE)

#transposing data to make first column rownames
aacolnamestrans <- setNames(data.frame(t(aacolnames[, -1])), aacolnames[ , 1])
colnames(aacolnamestrans)
#selecting the ids of aa data
aacoldatamatch <- select(nonadata2, c("PN0102.K.T.",   "PN0103.L.T",    "PN0104.BB.T", "PN0105.VVV.T",  "PN0110.Y.T",    "PN0112.DDD.T", "PN0114.T.T",    "PN0120.H.T",    "PN0123.D.T",  "PN0124.B.T",    "PN0153.L.T",    "PN0171.ZZ.T",  "PN0172.M.T",    "PN0173.T.T",    "PN0177.L.T",   "PN0178.L.T",    "PN0179.J.T",    "PN0179.J.T.repeat" ,"PN0220.CC.T",   "PN0221.H.T",    "PN0223.C.T", "PN0224.B.T",    "PN0225.D.T",    "PN0226.N.T",   "PN0230.A.T",    "PN0231.E.T",    "PN0235.J.T",   "PN0240.E.T",    "PN0241.K.T",    "PN0246.H.T",   "PN0248.1D.T",   "PN0249.1S.T",   "PN0251.IA.T",  "PN0255.IC.T",   "PN0256.2N.T",   "PN0260.IL.T",  "PN0261.II.T",   "PN0262.ID.T",   "PN0263.IH.T",  "PN0266.II.T",   "PN0268.IQ.T",   "PN0270.IF.T",  "PN0273.II.T",   "PN0275.IB.T",   "PN0277.IF.T",  "PN0279.2C.T",   "PN0281.IJ.T",   "PN0283.IC.T",  "PN0284.IJ.T",   "PN0286.4G.T",   "PN0289.ID.T",  "PN0299.2E.T",   "PN0301.IR.T",   "PN0303.IP.T",  "PNO282.2D.T",   "X2001.062.1",   "X2001.188.1",  "X2001.221.1.brca1.mut" ,    "X2001.228.1",   "X2001.239.1",  "X2001.263.1",   "X2002.161.1.brca2.mut" ,    "X2002.212.1",  "X2002.217.1",   "X2002.97.41",   "X2003.133.1",  "X2003.166.1",   "X2003.168.1.brca1.mutated", "X2003.21.1.brca1.mutated", "X2003.85.1", "X2004.101.1",   "X2004.229.1",  "X2004.261.1",   "X2004.400.1",   "X2004.406.1",  "X2004.418.1",   "X2004.87.1.brca1.mut",  "X2005.117.1",  "X2005.184.1.brca1.methyl" , "X2005.238.1",   "X2005.259.1",  "X2005.35.1",  "X2005.54.1",    "X2006.108.1.brca1.methyl", "X2006.148.1.brca1.methyl", "X2006.53.1",  "X2006.75.1", "X93.8676S", "X98.151.6", "X98.19284" ))

#for cau data
newsubsetcau2 <-subset(SECTECTNEWMETADATROW2, country == "Cau")
#print out rownames and copy and then paste notepad and modify
rownames(newsubsetcau2)
#selecting the ids of aa data
caucoldatamatch <- select(nonadata2, c("PN0101.Q.T", "PN0107.H.T",  "PN0109.M.T", "PN0117.EE.T", "PN0118.CC.T", "PN0119.M.T",  "PN0121.AA.T", "PN0159.Y.T", "PN0160.HH.T", "PN0163.B.T", "PN0164.CC.T", "PN0165.K.T", "PN0166.J.T", "PN0167.T.T", "PN0168.R.T", "PN0169.L.T", "PN0174.P.T", "PN0175.B.T", "PN0176.P.T", "PN0180.E.T", "PN0222.L.T", "PN0227.CC.T", "PN0228.M.T", "PN0236.G.T", "PN0237.G.T", "PN0238.J.T", "PN0239.L.T", "PN0242.B.T", "PN0244.Y.T",  "PN0245.J.T", "PN0247.G.T", "PN0250.IT.T", "PN0252.IQ.T", "PN0253.IL.T", "PN0254.5H.T", "PN0258.IA.T", "PN0259.IB.T", "PN0264.IH.T", "PN0274.IK.T", "PN0276.IC.T", "PN0280.2A.T",  "PN0285.ID.T", "PN0288.IW.T", "PN0290.3I.T", "PN0291.2U.T", "PN0292.IFS.T", "PN0296.IAA.T", "PN0297.IS.T", "PN0302.IP.T", "PN0306.T.T", "PN0308.Q.T", "X2001.70.1.brca1.mut", "X2002.228.1", "X2002.59.1.brca2.mut",  "X2003.130.1.brca1.mutated", "X2004.59.1.brca2.mut", "X2004.59.1.brca2.mut.1", "X2004.59.1.brca2.mut.2", "X2004_04_1.brca2.mut", "X2005.349.1.brca1.mutated", "X2005.89.1.brca1.mutated", "X2006.35.1", "X94.13.1.brca1.mut", "X95.42.1.brca1.mut",  "X98.11.6.brca2.mut", "X99.009.01.brca1.methyl", "X99.021.42.brca1.mutated", "X99.094.1.brca1.mut", "X99.144.1.brca1.mut", "X99.144.6.brca1.mutated",  "X99.229.6",))
#above the joining has been done
