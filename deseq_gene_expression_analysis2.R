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

library("dplyr")
library("tidyverse")
library("GEOquery")
library(deseq)
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

nolowreadnigs <- rowSums(counts(nigsubdatabasal)) >= 10

colnames(nigsubdata)