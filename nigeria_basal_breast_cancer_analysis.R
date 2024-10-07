#select conditon
conditionnig <- select(nigcompare, c(2))
conditionnig$tissue <- factor(conditionnig$tissue, levels = c("Normal", "Basal", "Her2", "LumA", "LumB"))
####################THE NO PROBLEM STARTED  FORM THE MERGEDNIGDATA2 DATA
ddsnig <- DESeqDataSetFromMatrix(mergednigdata2, conditionnig, ~ tissue)
ddsnig <- ddsnig[rowSums(counts(ddsnig)>2) >=4]
ddsnig <-DESeq(ddsnig)
resultsNames(ddsnig)
resnigbasal <- results(ddsnig)
#to show the number up and down regulated genes for all experimental samples
summary(resnigbasal)


#extracting results for "tissue_Basal_vs_Normal" 
basalresultsnignew <- results(ddsnig, name =  "tissue_Basal_vs_Normal" )
#resdf - this dataset in global environment side is the same as resdfbasal data set below
resdfnigbasal <- as.data.frame(basalresultsnignew)
summary(resdfnigbasal$log2FoldChange )
#to show the number up and down regulated genes for only basal samples
summary(basalresultsnignew)

#extracting results for "tissue_Her2_vs_Normal" 
her2resultsnignew <- results(ddsnig, name =  "tissue_Her2_vs_Normal" )
resdfnigher2 <- as.data.frame(her2resultsnignew)
summary(resdfnigher2$log2FoldChange )
#to show the number up and down regulated genes for only her2 samples
summary(her2resultsnignew)

#extracting results for "tissue_LumA_vs_Normal" 
LumAresultsnignew <- results(ddsnig, name =  "tissue_LumA_vs_Normal"  )
resdfnigLumA <- as.data.frame(LumAresultsnignew)
summary(resdfnigLumA$log2FoldChange )
#to show the number up and down regulated genes for only LumA samples
summary(LumAresultsnignew)

#extracting results for "tissue_LumB_vs_Normal" 
LumBresultsnignew <- results(ddsnig, name =  "tissue_LumB_vs_Normal"  )
resdfnigLumB <- as.data.frame(LumBresultsnignew)
summary(resdfnigLumB$log2FoldChange )
#to show the number up and down regulated genes for only LumB samples
summary(LumBresultsnignew)


########BACK TO ANALYING ALL SAMPLES##############################
#save the resultsnames of all sample experiment into a variable
result_namenignew <- resultsNames(ddsnig)
#to remove intercept from the results
result_namenignew2 <- result_namenignew[-1]
res_shrunkennignew <- lfcShrink(ddsnig, coef =2 )
res_shrunkendfnignew <- as.data.frame(res_shrunkennignew)
#res_shrunken <- lfcShrink(dds, coef = result_name, res = res)
reslfcbasal <- lfcShrink(basalresultsnignew, coef = 2)

(resOrderednigbasal<- resnigbasal[order(resnigbasal$padj),]) #Check this out later by followimg the tutorial because it looks irrelevant

nigbasal_vs_nignormal <- as.data.frame(resnigbasal$log2FoldChange) #Check this out later by followimg the tutorial because it looks irrelevant
as.data.frame(resOrderednigbasal$log2FoldChange) #Check this out later

head(nigbasal_vs_nignormal)

plotMA(res_shrunkennignew, ylim = c(-2, 2))
#plot pca
ddspca = estimateSizeFactors(ddsnig)
se <- SummarizedExperiment(log2(counts(ddspca, normalize = TRUE)+1), colData=colData (ddspca))
plotPCA(DESeqTransform(se), intgroup = "tissue")

#############################################################################################
#normalize the dds data
############################################################################################
ddsnig <-DESeq(ddsnig)

#all sample condition heatmap
normcountnignew <- counts(ddsnig, normalized = T)
heatmap1nignew <- as.data.frame(normcountnignew)
heatmap2nignew <- as.data.frame(resOrderednigbasal) #i have seen were the one i was thinking was irrelevant so check well later
#add yes or no column
#creating a sig col to show if the gene expression values were significant or not
heatmap2nignew$sig <- ifelse(heatmap2nignew$padj <= 0.05, "yes", 'no')
#gg plot
ggplot(heatmap2nignew, aes(x = log10 (baseMean), y = log2FoldChange, colour = sig)) + geom_point()

#basal and normal heatmap
#all sample condition heatmap
normcountnignew <- counts(ddsnig, normalized = T)
heatmap1nignew <- as.data.frame(normcountnignew)
heatmap2nignew <- as.data.frame(resOrderednigbasal)
#add yes or no column
#creating a sig col to show if the gene expression values were significant or not
heatmap2nignew$sig <- ifelse(heatmap2nignew$padj <= 0.05, "yes", 'no')
#changing rownmaes from numbers to gene names
row.names(heatmap2nignew) <- row.names(mergednigdata2)
row.names(heatmap1nignew) <- row.names(mergednigdata2)
#ma plot
ggplot(heatmap2nignew, aes(x = log10 (baseMean), y = log2FoldChange, colour = sig)) + geom_point()



#volcano plot
ggplot(heatmap2nignew, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) + geom_point()

#pheatmap
#subsetting gene expression values that were significate
signsnignew <- subset(heatmap2nignew, padj <= 0.05)

#merge and filter
allsignsnignew <- merge(heatmap1nignew, signsnignew , by = 0)


#trial with the molhonk guy
sigcountsnignew <- allsignsnignew[,2:159]
#changing the number rownames into gene names
row.names(sigcountsnignew) <- allsignsnignew$Row.names
pheatmap::pheatmap(log2(sigcountsnignew+1), scale = 'row')


#trying heatmap with both no and yes
signsnignewtry <- subset(heatmap2nignew)

#merge and filter
allsignsnignewtry <- merge(heatmap1nignew, signsnignewtry , by = 0)


#trial with the molhonk guy
sigcountsnignewtry <- allsignsnignewtry[,2:182]
row.names(sigcountsnignewtry) <- allsignsnignewtry$Row.names
pheatmap::pheatmap(log2(sigcountsnignewtry+1), scale = 'row')











nigheatmap = select(nigcoldatamatch, c("BC004.13", "BC016.06", "BC018.13", "BC020.11_206301581116_01_08", "BC020.12", "BC026.13", "BC031.13", "BC041.13_206422420217_01_11", "BC044.13.1", "BC044.13.2", "BC059.14", "BC059.14..repeat.", "BC120.13", "H1153.13", "H1318.10", "H1551.12_31020642253_01_11", "H1836.09", "H500.12_206422550217_01_05"))

nigheatmap2 = nigheatmap[1:110,]
nigheatmap3 = nigheatmap[100:110,]
pheatmap::pheatmap(log2(nigheatmap3+1))
rownames(nigheatmap3)


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