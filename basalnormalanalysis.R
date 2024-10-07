#filter all basal samples from the nigerian meta dataset
nigcomparebasalheatmap <- filter(nigcompare, tissue== "Basal")
#filter all normal samples from the nigerian meta dataset
nigcomparenormalheatmap <- filter(nigcompare, tissue == "Normal")
#combining the normal and the basal dataset
mergedbasalnormal <- combine(nigcomparebasalheatmap, nigcomparenormalheatmap)
#selecting all basal nigerian count data
nigcolcomparebaslheatmap <- select(nonadata22, c("BC004.13", "BC016.06", "BC018.13", "BC020.11_206301581116_01_08" , "BC020.12", "BC026.13", "BC031.13", "BC041.13_206422420217_01_11", "BC044.13.1",  "BC044.13.2" , "BC046.12_31020642251_01_02", "BC046.14"  , "BC047.06",  "BC047.10" , "BC047.14", "BC047.14_206422420217_01_01" , "BC051.08", "BC051.14", "BC051.14_206422550217_01_02" , "BC054.13" , "BC062.14" , "BC063.13"  ,   "BC066.08" , "BC066.14" , "BC075.13" , "BC080.13"  ,   "BC084.06", "BC084.13",    "BC085.08" , "BC088.13", "BC100.12_206422610217_01_04" , "BC122.13" ,  "BC123.13" , "BC130.13", "BC147.13" , "BC154.13", "BC161.13"  ,  "BC178.11_31020642253_01_05" , "H100.11",  "H1015.10", "H108.99",   "H1151.13", "H135.14"  ,  "H142.14" , "H1620.13", "H171.03", "H1873.07", "H2368.12_206422720217_01_09",  "H2374.10", "H2469.08", "H2599.03", "H2741.06", "H29.12_206422720217_01_08" , "H3077.08", "H694.13", "H76.13" , "H773.03", "H826.14", "H826.14_206422710217_H1626.14_06", "H840.13", "N010822", "N010828", "N010851",  "N010858", "N010862", "N010867",  "N010868", "TB170.FFPE", "TB170.FFPE.REPEAT",  "TB180", "TB245", "TB263.FFPE", "TB310.FFPE"
))
#combining both basal and normal counts data for nigeria
nigcolcomparebaslnormalheatmap <- select(nonadata22, c( "BC004.13", "BC016.06", "BC018.13", "BC020.11_206301581116_01_08" , "BC020.12", "BC026.13", "BC031.13", "BC041.13_206422420217_01_11", "BC044.13.1",  "BC044.13.2" , "BC046.12_31020642251_01_02", "BC046.14"  , "BC047.06",  "BC047.10" , "BC047.14", "BC047.14_206422420217_01_01" , "BC051.08", "BC051.14", "BC051.14_206422550217_01_02" , "BC054.13" , "BC062.14" , "BC063.13"  ,   "BC066.08" , "BC066.14" , "BC075.13" , "BC080.13"  ,   "BC084.06", "BC084.13",    "BC085.08" , "BC088.13", "BC100.12_206422610217_01_04" , "BC122.13" ,  "BC123.13" , "BC130.13", "BC147.13" , "BC154.13", "BC161.13"  ,  "BC178.11_31020642253_01_05" , "H100.11",  "H1015.10", "H108.99",   "H1151.13", "H135.14"  ,  "H142.14" , "H1620.13", "H171.03", "H1873.07", "H2368.12_206422720217_01_09",  "H2374.10", "H2469.08", "H2599.03", "H2741.06", "H29.12_206422720217_01_08" , "H3077.08", "H694.13", "H76.13" , "H773.03", "H826.14", "H826.14_206422710217_H1626.14_06", "H840.13", "N010822", "N010828", "N010851",  "N010858", "N010862", "N010867",  "N010868", "TB170.FFPE", "TB170.FFPE.REPEAT",  "TB180", "TB245", "TB263.FFPE", "TB310.FFPE", "BC059.14", "BC059.14..repeat.", "BC120.13", "H1153.13", "H1318.10", "H1551.12_31020642253_01_11", "H1836.09", "H500.12_206422550217_01_05"))


conditionnigbasal <- select(mergedbasalnormal, c(2))
conditionnigbasal$tissue <- factor(conditionnigbasal$tissue, levels = c("Normal", "Basal"))
ddsbasal <- DESeqDataSetFromMatrix(nigcolcomparebaslnormalheatmap, conditionnigbasal, ~ tissue)
ddsbasal <- ddsbasal[rowSums(counts(ddsbasal)>2) >=4]
ddsbasal <-DESeq(ddsbasal)
resultsNames(ddsbasal)
resnigbasalnormal <- results(ddsbasal)
#to show the number up and down regulated genes between basal and normal experimental samples
summary(resnigbasalnormal)
as.data.frame(resnigbasalnormal)

resnigbasalnormalnames <- resultsNames(ddsbasal)
#to remove intercept from the results
resnigbasalnormalnames2 <- resnigbasalnormalnames[-1]
res_shrunkenbasalnormal <- lfcShrink(ddsbasal, coef =2 )
res_shrunkendbasalnormal <- as.data.frame(res_shrunkenbasalnormal)

(resOrderedbasalnormal<- resnigbasalnormal[order(resnigbasalnormal$padj),]) #Check this out later by followimg the tutorial because it looks irrelevant

######################the code hashed below cause the rownames to be numbers instead of gene names as rownames
#nigbasal_vs_nignormalmain <- as.data.frame(resnigbasalnormal$log2FoldChange) #Check this out later by followimg the tutorial because it looks irrelevant
#while this code below with addition of rownames funtion used rownamaes of the resnigbasalnormal as rownames for the dataset
nigbasal_vs_nignormalmain <- as.data.frame(resnigbasalnormal$log2FoldChange, row.names = row.names(resnigbasalnormal) )

as.data.frame(resOrderedbasalnormal$log2FoldChange) #Check this out later

head(nigbasal_vs_nignormalmain)

plotMA(res_shrunkenbasalnormal, ylim = c(-2, 2))
#plot pca
ddspcabasalnormal = estimateSizeFactors(ddsbasal)
sebasalnormal <- SummarizedExperiment(log2(counts(ddspcabasalnormal, normalize = TRUE)+1), colData=colData (ddspcabasalnormal))
plotPCA(DESeqTransform(sebasalnormal), intgroup = "tissue")

# basal and normal heatmap
normcountbasalnormal <- counts(ddsbasal, normalized = T)
heatmap1basalnormal <- as.data.frame(normcountbasalnormal)
heatmap2basalnormal <- as.data.frame(resOrderedbasalnormal)
#add yes or no column
#creating a sig col to show if the gene expression values were significant or not
heatmap2basalnormal$sig <- ifelse(heatmap2basalnormal$padj <= 0.05, "yes", 'no')
#ma plot
ggplot(heatmap2basalnormal, aes(x = log10 (baseMean), y = log2FoldChange, colour = sig)) + geom_point()

#volcano plot
ggplot(heatmap2basalnormal, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) + geom_point()

#pheatmap
#subsetting gene expression values that were significate
signsbasalnormal <- subset(heatmap2basalnormal, padj <= 0.05)

#merge and filter
allsignsbasalnormal <- merge(heatmap1basalnormal, signsbasalnormal , by = 0)

#trial with the molhonk guy
#sa=tarting selection from 2 cos row name is 1, to 82 not 89 inorder to get rid of the unwanted colums e.g basemean, log2foldchange etc that are not counts
sigcountsbasalnormal <- allsignsbasalnormal[,2:82]
row.names(sigcountsbasalnormal) <- allsignsbasalnormal$Row.names
pheatmap::pheatmap(log2(sigcountsbasalnormal+1), scale = 'row' )



#for selecting only 8 sample from basal since normal has only 8 samples
sigcountsbasalnormal8 <- allsignsbasalnormal[,66:82]
row.names(sigcountsbasalnormal8) <- allsignsbasalnormal$Row.names
pheatmap::pheatmap(log2(sigcountsbasalnormal8+1), scale = 'row' )






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
sigcountsnignew <- allsignsnignew[,2:161]
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
