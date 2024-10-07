#View(aacoldatamatch)
#View(AAcompare)
#View(caucompare)
#View(caucoldatamatch)

#filter all basal samples from the nigerian meta dataset
aacomparebasalheatmap <- filter(AAcompare, tissue== "Basal")
#filter all normal samples from the nigerian meta dataset
aacomparenormalheatmap <- filter(AAcompare, tissue == "Normal")
#combining the normal and the basal dataset   ... N/B--- in the nig data set i used the function combine to combine both dat, but now i try using combine again and this what its telling me ---combine()` was deprecated in dplyr 1.0.0.,  Please use `vctrs::vec_c()` instead.
combinebasalnormal <- vctrs::vec_c(aacomparebasalheatmap, aacomparenormalheatmap)

aacomparebasalnormalheatmap <- select(aacoldatamatch, c("PN0105.VVV.T", "PN0110.Y.T", "PN0120.H.T", "PN0221.H.T", "PN0223.C.T", "PN0226.N.T", "PN0231.E.T", "PN0248.1D.T", "PN0256.2N.T", "PN0261.II.T", "PN0263.IH.T", "PN0266.II.T", "PN0277.IF.T", "PN0281.IJ.T", "PN0286.4G.T", "X2001.221.1.brca1.mut", "X2001.228.1", "X2002.212.1", "X2003.166.1", "X2003.21.1.brca1.mutated", "X2003.85.1", "X2004.87.1.brca1.mut", "X2005.117.1", "X2005.184.1.brca1.methyl" , "X2005.238.1", "X2005.35.1", "X2005.54.1", "X2006.108.1.brca1.methyl", "X2006.148.1.brca1.methyl", "X2006.75.1", "X93.8676S", "X98.151.6", "PN0172.M.T", "PN0179.J.T", "PN0179.J.T.repeat", "PN0220.CC.T", "PN0235.J.T", "PN0246.H.T", "PN0249.1S.T", "X2003.133.1"))
aacomparebasalnormalheatmap <- aacomparebasalnormalheatmap[1:110,]
conditionAAbasalnormal <- select(combinebasalnormal, c(2))
conditionAAbasalnormal$tissue <- factor(conditionAAbasalnormal$tissue, levels = c("Normal", "Basal"))

ddsaabasalnormal <- DESeqDataSetFromMatrix(aacomparebasalnormalheatmap, conditionAAbasalnormal, ~ tissue)
ddsaabasalnormal <- ddsaabasalnormal[rowSums(counts(ddsaabasalnormal)>2) >=4]
ddsaabasalnormal <-DESeq(ddsaabasalnormal)
resultsNames(ddsaabasalnormal)
resAAbasalnormal <- results(ddsaabasalnormal)
#to show the number up and down regulated genes between basal and normal experimental samples
summary(resAAbasalnormal)

resAAbasalnormalnames <- resultsNames(ddsaabasalnormal)
#to remove intercept from the results
resAAbasalnormalnames2 <- resAAbasalnormalnames[-1]
res_shrunkenAAbasalnormal <- lfcShrink(ddsaabasalnormal, coef =2 )
res_shrunkendAAbasalnormaldf <- as.data.frame(res_shrunkenAAbasalnormal)

(resOrderedAAbasalnormal<- resAAbasalnormal[order(resAAbasalnormal$padj),]) #Check this out later by followimg the tutorial because it looks irrelevant

AAbasal_vs_AAnormalmain <- as.data.frame(resAAbasalnormal$log2FoldChange, row.names = row.names(resAAbasalnormal) )


plotMA(res_shrunkenAAbasalnormal, ylim = c(-2, 2))
#plot pca
ddspcaAAbasalnormal = estimateSizeFactors(ddsaabasalnormal)
seAAbasalnormal <- SummarizedExperiment(log2(counts(ddspcaAAbasalnormal, normalize = TRUE)+1), colData=colData (ddspcaAAbasalnormal))
plotPCA(DESeqTransform(seAAbasalnormal), intgroup = "tissue")

# basal and normal heatmap
normcountAAbasalnormal <- counts(ddsaabasalnormal, normalized = T)
heatmap1AAbasalnormal <- as.data.frame(normcountAAbasalnormal)
heatmap2AAbasalnormal <- as.data.frame(resOrderedAAbasalnormal)
#add yes or no column
#creating a sig col to show if the gene expression values were significant or not
heatmap2AAbasalnormal$sig <- ifelse(heatmap2AAbasalnormal$padj <= 0.05, "yes", 'no')
#ma plot
ggplot(heatmap2AAbasalnormal, aes(x = log10 (baseMean), y = log2FoldChange, colour = sig)) + geom_point()

#volcano plot
ggplot(heatmap2AAbasalnormal, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) + geom_point()

#pheatmap
signsAAbasalnormal <- subset(heatmap2AAbasalnormal, padj <= 0.05)

#merge and filter
allsignsAAbasalnormal <- merge(heatmap1AAbasalnormal, signsAAbasalnormal , by = 0)

#trial with the molhonk guy
#sa=tarting selection from 2 cos row name is 1, to 82 not 89 inorder to get rid of the unwanted colums e.g basemean, log2foldchange etc that are not counts
sigcountsAAbasalnormal <- allsignsAAbasalnormal[,2:42]
row.names(sigcountsAAbasalnormal) <- allsignsAAbasalnormal$Row.names
pheatmap::pheatmap(log2(sigcountsAAbasalnormal+1), scale = 'row' )









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


