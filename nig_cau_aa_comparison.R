#select all nigeria data, cau data, aa data
#data name for all cau basal data ---caucomparebasalheatmap (and this is for coldata)
#for cau col data is below 
justcaubasalcount <- select(caucoldatamatch, c("PN0227.CC.T","PN0244.Y.T", "PN0247.G.T", "PN0258.IA.T", "PN0274.IK.T",  "PN0292.IFS.T", "PN0296.IAA.T", "X2003.130.1.brca1.mutated", "X2005.349.1.brca1.mutated", "X2005.89.1.brca1.mutated", "X2006.35.1", "X94.13.1.brca1.mut", "X95.42.1.brca1.mut", "X99.094.1.brca1.mut", "X99.144.1.brca1.mut", "X99.144.6.brca1.mutated", "X99.229.6",))

justcaubasalcount <- justcaubasalcount[1:110,]
#data name for all aa basal data --- aacomparebasalheatmap ( this for coldata)
justaabasalcount <- select(aacoldatamatch, c("PN0105.VVV.T", "PN0110.Y.T", "PN0120.H.T", "PN0221.H.T", "PN0223.C.T", "PN0226.N.T", "PN0231.E.T", "PN0248.1D.T", "PN0256.2N.T", "PN0261.II.T", "PN0263.IH.T", "PN0266.II.T", "PN0277.IF.T", "PN0281.IJ.T", "PN0286.4G.T", "X2001.221.1.brca1.mut", "X2001.228.1", "X2002.212.1", "X2003.166.1", "X2003.21.1.brca1.mutated", "X2003.85.1", "X2004.87.1.brca1.mut", "X2005.117.1", "X2005.184.1.brca1.methyl" , "X2005.238.1", "X2005.35.1", "X2005.54.1", "X2006.108.1.brca1.methyl", "X2006.148.1.brca1.methyl", "X2006.75.1", "X93.8676S", "X98.151.6"))
justaabasalcount <- justaabasalcount[1:110,]
#data name for all nig basal data ---nigcomparebasalheatmap
#below is only nig basal count data
justnigbasalcount <- select(nonadata22, c("BC004.13", "BC016.06", "BC018.13", "BC020.11_206301581116_01_08" , "BC020.12", "BC026.13", "BC031.13", "BC041.13_206422420217_01_11", "BC044.13.1",  "BC044.13.2" , "BC046.12_31020642251_01_02", "BC046.14"  , "BC047.06",  "BC047.10" , "BC047.14", "BC047.14_206422420217_01_01" , "BC051.08", "BC051.14", "BC051.14_206422550217_01_02" , "BC054.13" , "BC062.14" , "BC063.13"  ,   "BC066.08" , "BC066.14" , "BC075.13" , "BC080.13"  ,   "BC084.06", "BC084.13",    "BC085.08" , "BC088.13", "BC100.12_206422610217_01_04" , "BC122.13" ,  "BC123.13" , "BC130.13", "BC147.13" , "BC154.13", "BC161.13"  ,  "BC178.11_31020642253_01_05" , "H100.11",  "H1015.10", "H108.99",   "H1151.13", "H135.14"  ,  "H142.14" , "H1620.13", "H171.03", "H1873.07", "H2368.12_206422720217_01_09",  "H2374.10", "H2469.08", "H2599.03", "H2741.06", "H29.12_206422720217_01_08" , "H3077.08", "H694.13", "H76.13" , "H773.03", "H826.14", "H826.14_206422710217_H1626.14_06", "H840.13", "N010822", "N010828", "N010851",  "N010858", "N010862", "N010867",  "N010868", "TB170.FFPE", "TB170.FFPE.REPEAT",  "TB180", "TB245", "TB263.FFPE", "TB310.FFPE"))

justnigbasalcount <- justnigbasalcount[1:110,]

allcombinedbasalcount <- select(nonadata22, c("BC004.13", "BC016.06", "BC018.13", "BC020.11_206301581116_01_08" , "BC020.12", "BC026.13", "BC031.13", "BC041.13_206422420217_01_11", "BC044.13.1",  "BC044.13.2" , "BC046.12_31020642251_01_02", "BC046.14"  , "BC047.06",  "BC047.10" , "BC047.14", "BC047.14_206422420217_01_01" , "BC051.08", "BC051.14", "BC051.14_206422550217_01_02" , "BC054.13" , "BC062.14" , "BC063.13"  ,   "BC066.08" , "BC066.14" , "BC075.13" , "BC080.13"  ,   "BC084.06", "BC084.13",    "BC085.08" , "BC088.13", "BC100.12_206422610217_01_04" , "BC122.13" ,  "BC123.13" , "BC130.13", "BC147.13" , "BC154.13", "BC161.13"  ,  "BC178.11_31020642253_01_05" , "H100.11",  "H1015.10", "H108.99",   "H1151.13", "H135.14"  ,  "H142.14" , "H1620.13", "H171.03", "H1873.07", "H2368.12_206422720217_01_09",  "H2374.10", "H2469.08", "H2599.03", "H2741.06", "H29.12_206422720217_01_08" , "H3077.08", "H694.13", "H76.13" , "H773.03", "H826.14", "H826.14_206422710217_H1626.14_06", "H840.13", "N010822", "N010828", "N010851",  "N010858", "N010862", "N010867",  "N010868", "TB170.FFPE", "TB170.FFPE.REPEAT",  "TB180", "TB245", "TB263.FFPE", "TB310.FFPE", "PN0105.VVV.T", "PN0110.Y.T", "PN0120.H.T", "PN0221.H.T", "PN0223.C.T", "PN0226.N.T", "PN0231.E.T", "PN0248.1D.T", "PN0256.2N.T", "PN0261.II.T", "PN0263.IH.T", "PN0266.II.T", "PN0277.IF.T", "PN0281.IJ.T", "PN0286.4G.T", "X2001.221.1.brca1.mut", "X2001.228.1", "X2002.212.1", "X2003.166.1", "X2003.21.1.brca1.mutated", "X2003.85.1", "X2004.87.1.brca1.mut", "X2005.117.1", "X2005.184.1.brca1.methyl" , "X2005.238.1", "X2005.35.1", "X2005.54.1", "X2006.108.1.brca1.methyl", "X2006.148.1.brca1.methyl", "X2006.75.1", "X93.8676S", "X98.151.6","PN0227.CC.T","PN0244.Y.T", "PN0247.G.T", "PN0258.IA.T", "PN0274.IK.T",  "PN0292.IFS.T", "PN0296.IAA.T", "X2003.130.1.brca1.mutated", "X2005.349.1.brca1.mutated", "X2005.89.1.brca1.mutated", "X2006.35.1", "X94.13.1.brca1.mut", "X95.42.1.brca1.mut", "X99.094.1.brca1.mut", "X99.144.1.brca1.mut", "X99.144.6.brca1.mutated", "X99.229.6") )
allcombinedbasalcount <- allcombinedbasalcount[1:110,]
caucomparebasalheatmap
#selecting the basal from general metadata
basalallcombined <- filter(SECTECTNEWMETADATROW2, tissue== "Basal")
basalallcombined2 <- select (basalallcombined2, c(1))
basalallcombined3 <- select (basalallcombined2, c(1))
basalallcombined$tissue <- factor(basalallcombined$country, levels = c("Nig", "Cau", "AA"))
basalallcombined2$tissue <- factor(basalallcombined2$country, levels = c("Nig", "Cau", "AA"))
basalallcombined3$country <- factor(basalallcombined3$country, levels = c("Nig", "Cau", "AA"))

summary(basalallcombined3$country)
rownames(allbasalcondition) %in% colnames(allcombinedbasalcount)

ddsbasalall <- DESeqDataSetFromMatrix(allcombinedbasalcount, basalallcombined3, ~country)

ddsaabasalnormalcheckingsomething <- DESeqDataSetFromMatrix(aacomparebasalnormalheatmap, conditionAAbasalnormal, ~ tissue)

allbasalcondition <- vctrs::vec_c(nigcomparebasalheatmap, aacomparebasalheatmap, caucomparebasalheatmap)
allbasalcondition1<- select(allbasalcondition, c(1))
#rownames(allbasalcondition) <- rownames(allbasalcondition)

#allbasalcondition <- factor(allbasalcondition$country, levels = c("Nig", "Cau", "AA"))
#allbasalcondition2<- select(allbasalcondition, c(1))
allbasalcondition1$country <- factor(allbasalcondition1$country, levels = c("Nig", "AA", "Cau"))
#allbasalcondition1df <- as.data.frame(allbasalcondition1, row.names = row.names(allbasalcondition1df))
#rownames(allbasalcondition1df) <- rownames(allbasalcondition)
#summary(allbasalcondition2)
ddsbasalall <- DESeqDataSetFromMatrix(allcombinedbasalcount, allbasalcondition1, ~country)

ddsbasalall <- ddsbasalall[rowSums(counts(ddsbasalall)>2) >=4]
ddsbasalall <-DESeq(ddsbasalall)
resultsNames(ddsbasalall)
resbasalall <- results(ddsbasalall)
#--------------------------------------------------
#extracting results for "nig_vs_cau" 
basalallresultsnigcau <- results(ddsbasalall, name =  "country_Cau_vs_Nig"  )
resdfnigcau <- as.data.frame(basalallresultsnigcau)
summary(resdfnigcau$log2FoldChange )
#to show the number up and down regulated genes for only nig vs cau samples
summary(basalallresultsnigcau)
#extracting results for "nig_vs_aa" 
basalallresultsnigaa <- results(ddsbasalall, name =  "country_AA_vs_Nig"  )
resdfnigaa <- as.data.frame(basalallresultsnigaa)
summary(resdfnigaa$log2FoldChange )
#to show the number up and down regulated genes for only nig vs cau samples
summary(basalallresultsnigaa)

################nig_vs_aa#################################
#nigcomparebasal (countdata) & nigcomparebasalheatmap (coldata) | nigcomparebasal & caucomparebasal
nigAAbasaldata <- vctrs::vec_c(nigcomparebasalheatmap, aacomparebasalheatmap)
nigcaubasaldata <- vctrs::vec_c(nigcomparebasalheatmap, caucomparebasalheatmap)
#for countdata
nigcaubasalcoldata <- select(nonadata22, c("BC004.13", "BC016.06", "BC018.13", "BC020.11_206301581116_01_08" , "BC020.12", "BC026.13", "BC031.13", "BC041.13_206422420217_01_11", "BC044.13.1",  "BC044.13.2" , "BC046.12_31020642251_01_02", "BC046.14"  , "BC047.06",  "BC047.10" , "BC047.14", "BC047.14_206422420217_01_01" , "BC051.08", "BC051.14", "BC051.14_206422550217_01_02" , "BC054.13" , "BC062.14" , "BC063.13"  ,   "BC066.08" , "BC066.14" , "BC075.13" , "BC080.13"  ,   "BC084.06", "BC084.13",    "BC085.08" , "BC088.13", "BC100.12_206422610217_01_04" , "BC122.13" ,  "BC123.13" , "BC130.13", "BC147.13" , "BC154.13", "BC161.13"  ,  "BC178.11_31020642253_01_05" , "H100.11",  "H1015.10", "H108.99",   "H1151.13", "H135.14"  ,  "H142.14" , "H1620.13", "H171.03", "H1873.07", "H2368.12_206422720217_01_09",  "H2374.10", "H2469.08", "H2599.03", "H2741.06", "H29.12_206422720217_01_08" , "H3077.08", "H694.13", "H76.13" , "H773.03", "H826.14", "H826.14_206422710217_H1626.14_06", "H840.13", "N010822", "N010828", "N010851",  "N010858", "N010862", "N010867",  "N010868", "TB170.FFPE", "TB170.FFPE.REPEAT",  "TB180", "TB245", "TB263.FFPE", "TB310.FFPE","PN0227.CC.T","PN0244.Y.T", "PN0247.G.T", "PN0258.IA.T", "PN0274.IK.T",  "PN0292.IFS.T", "PN0296.IAA.T", "X2003.130.1.brca1.mutated", "X2005.349.1.brca1.mutated", "X2005.89.1.brca1.mutated", "X2006.35.1", "X94.13.1.brca1.mut", "X95.42.1.brca1.mut", "X99.094.1.brca1.mut", "X99.144.1.brca1.mut", "X99.144.6.brca1.mutated", "X99.229.6"))
nigaabasalcoldata <- select(nonadata22, c("BC004.13", "BC016.06", "BC018.13", "BC020.11_206301581116_01_08" , "BC020.12", "BC026.13", "BC031.13", "BC041.13_206422420217_01_11", "BC044.13.1",  "BC044.13.2" , "BC046.12_31020642251_01_02", "BC046.14"  , "BC047.06",  "BC047.10" , "BC047.14", "BC047.14_206422420217_01_01" , "BC051.08", "BC051.14", "BC051.14_206422550217_01_02" , "BC054.13" , "BC062.14" , "BC063.13"  ,   "BC066.08" , "BC066.14" , "BC075.13" , "BC080.13"  ,   "BC084.06", "BC084.13",    "BC085.08" , "BC088.13", "BC100.12_206422610217_01_04" , "BC122.13" ,  "BC123.13" , "BC130.13", "BC147.13" , "BC154.13", "BC161.13"  ,  "BC178.11_31020642253_01_05" , "H100.11",  "H1015.10", "H108.99",   "H1151.13", "H135.14"  ,  "H142.14" , "H1620.13", "H171.03", "H1873.07", "H2368.12_206422720217_01_09",  "H2374.10", "H2469.08", "H2599.03", "H2741.06", "H29.12_206422720217_01_08" , "H3077.08", "H694.13", "H76.13" , "H773.03", "H826.14", "H826.14_206422710217_H1626.14_06", "H840.13", "N010822", "N010828", "N010851",  "N010858", "N010862", "N010867",  "N010868", "TB170.FFPE", "TB170.FFPE.REPEAT",  "TB180", "TB245", "TB263.FFPE", "TB310.FFPE", "PN0105.VVV.T", "PN0110.Y.T", "PN0120.H.T", "PN0221.H.T", "PN0223.C.T", "PN0226.N.T", "PN0231.E.T", "PN0248.1D.T", "PN0256.2N.T", "PN0261.II.T", "PN0263.IH.T", "PN0266.II.T", "PN0277.IF.T", "PN0281.IJ.T", "PN0286.4G.T", "X2001.221.1.brca1.mut", "X2001.228.1", "X2002.212.1", "X2003.166.1", "X2003.21.1.brca1.mutated", "X2003.85.1", "X2004.87.1.brca1.mut", "X2005.117.1", "X2005.184.1.brca1.methyl" , "X2005.238.1", "X2005.35.1", "X2005.54.1", "X2006.108.1.brca1.methyl", "X2006.148.1.brca1.methyl", "X2006.75.1", "X93.8676S", "X98.151.6" ))
nigAAbasaldata$country <- factor(nigAAbasaldata$country, levels = c("Nig", "AA"))
nigcaubasaldata$country <- factor(nigcaubasaldata$country, levels = c("Nig", "Cau"))
ddsbasalnigvscau <- DESeqDataSetFromMatrix(nigcaubasalcoldata , nigcaubasaldata, ~country)
ddsbasalnigvsaa <- DESeqDataSetFromMatrix(nigaabasalcoldata , nigAAbasaldata, ~country)


ddsbasalnigvscau <- ddsbasalnigvscau[rowSums(counts(ddsbasalnigvscau)>2) >=4]
ddsbasalnigvsaa <- ddsbasalnigvsaa[rowSums(counts(ddsbasalnigvsaa)>2) >=4]

ddsbasalnigvscau <-DESeq(ddsbasalnigvscau)
ddsbasalnigvsaa <-DESeq(ddsbasalnigvsaa)

resultsNames(ddsbasalnigvscau)
resultsNames(ddsbasalnigvsaa)

resbasalnigvscau <- results(ddsbasalnigvscau)
resbasalnigvsaa <- results(ddsbasalnigvsaa)

resbasalnigvscaunames <- resultsNames(ddsbasalnigvscau)
resbasalnigvsaanames <- resultsNames(ddsbasalnigvsaa)


#to show the number up and down regulated genes between basal and normal experimental samples
summary(ddsbasalnigvscau)
summary(ddsbasalnigvsaa)
resbasalnigvscau2 <- resbasalnigvscau[-1]
resbasalnigvsaa2 <- resbasalnigvsaa[-1]

resbasalnigvscaunames2 <- resbasalnigvscaunames[-1]
resbasalnigvsaanames2 <- resbasalnigvsaanames[-1]


res_shrunkendbasalnigvscau<- lfcShrink(ddsbasalnigvscau, coef =2 )
res_shrunkendbasalnigvsaa<- lfcShrink(ddsbasalnigvsaa, coef =2 )


res_shrunkendbasalnigvscaudf <- as.data.frame(res_shrunkendbasalnigvscau)
res_shrunkendbasalnigvsaadf <- as.data.frame(res_shrunkendbasalnigvsaa)

(resOrderedbasalnigvscau <- resbasalnigvscau [order(resbasalnigvscau$padj),]) #Check this out later by followimg the tutorial because it looks irrelevant
(resOrderedbasalnigvsaa <- resbasalnigvsaa[order(resbasalnigvsaa$padj),]) #Check this out later by followimg the tutorial because it looks irrelevant

nigbasal_vs_caubasal <- as.data.frame(resbasalnigvscau$log2FoldChange, row.names = row.names(resbasalnigvscau) )
nigbasal_vs_AAbasal <- as.data.frame(resbasalnigvsaa$log2FoldChange, row.names = row.names(resbasalnigvsaa) )


plotMA(res_shrunkendbasalnigvscau, ylim = c(-2, 2))
plotMA(res_shrunkendbasalnigvsaa, ylim = c(-2, 2))




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
sigcountsAAbasalnormal <- allsignsAAbasalnormal[,26:41]
row.names(sigcountsAAbasalnormal) <- allsignsAAbasalnormal$Row.names
pheatmap::pheatmap(log2(sigcountsAAbasalnormal+1), scale = 'row' )











#----------------------------------------------------------
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
sigcountsAAbasalnormal <- allsignsAAbasalnormal[,26:41]
row.names(sigcountsAAbasalnormal) <- allsignsAAbasalnormal$Row.names
pheatmap::pheatmap(log2(sigcountsAAbasalnormal+1), scale = 'row' )


################comparing AA and NIG  Log2foldchange######################################################
nigaalog2foldchange <- merge(nigbasal_vs_nignormalmain, AAbasal_vs_AAnormalmain, by = 0)



