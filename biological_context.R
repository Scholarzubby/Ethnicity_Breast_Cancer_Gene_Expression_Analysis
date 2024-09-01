#biological context 
#code for installing annotation dbi

#code for installing org.Hs.eg.db 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

columns(org.Hs.eg.db)

foldchanges <- as.data.frame(res$log2folchange, row.names = row.names(res))
res$symbol <- mapids(org.Hs.eg.db, keys = row.names(res),
                     column = "SYMBOL"
                     keytypes = "ENSEMBL"
                     multival = "first")
res$entrez = mapids(org.Hs.eg.db, keys = row.names(res),
                    column = "ENTREZID"
                    keytypes = "ENSEMBL"
                    multival = "first")
res$name = mapids(org.Hs.eg.db, keys = row.names(res),
                  column = "GENENAME"
                  keytypes = "ENSEMBL"
                  multival = "first")



install.packages("org.Hs.eg.db", repos="http://bioconductor.org/packages/3.2/bioc")