#### module enrichment #####
#module load R/4.0.3-Python-3.8.5-Anaconda3-2020.11
# on the Kimel cluster
module load R/4.2.1
#install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#install.packages('./HDO.db_0.99.1.tar.gz', repos = NULL, type="source")
#BiocManager::install("clusterProfiler")
#update.packages(oldPkgs = "ggplot2")


install.packages('./rlang_1.0.6.tar.gz', repos = NULL, type="source")

library(rlang)
library(clusterProfiler)
library(DOSE)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(enrichplot)
library(gridExtra)
library(rms)
library(data.table)
library(ggraph)

​


genelist=read.csv("GO_forR.csv")
####

ensg <- genelist$ensg

Entrez=as.vector((mapIds(org.Hs.eg.db, keys=ensg, column="ENTREZID", keytype="ENSEMBL", multiVals="first")))
group=genelist$trait
group=group[!is.na(Entrez)]
Entrez=na.omit(Entrez)
mydf <- data.frame(Entrez, group)

#universe <- unique(as.vector(na.omit(mapIds(org.Hs.eg.db,
#                                            keys=names(net$colors), 
#                                            column="ENTREZID",
#                                            keytype="ENSEMBL",
#                                            multiVals="first"))))




ck1 <- compareCluster(Entrez~group, data=mydf, fun = "enrichGO", OrgDb='org.Hs.eg.db', ont="BP", minGSSize = 20, maxGSSize = 300)
summary(ck1)
#, universe=universe)
saveRDS(ck1, "compareClusterOutout_CT_GO_GWAS.rds")

ckr <- setReadable(ck1, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
cks <- simplify(ckr)

​
pdf("./CT_module_enrichment_dotplot.pdf",w=10/1.1,h=11/1.1)
clusterProfiler::dotplot(cks,font.size=8)
dev.off()

​
