library(rtracklayer)
session <- browserSession()
genome(session)<-"hg38"
library(ChIPpeakAnno)
cpg <- session[["CpG Islands"]]
cpg_df <- as.data.frame(cpg)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg38_genes <- as.data.frame(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
annoData <- annoGR(TxDb.Hsapiens.UCSC.hg38.knownGene)
annoData_df<-as.data.frame(annoData)
cpg.anno <- annotatePeakInBatch(cpg, AnnotationData=annoData, output=c("upstream"),FeatureLocForDistance=c("TSS"),maxgap = 1000L,ignore.strand=FALSE)
cpg_anno_df<-as.data.frame(cpg.anno)
cpg_anno_df<-na.omit(cpg_anno_df)

library(AnnotationDbi)
library(org.Hs.eg.db)
ENTREZID_org <- keys(org.Hs.eg.db, keytype = "ENTREZID")
gene_dataframe_orgDb <- AnnotationDbi::select(org.Hs.eg.db, keys=ENTREZID_org, columns=c("ENSEMBL", "SYMBOL"), keytype="ENTREZID")
colnames(gene_dataframe_orgDb) <- c("Entrez", "Ensembl", "HGNC")
colnames(gene_dataframe_orgDb) <- paste(colnames(gene_dataframe_orgDb), "orgDb", sep = "_")

merged_gene_df<-merge(cpg_anno_df, gene_dataframe_orgDb, by.x="feature", by.y="Entrez_orgDb")
merged_gene_df<-merged_gene_df[merged_gene_df$shortestDistance<1500L,]
gene_names_with_CGI<-as.data.frame(unique(merged_gene_df$HGNC_orgDb))

data_path = "~/PycharmProjects/TCGA_Analyze/snv_dat/"
write.table(gene_names_with_CGI[ , 1], file= paste(data_path, "gene_names_with_CGI.txt",sep = ""),row.names = FALSE,col.names=FALSE,qmethod = "double", sep="\t")
