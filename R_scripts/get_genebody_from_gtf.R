library(GenomicFeatures)
data_path = "~/PycharmProjects/TCGA_Analyze/data/GRCh38/"
txdb <- makeTxDbFromGFF(data_path + "Homo_sapiens.GRCh38.90.gtf", format="gtf")
genes <- genes(txdb)
gene_df <- as.data.frame(genes)
ensembl_ids<-gene_df[,6]
library(biomaRt)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
attributs_df<- as.data.frame(listAttributes(mart)) #list all available attributes

genes_table <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol", "description"), values= ensembl_ids, mart= mart)
merged_gene_df<-merge(gene_df, genes_table, by.x="gene_id", by.y="ensembl_gene_id")
write.table(merged_gene_df[ ,1:7], file= paste(data_path, "human_gene_bodys.tsv",sep = ""), col.names=F, sep="\t")