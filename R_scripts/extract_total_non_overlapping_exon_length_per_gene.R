# First, import the GTF-file that you have also used as input for htseq-count
library(GenomicFeatures)
gtf_path = file.path("~/PycharmProjects/TCGA_Analyze/data/GRCh38/","gencode.v22.annotation.gtf")
txdb <- makeTxDbFromGFF(gtf_path, format="gtf")
# then collect the exons per gene id
exons_list_per_gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic_gene_sizes <- lapply(exons_list_per_gene,function(x){sum(width(reduce(x)))})
exonic_gene_sizes_df <- as.data.frame(exonic_gene_sizes)
