library(MIND)
library(Matrix)
library(data.table)
library(tidyverse)
library(biomaRt)

## bulk data
bulk = fread("../bulk_target_reID_FPKM_gene.tsv", header = T)
bulk_geneID = as.data.frame(bulk) %>% dplyr::select(gene_id) %>% separate(gene_id, c("Gene_ID","b"), '\\.') %>% dplyr::select(-b) %>% as.vector()

ensembl <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl', GRCh=37)
gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),filters='ensembl_gene_id', mart=ensembl, values = as.vector(bulk_geneID$Gene_ID))
colnames(gene_symbol) = c("Gene_ID","gene_name")

bulk = bulk %>% as_tibble() %>% dplyr::select(-gene_id) %>% separate(tracking_id, c("Gene_ID","b"), '\\.') %>% dplyr::select(-b) %>%
  left_join(gene_symbol, by = "Gene_ID") %>% dplyr::select(-Gene_ID) %>% drop_na() %>%
  arrange(gene_name) %>% dplyr::filter(!duplicated(gene_name)) %>%
  column_to_rownames(var="gene_name")
bulk = as.matrix(bulk)

## QC of bulk data

f_express = function(x){
sum(x > 0)
}
n = dim(bulk)[2]
n_express = apply(bulk, 1, f_express)
index = (n_express > n*0.8)

bulk = bulk[index,]

frac = read.table("../Est.prop.new_5cells_48inRefs_decon_416inBulk_reID.txt", row.names = 1, header = T)
frac = as.matrix(frac)

## reference
sc_ref <- readMM("filtered_count_matrix.mtx")
rownames(sc_ref) <- readLines("filtered_gene_row_names.txt")
sc_meta = fread("filtered_column_metadata.txt", header = T)
sc_meta = as.data.frame(sc_meta)

index = rownames(sc_ref) %in% rownames(bulk)
sc_ref = as.matrix(sc_ref[index,])

colnames(sc_ref) = make.names(as.vector(sc_meta$projid), unique = TRUE)
colnames(sc_meta) = c("TAG","sample","tsne1","tsne2","clusterIndex","cell_type","cells")

save.image(file = "snref_bulk.RData")


