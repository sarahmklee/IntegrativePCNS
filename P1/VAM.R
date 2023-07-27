# ===============================================================
# Full Vireo - VAM
# SL
# 2/28/22
# ===============================================================

# Load packages
library(Seurat)
library(patchwork)
library(tidyverse)
library(msigdbr)
library(VAM)

# Load data
load("~/OneDrive - Dartmouth College/PCNS/snRNAseq/VireoData/Processed_Combined/Annotated_RPCAintegrated_FullDataset_SeuratObj.rdata")
data <- full.integrated_v2
rm(full.integrated_v2)

# Get gene indices
# Read in ensembl id and gene names from features.tsv file
feature.data = read.delim("~/OneDrive - Dartmouth College/PCNS/snRNAseq/PCNS_Pool13/filtered_feature_bc_matrix/features.tsv.gz", 
                          header = FALSE, stringsAsFactors = FALSE)

ensembl.ids = feature.data[,1]
gene.names = feature.data[,2]

# Match genes left after quality control
genes.after.qc = rownames(data@assays$RNA@counts)

indices.to.keep = unlist(sapply(genes.after.qc, function(x){which(gene.names == x[1])}))

ensembl.ids = ensembl.ids[indices.to.keep]
gene.names = gene.names[indices.to.keep]

# Load msigdbr of interest
all_gene_sets = msigdbr(species = "Homo sapiens")

# Cell Type signatures
# Subset
C8_pathways <- all_gene_sets[(all_gene_sets$gs_cat %in% "C8"), ]
C8_brain_pathways <- C8_pathways[grepl("BRAIN", C8_pathways$gs_name), ]

# Replace "_" with "-"
C8_brain_pathways$gs_name <- gsub("_", "-", C8_brain_pathways$gs_name)
head(C8_brain_pathways)

# Hall mark signatures
# Keep hallmark sets
hallmark <- all_gene_sets[(all_gene_sets$gs_cat %in% "H"), ]

# Replace "_" with "-"
hallmark$gs_name <- gsub("_", "-", hallmark$gs_name)
head(hallmark)

# Descartes fetal cerebellum
# Subset
C8_d_cerebellum_pathways <- C8_pathways[grepl("DESCARTES_FETAL_CEREBELLUM", C8_pathways$gs_name), ]

# Replace "_" with "-"
C8_d_cerebellum_pathways$gs_name <- gsub("_", "-", C8_d_cerebellum_pathways$gs_name)
head(C8_d_cerebellum_pathways)

# Descarates fetal cerebrum
# Subset
C8_d_cerebrum_pathways <- C8_pathways[grepl("DESCARTES_FETAL_CEREBRUM", C8_pathways$gs_name), ]

# Replace "_" with "-"
C8_d_cerebrum_pathways$gs_name <- gsub("_", "-", C8_d_cerebrum_pathways$gs_name)
head(C8_d_cerebrum_pathways)

# Descartes fetal main
# Subset
C8_mainfetal_pathways <- C8_pathways[grepl("DESCARTES_MAIN_FETAL", C8_pathways$gs_name), ]

# Replace "_" with "-"
C8_mainfetal_pathways$gs_name <- gsub("_", "-", C8_mainfetal_pathways$gs_name)
head(C8_mainfetal_pathways)

# Durante adult olfactory
# Subset
C8_dur_adult_olf_pathways <- C8_pathways[grepl("DURANTE_ADULT_OLFACTORY", C8_pathways$gs_name), ]

# Replace "_" with "-"
C8_dur_adult_olf_pathways$gs_name <- gsub("_", "-", C8_dur_adult_olf_pathways$gs_name)
head(C8_dur_adult_olf_pathways)

# Fan embryonic cortex
# Subset
C8_fan_emb_ctx_pathways <- C8_pathways[grepl("FAN_EMBRYONIC_CTX", C8_pathways$gs_name), ]

# Replace "_" with "-"
C8_fan_emb_ctx_pathways$gs_name <- gsub("_", "-", C8_fan_emb_ctx_pathways$gs_name)
head(C8_fan_emb_ctx_pathways)

# Subset
C8_zhong_PFC_pathways <- C8_pathways[grepl("ZHONG_PFC", C8_pathways$gs_name), ]

# Replace "_" with "-"
C8_zhong_PFC_pathways$gs_name <- gsub("_", "-", C8_zhong_PFC_pathways$gs_name)
head(C8_zhong_PFC_pathways)



# Get list of pathways
brain_sets <- unlist(unique(C8_brain_pathways$gs_name))

# Create gene collection 
gene_set_id_list = list()

# Create gene collection
for(i in 1:length(brain_sets)){
  # Get list of genes
  set_gene_ids = unlist(C8_brain_pathways$ensembl_gene[C8_brain_pathways$gs_name %in% brain_sets[i]])
  
  # Create gene collection 
  gene_set_id_list[[i]] = set_gene_ids
  names(gene_set_id_list)[i] = brain_sets[i]
}

# Create the list of gene indices required by vam
gene_set_collection <- createGeneSetCollection(gene.ids = ensembl.ids, gene.set.collection = gene_set_id_list)

# Change default assay
DefaultAssay(data) <- "RNA"

# Run VAM
data <- vamForSeurat(data, gene.set.collection = gene_set_collection, 
                     center = F, gamma = T, sample.cov = F, return.dist = T)

DefaultAssay(object = data) = "VAMcdf"


# C8 Durante cerebellum
# Make gene set collection
# Get list of pathways
d_cerebellum_sets <- unlist(unique(C8_d_cerebellum_pathways$gs_name))

# Create gene collection 
gene_set_id_list = list()

# Create gene collection
for(i in 1:length(d_cerebellum_sets)){
  # Get list of genes
  set_gene_ids = unlist(C8_d_cerebellum_pathways$ensembl_gene[C8_d_cerebellum_pathways$gs_name %in% d_cerebellum_sets[i]])
  
  # Create gene collection 
  gene_set_id_list[[i]] = set_gene_ids
  names(gene_set_id_list)[i] = d_cerebellum_sets[i]
}

# Create the list of gene indices required by vam
gene_set_collection <- createGeneSetCollection(gene.ids = ensembl.ids, gene.set.collection = gene_set_id_list)

# Change default assay
DefaultAssay(data) <- "RNA"

# Run VAM
data <- vamForSeurat(data, gene.set.collection = gene_set_collection, 
                     center = F, gamma = T, sample.cov = F, return.dist = T)

DefaultAssay(object = data) = "VAMcdf"


# C8 Durante cerebrum
# Make gene set collection
# Get list of pathways
d_cerebrum_sets <- unlist(unique(C8_d_cerebrum_pathways$gs_name))

# Create gene collection 
gene_set_id_list = list()

# Create gene collection
for(i in 1:length(d_cerebrum_sets)){
  # Get list of genes
  set_gene_ids = unlist(C8_d_cerebrum_pathways$ensembl_gene[C8_d_cerebrum_pathways$gs_name %in% d_cerebrum_sets[i]])
  
  # Create gene collection 
  gene_set_id_list[[i]] = set_gene_ids
  names(gene_set_id_list)[i] = d_cerebrum_sets[i]
}

# Create the list of gene indices required by vam
gene_set_collection <- createGeneSetCollection(gene.ids = ensembl.ids, gene.set.collection = gene_set_id_list)

# Change default assay
DefaultAssay(data) <- "RNA"
# Run VAM
data <- vamForSeurat(data, gene.set.collection = gene_set_collection, 
                     center = F, gamma = T, sample.cov = F, return.dist = T)
DefaultAssay(object = data) = "VAMcdf"


# C8 Durante adult olfactory
# Make gene set collection
# Get list of pathways
d_ad_olf_sets <- unlist(unique(C8_dur_adult_olf_pathways$gs_name))

# Create gene collection 
gene_set_id_list = list()

# Create gene collection
for(i in 1:length(d_ad_olf_sets)){
  # Get list of genes
  set_gene_ids = unlist(C8_dur_adult_olf_pathways$ensembl_gene[C8_dur_adult_olf_pathways$gs_name %in% d_ad_olf_sets[i]])
  
  # Create gene collection 
  gene_set_id_list[[i]] = set_gene_ids
  names(gene_set_id_list)[i] = d_ad_olf_sets[i]
}

# Create the list of gene indices required by vam
gene_set_collection <- createGeneSetCollection(gene.ids = ensembl.ids, gene.set.collection = gene_set_id_list)

# Change default assay
DefaultAssay(data) <- "RNA"
# Run VAM
data <- vamForSeurat(data, gene.set.collection = gene_set_collection, 
                     center = F, gamma = T, sample.cov = F, return.dist = T)

DefaultAssay(object = data) = "VAMcdf"


# C8 Fan Embryonic Ctx
# Make gene set collection
# Get list of pathways
fan_emb_ctx_sets <- unlist(unique(C8_fan_emb_ctx_pathways$gs_name))

# Create gene collection 
gene_set_id_list = list()

# Create gene collection
for(i in 1:length(fan_emb_ctx_sets)){
  # Get list of genes
  set_gene_ids = unlist(C8_fan_emb_ctx_pathways$ensembl_gene[C8_fan_emb_ctx_pathways$gs_name %in% fan_emb_ctx_sets[i]])
  
  # Create gene collection 
  gene_set_id_list[[i]] = set_gene_ids
  names(gene_set_id_list)[i] = fan_emb_ctx_sets[i]
}

# Create the list of gene indices required by vam
gene_set_collection <- createGeneSetCollection(gene.ids = ensembl.ids, gene.set.collection = gene_set_id_list)

# Change default assay
DefaultAssay(data) <- "RNA"

# Run VAM
data <- vamForSeurat(data, gene.set.collection = gene_set_collection, 
                     center = F, gamma = T, sample.cov = F, return.dist = T)

DefaultAssay(object = data) = "VAMcdf"


# C8 zhong pfc
# Make gene set collection
# Get list of pathways
zhong_pfc_sets <- unlist(unique(C8_zhong_PFC_pathways$gs_name))

# Create gene collection 
gene_set_id_list = list()

# Create gene collection
for(i in 1:length(zhong_pfc_sets)){
  # Get list of genes
  set_gene_ids = unlist(C8_zhong_PFC_pathways$ensembl_gene[C8_zhong_PFC_pathways$gs_name %in% zhong_pfc_sets[i]])
  
  # Create gene collection 
  gene_set_id_list[[i]] = set_gene_ids
  names(gene_set_id_list)[i] = zhong_pfc_sets[i]
}

# Create the list of gene indices required by vam
gene_set_collection <- createGeneSetCollection(gene.ids = ensembl.ids, gene.set.collection = gene_set_id_list)

# Change default assay
DefaultAssay(data) <- "RNA"

# Run VAM
data <- vamForSeurat(data, gene.set.collection = gene_set_collection, 
                     center = F, gamma = T, sample.cov = F, return.dist = T)

DefaultAssay(object = data) = "VAMcdf"


# C8 main fetal
# Make gene set collection
# Get list of pathways
mainfetal_sets <- unlist(unique(C8_mainfetal_pathways$gs_name))

# Create gene collection 
gene_set_id_list = list()

# Create gene collection
for(i in 1:length(mainfetal_sets)){
  # Get list of genes
  set_gene_ids = unlist(C8_mainfetal_pathways$ensembl_gene[C8_mainfetal_pathways$gs_name %in% mainfetal_sets[i]])
  
  # Create gene collection 
  gene_set_id_list[[i]] = set_gene_ids
  names(gene_set_id_list)[i] = mainfetal_sets[i]
}

# Create the list of gene indices required by vam
gene_set_collection <- createGeneSetCollection(gene.ids = ensembl.ids, gene.set.collection = gene_set_id_list)

# Change default assay
DefaultAssay(data) <- "RNA"

# Run VAM
data <- vamForSeurat(data, gene.set.collection = gene_set_collection, 
                     center = F, gamma = T, sample.cov = F, return.dist = T)

DefaultAssay(object = data) = "VAMcdf"
