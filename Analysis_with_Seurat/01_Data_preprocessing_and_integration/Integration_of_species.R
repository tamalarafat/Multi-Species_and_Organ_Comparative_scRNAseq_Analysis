# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###
# Load data tables
###

###
# Protoplasting induced genes
###

# Load the table containing the list of protoplasting-induced genes.
PP_genes_table = read.csv("/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/Protoplasting_genes/Ox_Co0_leaf_protoplast_v12_final_August_2022_ortho.csv")

# Gene IDs - protoplasting-induced genes
PP_genes = PP_genes_table$GeneID

###
# Orthologues table
###

ortho_table = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Orthologues_n_correspondence/Orthos_table.csv")

###
# WT C. hirsuta
###

# Load data - WT OX 1st Experiment (leaf 5 and 6)
OX_data_1E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_OX_RNA_1ST_2_Newest/filtered_feature_bc_matrix/")

# Load data - WT OX 2nd Experiment (leaf 6 and 7)
OX_data_2E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_2nd_ALL_2_Newest/filtered_feature_bc_matrix/")

# Load data - WT OX 3rd Experiment (leaf 5 and 6)
OX_data_3E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_3rd_ALL_3000_Newest/filtered_feature_bc_matrix/")

# Load data - WT OX 7th Experiment (leaf 6 and 7)
OX_data_7E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_7th_ALL_2_Newest/filtered_feature_bc_matrix/")

# Convert the gene ids in the data table to ortho gene ids
OX_DF_1 = prepare_ortho_data(input_data = OX_data_1E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

OX_DF_2 = prepare_ortho_data(input_data = OX_data_2E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

OX_DF_3 = prepare_ortho_data(input_data = OX_data_3E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

OX_DF_7 = prepare_ortho_data(input_data = OX_data_7E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

###
# A. thaliana - SAM (Apex)
###

# WT COL data 1st Experiment - leaf 5 and 6
COL_SAM <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col0_RNA_12th_SAM_New/filtered_feature_bc_matrix/", gene.column = 1)


# All gene IDs - Arabidopsis Thaliana
thaliana_genes = rownames(COL_SAM)

# extracting the Cardamine IDs that are present in orthologues table 
thaliana_ortho_genes = as.character(ortho_table$A.thaliana.TAIR10)

# not all the thaliana ids are present in the ortho data - 51 thaliana genes are missing in the thaliana data
thaliana_ortho_genes = intersect(thaliana_genes, thaliana_ortho_genes)

# SAM data
COL_SAM <- COL_SAM[thaliana_ortho_genes, ]

# remove the missing genes from the data
OX_DF_1 <- OX_DF_1[thaliana_ortho_genes, ]
OX_DF_2 <- OX_DF_2[thaliana_ortho_genes, ]
OX_DF_3 <- OX_DF_3[thaliana_ortho_genes, ]
OX_DF_7 <- OX_DF_7[thaliana_ortho_genes, ]

# All gene IDs - all datasets
ortho_genes = rownames(COL_SAM)

# Remove protoplasting-induced genes from the total set of hirsuta genes
genes_to_keep = setdiff(ortho_genes, PP_genes)


##### Remove the protoplasting induced genes
OX_DF_1 <- OX_DF_1[genes_to_keep, ]
OX_DF_2 <- OX_DF_2[genes_to_keep, ]
OX_DF_3 <- OX_DF_3[genes_to_keep, ]
OX_DF_7 <- OX_DF_7[genes_to_keep, ]

# SAM data
COL_SAM <- COL_SAM[genes_to_keep, ]


###
# OX - 1 E
###

# First replicate - OX 1E - total cells 6640; filter out genes that are not detected in at least 13 cells
OX_1E <- CreateSeuratObject(counts = OX_DF_1, project = "OX_1E", min.features = 200)

# Add metadata information to the seurat object
OX_1E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-1", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_1E <- subset(OX_1E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_1E[["percent.mt"]] <- PercentageFeatureSet(OX_1E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_1E[["percent.pt"]] <- PercentageFeatureSet(OX_1E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_1E <- subset(OX_1E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_1E <- NormalizeData(OX_1E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_1E <- FindVariableFeatures(OX_1E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
OX_W1 <- GetAssayData(OX_1E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W1) <- paste("O1", colnames(OX_W1), sep = "_")


###
# OX - 2 E
###

# First replicate - OX 2E - total cells 10760; filter out genes that are not detected in at least 21 cells
OX_2E <- CreateSeuratObject(counts = OX_DF_2, project = "OX_2E", min.features = 200)

# Add metadata information to the seurat object
OX_2E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-2", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_2E <- subset(OX_2E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_2E[["percent.mt"]] <- PercentageFeatureSet(OX_2E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_2E[["percent.pt"]] <- PercentageFeatureSet(OX_2E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_2E <- subset(OX_2E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_2E <- NormalizeData(OX_2E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_2E <- FindVariableFeatures(OX_2E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
OX_W2 <- GetAssayData(OX_2E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W2) <- paste("O2", colnames(OX_W2), sep = "_")


###
# OX - 3 E
###

# First replicate - OX 3E - total cells 4100; filter out genes that are not detected in at least 8 cells
OX_3E <- CreateSeuratObject(counts = OX_DF_3, project = "OX_3E", min.features = 200)

# Add metadata information to the seurat object
OX_3E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-3", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_3E <- subset(OX_3E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_3E[["percent.mt"]] <- PercentageFeatureSet(OX_3E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_3E[["percent.pt"]] <- PercentageFeatureSet(OX_3E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_3E <- subset(OX_3E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_3E <- NormalizeData(OX_3E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_3E <- FindVariableFeatures(OX_3E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
OX_W3 <- GetAssayData(OX_3E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W3) <- paste("O3", colnames(OX_W3), sep = "_")

###
# OX - 7 E
###

# First replicate - OX 7E - total cells 9090; filter out genes that are not detected in at least 18 cells
OX_7E <- CreateSeuratObject(counts = OX_DF_7, project = "OX_7E", min.features = 200)

# Add metadata information to the seurat object
OX_7E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-7", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_7E <- subset(OX_7E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_7E[["percent.mt"]] <- PercentageFeatureSet(OX_7E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_7E[["percent.pt"]] <- PercentageFeatureSet(OX_7E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_7E <- subset(OX_7E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_7E <- NormalizeData(OX_7E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_7E <- FindVariableFeatures(OX_7E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
OX_W7 <- GetAssayData(OX_7E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W7) <- paste("O7", colnames(OX_W7), sep = "_")


###
# Intersection of highly variable genes - Cardamine hirsuta
###

# Lets find common variable features for the replicates of C. hirsuta
ox_hvgs = Reduce(intersect, list(OX_1E@assays$RNA@var.features, OX_2E@assays$RNA@var.features, OX_3E@assays$RNA@var.features, OX_7E@assays$RNA@var.features))

fileGenerator(ox_hvgs, "Shared_highly_variable_genes_between_replicates_CH.txt")


#### SAM
# First replicate - SAM 1E - total cells 4850; filter out genes that are not detected in at least 13 cells
SAM_1E <- CreateSeuratObject(counts = COL_SAM, project = "SAM_1E", min.features = 200)

# Add metadata information to the seurat object
SAM_1E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-SAM-1", "WT", "Apex")

# Remove cells with a total count more than 110000
SAM_1E <- subset(SAM_1E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
SAM_1E[["percent.mt"]] <- PercentageFeatureSet(SAM_1E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
SAM_1E[["percent.pt"]] <- PercentageFeatureSet(SAM_1E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
SAM_1E <- subset(SAM_1E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
SAM_1E <- NormalizeData(SAM_1E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
SAM_1E <- FindVariableFeatures(SAM_1E, selection.method = "vst", nfeatures = 2000)

# Select HVGs for SAM_1E
SAM_HVGs = VariableFeatures(SAM_1E)

fileGenerator(SAM_HVGs, "SAM_2000_HVGs_Seurat.txt")

###
# Combine the two sets of highly variable genes - Arabidopsis thaliana and Cardamine hirsuta
###

HVGs_combined = union(ox_hvgs, SAM_HVGs)



###
# Integrate the datasets
###

# Let's create a list to store all the seurat objects
seurat.objects <- list(O1 = OX_1E, O2 = OX_2E, O3 = OX_3E, O7 = OX_7E, S1 = SAM_1E)

# Integration of the replicates - find features from the datasets to anchor cells from different sources
anchFeatures <- SelectIntegrationFeatures(object.list = seurat.objects)

ingAnchors <- FindIntegrationAnchors(object.list = seurat.objects, dims = 1:50, anchor.features = HVGs_combined)

# To keep records of all the genes in the integrated assay, create a feature set with all of the genes from different replicates.
features_integrated <- unique(c(rownames(OX_1E), rownames(OX_2E), rownames(OX_3E), rownames(OX_7E), rownames(SAM_1E)))

# Integrate the replicates
integrated.data <- IntegrateData(anchorset = ingAnchors, dims = 1:50, verbose = T, features.to.integrate = features_integrated)

# Setting the default assay to "integrated"
DefaultAssay(integrated.data) <- "integrated"

# Gene level scaling - standardization
integrated.data <- ScaleData(integrated.data, verbose = FALSE)

# Run PCA
integrated.data <- RunPCA(integrated.data, npcs = 50, verbose = FALSE)

# Run UMAP and tSNE
integrated.data <- RunUMAP(integrated.data, reduction = "pca", dims = 1:50, n.components = 2)

integrated.data <- RunTSNE(integrated.data, reduction = "pca", dims = 1:50, dim.embed = 2)

# Find neighbours and clusters
integrated.data <- FindNeighbors(integrated.data, reduction = "pca", dims = 1:50)

for (i in seq(0.1, 1.2, 0.1)) {
  integrated.data <- FindClusters(integrated.data, resolution = i, n.start = 50, n.iter = 50)
}

# Save the integrated file as output
save(integrated.data, file = "integrated_wt_hirsta_apex_seurat.RData")

writeLines(capture.output(sessionInfo()), "Session_info_integrated_wt_hirsuta_apex_seurat.txt")
