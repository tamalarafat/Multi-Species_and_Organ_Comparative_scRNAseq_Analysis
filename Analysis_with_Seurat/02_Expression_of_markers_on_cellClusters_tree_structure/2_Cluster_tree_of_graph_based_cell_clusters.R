# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)


# Load the known cell type markers file
known_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Thaliana_known_cell_type_markers.csv")

# Load
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_Cardamine_and_Arabidopsis_Apex/Seurat_Analyses/Seurat_Analysis_strategy_2/01_Integration_of_data/integrated_wt_hirsta_apex_seurat.RData")

DefaultAssay(integrated.data) <- "RNA"

# Storing directory
storing_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_Cardamine_and_Arabidopsis_Apex/Seurat_Analyses/Seurat_Analysis_strategy_2"

clus_tree_of_cell_clusters(seuratObject = integrated.data, marker_file = known_markers, query_pattern = "integrated_snn_res.", gene_ID_column = "AT_ID", gene_name_column = "Gene_Name", store_dir = storing_dir)
