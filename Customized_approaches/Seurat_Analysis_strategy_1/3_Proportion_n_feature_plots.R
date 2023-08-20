# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

# Load
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_Cardamine_and_Arabidopsis_Apex/Seurat_Analyses/Seurat_Analysis_strategy_1/01_Integration_of_data/integrated_wt_hirsta_apex_seurat.RData")

# Storing directory
storing_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_Cardamine_and_Arabidopsis_Apex/Seurat_Analyses/Seurat_Analysis_strategy_1"

cell_type_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Curated_markers_thaliana.csv")

resolution_explorer(seuratObject = integrated.data, store_dir = storing_dir, store_folder = "Results_for_a_range_of_resolution")
