# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_Cardamine_and_Arabidopsis_Apex/Liger_Analyses/Liger_Analysis_strategy_2/On_coefficient/Seurat_objects/seurat_object_of_K_50.RData")

Idents(integrated.data) <- "RNA_snn_res.0.2"

cluster_degs = FindConservedMarkers(integrated.data, ident.1 = "13", grouping.var = "Species", test.use = "wilcox", only.pos = TRUE)
cluster_degs$gene_ID = rownames(cluster_degs)

save(cluster_degs, file = "Cluster_13_all_degs.RData")