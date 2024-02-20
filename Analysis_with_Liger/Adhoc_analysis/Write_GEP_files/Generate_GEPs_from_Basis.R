# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
projects_dir <- "/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/"

# Functions - Markers identification
list_1 <- list.files(paste0(projects_dir, "Library_handler"), pattern = "*.R$", full.names = TRUE) 
sapply(list_1, source, .GlobalEnv)

# Functions - Data manipulation
list_2 <- list.files(paste0(projects_dir, "Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE) 
sapply(list_2, source, .GlobalEnv)

# Functions - Library and packages handler
list_3 <- list.files(paste0(projects_dir, "Functions_marker_identification"), pattern = "*.R$", full.names = TRUE) 
sapply(list_3, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the file containing grouping of the genes from the basis matrix for each factorization
load("/home/ytamal2/Documents/2023/PhD_projects_Yasir/comparative_study_Cardamine_and_Arabidopsis_Apex/Analysis_with_Liger/Analysis_objects/Gene_clusters/Basis.RData")

# Directory path -  store the outputs on this directory
storing_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/comparative_study_Cardamine_and_Arabidopsis_Apex/Analysis_with_Liger/Analysis_output"

Basis_to_GEP_generator(input_data = Basis, 
                       Factor_ID = "50", 
                       store_GEPs = TRUE, 
                       return_GEPs = FALSE, 
                       generate_ortho = FALSE, 
                       store_dir = storing_dir)
