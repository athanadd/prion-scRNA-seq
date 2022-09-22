# Load packages
library("Seurat")
library("openxlsx")
 
# Samples info
samples <- read.xlsx("samples.xlsx", detectDates = T)
 
# Create the directory for the Seurat objects
dir.create('./seurat_objects', showWarnings = F)
 
# Function to create the Seurat objects and add metadata
create_seurat_object <- function(samples) {
  timepoint <- samples[1]
  animal <- samples[2]
  inocula <- samples[3]
  date <- samples[4]
  
  project_folder <- paste0("mice_", timepoint)
  sample_folder <- paste0(animal, '_DGE_filtered')
  dge_path <- file.path('..',project_folder, 'splitseq_pipeline', 'libs_merged', sample_folder)
  
  # Read the 3 files as a sparse matrix
  sr.data <- Seurat::ReadMtx(
    mtx = file.path(dge_path, "DGE.mtx"),
    cells = file.path(dge_path, "cell_metadata.csv"),
    features = file.path(dge_path, "genes.csv"),
    feature.column = 2,
    cell.sep = ",",
    feature.sep = ",",
    mtx.transpose = T,
    skip.cell = 1,
    skip.feature = 1
  )
  
  # Create the Seurat object
  sr <- CreateSeuratObject(counts = sr.data, project = "mouse_sc", min.cells = 3)
  
  # Add metadata
  sr[["timepoint"]] <- timepoint
  sr[["animal"]] <- animal
  sr[["inocula"]] <- inocula
  sr[["date"]] <- date
  
  # Set identities for each cell
  Idents(sr) <- paste(timepoint, inocula, animal, sep = "_")
  
  # Save the object
  saveRDS(sr, file.path("seurat_objects", paste0(animal, ".rds")))
  
  # Return the object to be saved in a list
  sr
}
 
# Run the function to save the objects in the directory
sr_objects <- apply(samples, 1, create_seurat_object)
 
# Now we need to merge the objects in a new object
sr_merged <- merge(x = sr_objects[[1]], y = sr_objects[-1])
 
# Save the merged object
saveRDS(sr_merged, "./seurat_objects/sr_merged.rds")