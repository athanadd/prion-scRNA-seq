# Load libraries
library("Seurat")
library("dplyr")
library("ggplot2")
library("R.matlab")
library("openxlsx")
 
# Load Matlab object
data <- readMat("./GSM3017261_150000_CNS_nuclei.mat")
 
# Keep the matrix that contains the expression values
expression.mat <- t(data$DGE)
 
# Rename the clusters to remove trailing space
cluster.assignment <- sapply(data$cluster.assignment[,1], trimws)
organ <- sapply(data$sample.type[,1], trimws)
barcode <- paste0("Cell-", data$barcodes[1,])
 
# Prepare the cell metadata
coldata <- data.frame(barcode = barcode, organ = organ, cluster_full_name = cluster.assignment)
coldata$cluster_number <- unlist(lapply(sapply(coldata$cluster_full_name, strsplit, split = " ", fixed = T), function(x) as.integer(x[1])))
 
# Add extra info to the coldata object
extra.info <- read.xlsx("./splitseq_clusters_no_unknown.xlsx")
extra.info$cluster_full_name <- NULL
coldata <- coldata %>% left_join(extra.info, by = "cluster_number")
row.names(coldata) <- coldata$barcode
 
# Prepare the genes
genes <- sapply(data$genes[,1], trimws)
 
# Add info to matrix
colnames(expression.mat) <- coldata$barcode
rownames(expression.mat) <- genes
 
# Create Seurat object
sr <- CreateSeuratObject(counts = expression.mat, project = "splitseq_paper", min.cells = 3)
sr <- AddMetaData(sr, coldata)
 
# Cleanup
rm(coldata, data, expression.mat, extra.info, cluster.assignment, genes, organ, barcode)
 
# Remove clusters from tissue not relevant to my study
sr <- subset(sr, subset = keep == "yes")
sr$keep <- NULL
 
# Split the object and keep the P2 and P11 brain only
sr.p2 <- subset(sr, subset = organ == "p2_brain")
sr.p11 <- subset(sr, subset = organ == "p11_brain")
sr.list <- list(sr.p2, sr.p11)
 
# Cleanup
rm(sr, sr.p2, sr.p11)
 
# Normalize datasets individually by SCTransform()
sr.list <- lapply(X = sr.list, FUN = SCTransform, method = "glmGamPoi")
 
# Select the integration features
features <- SelectIntegrationFeatures(object.list = sr.list, nfeatures = 3000)
 
# Run the PrepSCTIntegration() function prior to identifying anchors
sr.list <- PrepSCTIntegration(object.list = sr.list, anchor.features = features)
 
# When running FindIntegrationAnchors(), and IntegrateData(),
# set the normalization.method parameter to the value SCT.
int.anchors <- FindIntegrationAnchors(
  object.list = sr.list,
  normalization.method = "SCT",
  anchor.features = features)
 
sr_integrated <- IntegrateData(
  anchorset = int.anchors,
  normalization.method = "SCT")
 
# Cleanup
sr_ref <- sr_integrated
rm(sr.list, features, int.anchors, sr_integrated)
 
# Run PCA and UMAP on the data
sr_ref <- sr_ref %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30, return.model = TRUE)
 
# Save the pre-processed integrated reference
saveRDS(sr_ref, "./sr_integrated_reference.rds")
