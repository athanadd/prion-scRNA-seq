# Load packages
library("Seurat")
library("openxlsx")
library("ggplot2")
library("ggrepel")
library("RColorBrewer")
library("tidyverse")
library("ensembldb")
library("AnnotationHub")
library("scProportionTest")
library("clusterProfiler")
library("cowplot")
 
# Set the time point for the whole script
TIMEPOINT <- "20dpi"
 
#### Create Seurat object ####
 
# Samples info
samples <- read.xlsx("samples.xlsx", detectDates = T)
 
# Subset the samples
samples <- subset(samples, subset = timepoint == TIMEPOINT)
 
# Function to load the Seurat objects and add metadata
load_seurat_object <- function(samples) {
  timepoint <- samples[1]
  animal <- samples[2]
  inocula <- samples[3]
  date <- samples[4]
  
  # Create the Seurat object
  sr <- readRDS(paste0("./seurat_objects/", animal, ".rds"))
  
  # Add metadata
  sr[["timepoint"]] <- timepoint
  sr[["animal"]] <- animal
  sr[["inocula"]] <- inocula
  sr[["date"]] <- date
  
  # Set identities for each cell
  Idents(sr) <- paste(timepoint, inocula, animal, sep = "_")
  
  # Return the object to be saved in a list
  sr
}
 
# Run the function to save the objects in the directory
sr_objects <- apply(samples, 1, load_seurat_object)
 
# Now we need to merge the objects in a new object
sr_merged <- merge(x = sr_objects[[1]], y = sr_objects[-1])
 
# Cleanup
sr <- sr_merged
rm(samples, load_seurat_object, sr_objects, sr_merged)
 
 
#### Rename features ####
 
## Load the annotation resource.
ah <- AnnotationHub()
 
# fetch one of the databases
# ahDb <- query(ah, pattern = c("Mus musculus", "EnsDb", 104))
ahEdb <- ah[["AH95775"]]
 
# Create one vector with the Ensembl IDs of all the genes from the experiment
ensembl.genes <- row.names(sr[["RNA"]])
 
# Convert the Ensembl IDs to Gene symbols
gene_ids <- ensembldb::select(ahEdb, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
 
# Some Ensembl IDs don't have corresponding gene symbols and will be removed
empty_genes <- which(gene_ids$SYMBOL == "")
gene_ids <- gene_ids[-empty_genes,]
sr <- sr[-empty_genes,]
 
# There might be duplicate names in the symbols. Add a suffix to make them unique
unique_gene_ids <- make.unique(gene_ids$SYMBOL)
 
# Replace underscores with dashes because underscores are not allowed in Seurat
unique_gene_ids <- gsub("_", "-", unique_gene_ids, fixed = T, )
 
# Function to remove rows that could not be matched and rename
# the RNA assay slot of the Seurat object
RenameGenesSeurat <- function(obj, newnames) {
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    
    # Rename features
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    
  } else {
    stop("Unequal gene sets: nrow(RNA) != nrow(newnames)")
  }
  obj@assays$RNA <- RNA
  
  # Fix the row.names in meta.features
  row.names(obj[["RNA"]]@meta.features) <- row.names(obj[["RNA"]])
  return(obj)
}
 
# Prepare the renamed object
sr_renamed <- RenameGenesSeurat(sr, unique_gene_ids)
 
# Save
saveRDS(sr_renamed, paste0("./seurat_objects/", TIMEPOINT, "_sr_renamed.rds"))
 
# Cleanup
sr <- sr_renamed
rm(ah, ahDb, ahEdb, gene_ids, empty_genes, ensembl.genes, unique_gene_ids, RenameGenesSeurat, sr_renamed)
 
#### QC ####
 
# Plot number of features and counts
VlnPlot(sr, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0)
 
# Calculate mitochondrial genes percentage
sr[["percent.mt"]] <- PercentageFeatureSet(sr, pattern = "^mt-")
 
# Filter cells with fewer than 200 expressed genes or more than 2500
# Filter cells that have >1% mitochondrial counts
sr_qc <- subset(sr, subset = nFeature_RNA > 250 & nFeature_RNA < 2500 & percent.mt < 1)
 
# Add cell cycle genes information
# Basic function to convert human to mouse gene names
convert_genes_human_to_mouse <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  mouse_genes <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(mouse_genes))
  
  mouse_genes
}
 
s.genes <- convert_genes_human_to_mouse(cc.genes.updated.2019$s.genes)
g2m.genes <- convert_genes_human_to_mouse(cc.genes.updated.2019$g2m.genes)
 
# Check if cells separate by cell cycle phase
dir.create("plots/cell_cycle", recursive = T, showWarnings = F)
 
sr_phase <- sr_qc
sr_phase <- CellCycleScoring(sr_phase, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
sr_phase <- NormalizeData(sr_phase)
sr_phase <- FindVariableFeatures(sr_phase)
sr_phase <- ScaleData(sr_phase, features = rownames(sr_phase))
sr_phase <- RunPCA(sr_phase, features = c(s.genes, g2m.genes))
DimPlot(sr_phase, shuffle = TRUE) +
  ggtitle(paste0(TIMEPOINT, " cell cycle PCA"))
ggsave(paste0("plots/cell_cycle/", TIMEPOINT, "_PCA.png"), width = 8, height = 4)
 
# Assign cell cycle scores
sr_qc <- CellCycleScoring(sr_qc, s.features = s.genes, g2m.features = g2m.genes)
 
# Plot number of features and counts
VlnPlot(sr_qc, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0)
 
# Save
saveRDS(sr_qc, paste0("./seurat_objects/", TIMEPOINT, "_sr_renamed_qc.rds"))
 
# Cleanup
sr <- sr_qc
rm(sr_qc, s.genes, g2m.genes, sr_phase, convert_genes_human_to_mouse)
 
 
#### Integrate datasets ####
 
# Split the dataset into a list of two seurat objects based on inocula
sr.list <- SplitObject(sr, split.by = "inocula")
 
# Remove the PBS group
#sr.list <- sr.list[-1]
 
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
 
# Save
saveRDS(sr_integrated, paste0("./seurat_objects/", TIMEPOINT, "_sr_renamed_qc_integrated.rds"))
 
# Cleanup
sr <- sr_integrated
rm(sr.list, features, int.anchors, sr_integrated)
 
#### Annotation (reference integration) ####
 
# Load pre-processed reference
sr_ref <- readRDS("../splitseq_paper_reference/sr_integrated_reference.rds")
 
query <- sr
rm(sr)
 
# Find the transfer anchors between the two datasets
query.anchors <- FindTransferAnchors(
  reference = sr_ref,
  query = query,
  dims = 1:50,
  reference.reduction = "pca",
  normalization.method = "SCT"
)
 
# Make a vector with the metadata to transfer from the reference to the query
labels_to_transfer <- list(
  cluster_number = "cluster_number"
)
 
# Do the transfer
query <- TransferData(
  reference = sr_ref,
  query = query,
  anchorset = query.anchors,
  refdata = labels_to_transfer,
  dims = 1:50
)
 
# Add extra info to the coldata object
extra.info <- read.xlsx("../splitseq_paper_reference/splitseq_clusters.xlsx")
 
df <- data.frame(cluster_number = as.integer(query$predicted.cluster_number)) %>%
  left_join(extra.info, by = "cluster_number")
 
df$cluster_number <- NULL
df$keep <- NULL
 
row.names(df) <- colnames(query)
 
query <- AddMetaData(query, df)
 
# Calculate mapping score and add to metadata
query <- AddMetaData(
  object = query,
  metadata = MappingScore(anchors = query.anchors),
  col.name = "mapping.score"
)
 
# Cleanup
rm(df, extra.info, labels_to_transfer, query.anchors)
 
# Assess the score of the predictions
ggplot() + aes(query$predicted.cluster_number.score) + geom_histogram()
ggplot() + aes(query$mapping.score) + geom_histogram()
 
# Filter the query object based on label transfer quality
query.filt <- subset(query, region != "Olfactory Bulb")
 
# Number of cells in each cluster
cells_per_cluster <- query.filt@meta.data %>%
  group_by(cluster_full_name) %>%
  count() %>%
  arrange(n)
 
# Select the clusters with fewer than 100 cells
clusters_to_keep <- cells_per_cluster[cells_per_cluster$n > 100,]
clusters_to_keep <- clusters_to_keep$cluster_full_name
 
# Filter query to remove clusters with fewer than 100 cells
query.filt <- subset(query.filt, subset = cluster_full_name %in% clusters_to_keep)
 
# Re-cluster the query
query.filt <- query.filt %>%
  SCTransform(method = "glmGamPoi") %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)
 
# Save
saveRDS(query, paste0("./seurat_objects/", TIMEPOINT, "_sr_renamed_qc_integrated_annotated.rds"))
saveRDS(query.filt, paste0("./seurat_objects/", TIMEPOINT, "_sr_renamed_qc_integrated_annotated_filtered.rds"))
 
 
# Cleanup
sr <- query.filt
rm(query, query.filt, clusters_to_keep)
 
# Sanity check after annotation
# Check that the annotated clusters use the gene markers
dir.create("plots/marker_genes", showWarnings = F, recursive = T)
 
# Relevel clusters by cluster id so that the plots are nicer
sr$cluster_full_name <- factor(sr$cluster_full_name,
                               levels = unique(sr$cluster_full_name)[order(as.integer(str_extract(unique(sr$cluster_full_name), "^\\d+")),
                                                                           decreasing = T)])
 
VlnPlot(sr, features = c("Gria1", "Snhg11", "Mbp", "Plp1", "Vcan", "Dock8", "Flt1", "Slc1a2", "Plpp3", "Dnah11"),
        pt.size = 0, stack = T, group.by = "cluster_full_name")
 
ggsave(paste0("plots/marker_genes/", TIMEPOINT, "_marker_genes.png"), width = 12, height = 8)
 
# Astro: Aqp4, Slc1a2, Plpp3, Gja1
# Oligodendrocytes: Mbp, Plp1
# Oligodendrocyte Precursor Cells: Vcan & Mbp, Pdgfra
# Endothelial/smooth muscle Cells: Rgs5, Flt1, Ly6c1, Pltp
# Microglia/macrophages: Dock2, Dock8, Csf1r, P2ry12
# Ependymal cells: Dnah11
# Neurons: Gria1, Snhg11?
 
dir.create("./cluster_metrics", showWarnings = F)
 
# Number of cells in each cluster
cells_per_cluster <- sr@meta.data %>%
  group_by(cluster_full_name) %>%
  count() %>%
  arrange(n)
write.table(cells_per_cluster, paste0("./cluster_metrics/", TIMEPOINT, "_cells_per_cluster.tsv"))
 
# Number of cells in each group
cells_per_group <- data.frame(table(sr$group, sr$inocula, sr$animal))
names(cells_per_group) <- c("group", "inocula", "animal", "n_cells")
write.table(cells_per_group, paste0("./cluster_metrics/", TIMEPOINT, "_cells_per_group.tsv"))
 
# Calculate mean and sd of number of features per cluster
nFeatures_per_cluster <- sr@meta.data %>%
  group_by(cluster_full_name) %>%
  summarise_at(vars(nFeature_RNA ),list(mean = ~round(mean(.),0), median = median, sd = ~round(sd(.),0))) %>%
  arrange(median)
 
# Calculate number of counts per cluster
nCounts_per_cluster <- sr@meta.data %>%
  group_by(cluster_full_name) %>%
  summarise_at(vars(nCount_RNA ),list(mean = ~round(mean(.),0), median = median, sd = ~round(sd(.),0))) %>%
  arrange(median)
 
extra_metrics <- nCounts_per_cluster %>%
  left_join(nFeatures_per_cluster, by = "cluster_full_name")
colnames(extra_metrics) <- c("Cluster name", "Counts mean", "Counts median", "Counts SD",
                             "Features mean", "Features median", "Features SD")
write.xlsx(extra_metrics, paste0("./cluster_metrics/", TIMEPOINT, "_extra_metrics.xlsx"), overwrite = T)
 
# Plots
dir.create("./plots/reduced_dimensions", showWarnings = F, recursive = T)
 
DimPlot(sr, group.by = "cluster_full_name", label = T, repel = T) +
  ggtitle(TIMEPOINT)
ggsave(paste0("./plots/reduced_dimensions/", TIMEPOINT, "_UMAP.png"), width = 16, height = 10)
 
DimPlot(sr, group.by = "cluster_full_name", split.by = "inocula") +
  ggtitle(TIMEPOINT)
ggsave(paste0("./plots/reduced_dimensions/", TIMEPOINT, "_split_inocula_UMAP.png"), width = 16, height = 10)
 
 
#### Cell type proportions ####
 
pt <- table(sr$group, sr$inocula)
pt <- as.data.frame(pt)
colnames(pt) <- c("Cell type", "Experimental group", "Frequency")
pt$`Cell type` <- as.character(pt$`Cell type`)
 
dir.create("./plots/celltype_proportions", showWarnings = F, recursive = T)
 
myColors <- brewer.pal(10, "Set3")
names(myColors) <- c("Migrating Interneurons",
                     "Cortical Neurons",
                     "Medium Spiny Neurons",
                     "Astrocytes",
                     "OPC",
                     "Oligodendrocytes",
                     "VLMC",
                     "Ependymal",
                     "Immune",
                     "Vascular")
 
ggplot(pt, aes(x = `Experimental group`, y = Frequency, fill = `Cell type`)) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(name = "Cell type", values = myColors) +
  ylab("Proportion") +
  ggtitle(TIMEPOINT)
ggsave(paste0("./plots/celltype_proportions/", TIMEPOINT, "_cell_proportions.png"), width = 8, height = 4)
 
## Plots using the permutation test
prop_test <- sc_utils(sr)
prop_test <- permutation_test(
  prop_test,
  cluster_identity = "group",
  sample_1 = "CD1",
  sample_2 = "RML",
  sample_identity = "inocula"
)
permutation_plot(prop_test, log2FD_threshold = log2(1.2)) +
  ylab("Log2-fold difference in cell numbers") +
  xlab("Cell type") +
  ggtitle(TIMEPOINT)
ggsave(paste0("./plots/celltype_proportions/", TIMEPOINT, "_cell_proportions_scPropTest.png"), width = 8, height = 4)
 
# Repeat for CD1 vs PBS
prop_test <- sc_utils(sr)
prop_test <- permutation_test(
  prop_test,
  cluster_identity = "group",
  sample_1 = "PBS",
  sample_2 = "CD1",
  sample_identity = "inocula"
)
permutation_plot(prop_test, log2FD_threshold = log2(1.2)) +
  ylab("Log2-fold difference in cell numbers") +
  xlab("Cell type") +
  ggtitle(TIMEPOINT)
ggsave(paste0("./plots/celltype_proportions/", TIMEPOINT, "_PBS_vs_CD1_cell_proportions_scPropTest.png"), width = 8, height = 4)
 
# Test only groups of neurons
sr_neurons <- subset(sr, subset = group %in% c("Medium Spiny Neurons", "Cortical Neurons", "Migrating Interneurons"))
 
prop_test <- sc_utils(sr_neurons)
prop_test <- permutation_test(
  prop_test,
  cluster_identity = "cluster_full_name",
  sample_1 = "CD1",
  sample_2 = "RML",
  sample_identity = "inocula"
)
# Reorder the data to have the clusters in order for the plot
prop_test@results$permutation$clusters <- 
  factor(prop_test@results$permutation$clusters,
         levels = prop_test@results$permutation$clusters[order(as.integer(str_extract(prop_test@results$permutation$clusters, "^\\d+")), decreasing = T)])
 
permutation_plot(prop_test, log2FD_threshold = log2(1.2), order_clusters = F) +
  ylab("Log2-fold difference in cell numbers") +
  xlab("Cell cluster") +
  ggtitle(TIMEPOINT)
ggsave(paste0("./plots/celltype_proportions/", TIMEPOINT, "_neurons_cell_proportions_scPropTest.png"), width = 8, height = 4)
 
rm(prop_test, pt, myColors, sr_neurons)
 
#### DGE ####
# Function to perform DGE between clusters and two conditions
get_DEGs <- function(cluster, condition1, condition2, seurat_obj){
  genes <- tryCatch(
    {
      FindMarkers(seurat_obj,
                       ident.1 = paste0(cluster, "_", condition1),
                       ident.2 = paste0(cluster, "_", condition2)) %>%
      rownames_to_column(var = "gene")
    }, error = function(cond) return (NULL))
  
  if(!is.null(genes) && nrow(genes) > 0) {
    cbind(cluster = cluster, genes)
  }
}
 
# Prepare a vector of all all cluster names
all_clusters <- unique(as.character(sr$cluster_full_name))
 
# Add new Idents to Seurat object
Idents(sr) <- paste0(sr$cluster_full_name, "_", sr$inocula)
 
# Run DGE on all clusters
degs <- map_dfr(all_clusters, get_DEGs, condition1 = "RML", condition2 = "CD1", seurat_obj = sr)
 
# Plot the distribution of the p-values
ggplot(degs, aes(p_val_adj)) + geom_histogram()
 
# Keep DEGs with adjusted p-values < 0.05
degs.filtered <- subset(degs, subset = p_val_adj < 0.05)
 
# Add info if DEG is unique in each cluster
add_unique_info <- function(row, degs.filtered) {
  current_cluster <- row["cluster"]
  genes <- subset(degs.filtered, subset = cluster != current_cluster)[, "gene"]
  gene <- row["gene"]
  
  return(!gene %in% genes)
}
degs.filtered$gene_unique <- apply(degs.filtered, 1, add_unique_info, degs.filtered = degs.filtered)
 
# Number of DEGs in each cluster
table(degs.filtered$cluster)
 
# Bar chart to visualize the number of DEGs in each cluster
dir.create("./plots/DGE", showWarnings = F, recursive = T)
ggplot(as.data.frame(table(degs.filtered$cluster)), aes(Var1, Freq)) + geom_col() + coord_flip() +
  xlab("Clusters") + ylab("Number of DEGs") + ggtitle(paste0(TIMEPOINT, " number of DEGs (adj_p_val < 0.05)"))
ggsave(paste0("./plots/DGE/", TIMEPOINT, "_number_of_DEGs_by_cluster.png"), width = 12, height = 8)
 
 
# Plot number of DEGs vs number of cells in cluster
tb <- cells_per_cluster %>% right_join(as.data.frame(table(degs.filtered$cluster)), by = c("cluster_full_name" = "Var1"))
 
ggplot(tb, aes(n, Freq, label=cluster_full_name)) +
  geom_point() +
  xlab("number of cells") +
  ylab("number of DEGs") +
  geom_text_repel(max.overlaps = 20)
ggsave(paste0("./plots/DGE/", TIMEPOINT, "_number_of_DEGs_vs_cells.png"), width = 12, height = 12)
 
# Save gene list
dir.create("./DGE_gene_lists", showWarnings = F, recursive = T)
write.xlsx(degs.filtered, paste0("./DGE_gene_lists/", TIMEPOINT, "_DEGs_by_cluster.xlsx"), overwrite = T)
 
# Compare CD1 vs PBS
degs_contr <- map_dfr(all_clusters, get_DEGs, condition1 = "CD1", condition2 = "PBS", seurat_obj = sr)
degs.filtered_contr <- subset(degs_contr, subset = p_val_adj < 0.05)
write.xlsx(degs.filtered_contr, paste0("./DGE_gene_lists/", TIMEPOINT, "_CD1_vs_PBS_DEGs_by_cluster.xlsx"), overwrite = T)