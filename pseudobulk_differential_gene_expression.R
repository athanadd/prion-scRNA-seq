# Load libraries
library("Seurat")
library("ggplot2")
library("tidyverse")
library("openxlsx")
library("SingleCellExperiment")
library("DESeq2")
 
#### Pseudobulk analysis ####
timepoints <- c("20dpi", "40dpi", "80dpi", "120dpi", "end")
 
degs_list <- list()
 
for (i in seq_along(timepoints)) {
 
  # Set the timepoint variable
  TIMEPOINT <- timepoints[i]
  
  # Load the Seurat file for the timepoint
  sr <- readRDS(paste0("../seurat_objects/", TIMEPOINT, "_sr_renamed_qc_integrated_annotated_filtered.rds"))
  
  # Keep only the CD1 and RML samples
  sr_RML_CD1 <- subset(sr, subset = inocula %in% c("RML", "CD1"))
  
  # Cleanup
  rm(sr)
  
  # Function to run DESeq2 for each cluster and generate
  # relevant plots
  run_DESeq <- function(current_cluster, seurat_obj) {
    
    print(paste0("Working on cluster: ", current_cluster))
    
    # Subset again to select the cluster of interest
    sr_cluster <- subset(seurat_obj, subset = cluster_full_name == current_cluster)
    
    # Convert Seurat object to SingleCellExperiment
    sce <- as.SingleCellExperiment(sr_cluster)
    
    # Convert characters to factors
    sce$animal <- factor(sce$animal)
    sce$inocula <- factor(sce$inocula, levels = c("CD1", "RML"))
    
    # Count aggregation to sample level
    sce_agg <- Matrix.utils::aggregate.Matrix(t(counts(sce)), 
                           groupings = sce$animal,
                           fun = "sum") 
    
    # Transpose the matrix
    sce_agg <- t(sce_agg)
    
    # Prepare the metadata
    sce_metadata <- data.frame(animal = as.numeric(colnames(sce_agg))) %>%
      left_join(read.xlsx("../samples.xlsx"), by = "animal") %>%
      column_to_rownames("animal")
    
    # Build the DESeq2 object
    dds <- DESeqDataSetFromMatrix(sce_agg, 
                                  colData = sce_metadata, 
                                  design = ~ inocula)
    
    
    # Transform counts for data visualization
    rld <- rlog(dds, blind=TRUE)
    
    # Plot PCA
    dir.create("PCA_plots", showWarnings = F)
    
    pca_plot <- DESeq2::plotPCA(rld, intgroup = "inocula")
    ggsave(paste0("./PCA_plots/", TIMEPOINT, "_", gsub("/", "-", current_cluster, fixed = T), ".png"), pca_plot, width = 10, height = 10)
    
    # Run the DESeq2 pipeline
    dds <- DESeq(dds)
    
    # Get the results
    res <- results(dds, 
                   contrast = c("inocula", "RML", "CD1"),
                   alpha = 0.05)
    
    # Shrink lfc
    res <- lfcShrink(dds, 
                     coef = "inocula_RML_vs_CD1",
                     res = res)
    
    # Significant DE genes
    res_sig <- data.frame(res) %>%
      filter(padj < 0.05) %>%
      arrange(padj) %>%
      rownames_to_column("gene")
    
    
    # Heatmap of the significant genes
    if (nrow(res_sig) >= 2) {
      save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
        png(filename, width = width, height = height, res = res)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
      }
      
      dir.create("DEGs_heatmaps", showWarnings = F)
      
      # Extract normalized counts for only the significant genes
      sig_norm <- data.frame(counts(dds, normalized = TRUE)) %>%
        rownames_to_column(var = "gene") %>%
        dplyr::filter(gene %in% res_sig$gene) %>%
        select(-gene)
      
      hm_anno <- sce_metadata[,"inocula", drop = F]
      row.names(hm_anno) <- colnames(sig_norm)
      
      # Run pheatmap using the metadata data frame for the annotation
      hm <- pheatmap::pheatmap(sig_norm, 
               color = RColorBrewer::brewer.pal(6, "YlOrRd"),
               border_color = NA,
               cluster_rows = T, 
               show_rownames = F,
               annotation = hm_anno, 
               scale = "row")   
      
      save_pheatmap_png(hm, paste0("./DEGs_heatmaps/", TIMEPOINT, "_", gsub("/", "-", current_cluster, fixed = T), ".png"))
    }
    
    # Return the results
    if(!is.null(res_sig) && nrow(res_sig) > 0) {
      cbind(cluster = current_cluster, res_sig)
    }
  }
  
  # Run DGE on all clusters
  all_clusters <- unique(sr_RML_CD1$cluster_full_name)
  degs <- map_dfr(all_clusters, run_DESeq, seurat_obj = sr_RML_CD1)
  
  degs_list[[TIMEPOINT]] <- degs
}
 
# Cleanup
rm(i, degs, sr_RML_CD1)
 
# Filter out CD1 vs PBS genes
genes_to_exclude <- c("Calm1", "Cdk8", "Cmss1", "Malat1", "mt-Rnr1", "mt-Rnr2", "Rn18s")
list_subt <- lapply(degs_list, function(x) x[!x$gene %in% genes_to_exclude,])
 
names(list_subt) <- timepoints
 
# Save as xlsx
write.xlsx(degs_list, "DEGs_DESeq2.xlsx", overwrite = TRUE)
write.xlsx(list_subt, "DEGs_DESeq2_subtracted_v2.xlsx", overwrite = TRUE)
 
#### PCA plots accross all time points ####
 
dir.create("PCA_plots/all_timepoints", showWarnings = F, recursive = T)
 
timepoints <- c("20dpi", "40dpi", "80dpi", "120dpi", "end")
 
load_sr_objects <- function(TIMEPOINT) {
  
  # Load the Seurat file for the timepoint
  sr <- readRDS(paste0("../seurat_objects/", TIMEPOINT, "_sr_renamed_qc_integrated_annotated_filtered.rds"))
  
  # Keep only the CD1 and RML samples
  sr_RML_CD1 <- subset(sr, subset = inocula %in% c("RML", "CD1"))
  
  # Cleanup
  rm(sr)
  
  return(sr_RML_CD1)
}
 
# Load all time points in memory
srs <- sapply(timepoints, load_sr_objects)
 
# Merge to create a combined object
sr_merged <- merge(srs[[1]], y = srs[-1], project = "mouse_sc")
 
# Remove assays and save the merged
DefaultAssay(sr_merged) <- "RNA"
sr_merged[["SCT"]] <- NULL
sr_merged[["integrated"]] <- NULL
sr_merged[["prediction.score.cluster_number"]] <- NULL
 
# Remove spurious sample 828719
sr_merged <- subset(sr_merged, subset = animal != "828719")
# Remove 40 and 80 dpi because there are no interesting transcriptomic changes
sr_merged <- subset(sr_merged, subset = timepoint %in% c("20dpi", "120dpi", "end"))
saveRDS(sr_merged, "sr_merged_for_PCA_plots.rds")
 
# Cleanup
rm(srs, load_sr_objects)
 
# Load the file for future use
#sr_merged <- readRDS("sr_merged_for_PCA_plots.rds")
 
# Prepare a vector with all clusters
all_clusters <- unique(sr_merged$cluster_full_name)
 
# Modified function from the DESeq2 visualisation functions
# to allow specification of different labeling for the timepoints
# adapted from https://github.com/mikelove/DESeq2/blob/master/R/plots.R
plotPCA_custom <- function(object, ntop=500, returnData=FALSE) {
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  intgroup.df <- as.data.frame(colData(object)[, c("inocula", "timepoint"), drop=FALSE])
  intgroup.df$animal <- row.names(colData(object))
  intgroup.df$inocula <- factor(intgroup.df$inocula, levels = c("CD1", "RML"))
  intgroup.df$timepoint <- factor(intgroup.df$timepoint, levels = timepoints)
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], intgroup.df, name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  # Set custom shapes
  custom_shapes <- c("20dpi" = 15, "120dpi" = 16, "end" = 17)
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="inocula", shape = "timepoint", label="animal")) +
    geom_point(size=3) +
    #geom_text_repel(max.overlaps = 20) +
    scale_shape_manual(values = custom_shapes) +
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed()
}
 
# Function to run DESeq2 for each cluster and generate
# relevant plots
generate_PCA_plot <- function(current_cluster, seurat_obj) {
  
  print(paste0("Working on cluster: ", current_cluster))
  
  # Subset to select the cluster of interest
  sr_cluster <- subset(seurat_obj, subset = cluster_full_name == current_cluster)
  
  # Convert Seurat object to SingleCellExperiment
  sce <- as.SingleCellExperiment(sr_cluster)
  
  # Convert characters to factors
  sce$animal <- factor(sce$animal)
  
  # Count aggregation to sample level
  sce_agg <- Matrix.utils::aggregate.Matrix(t(counts(sce)), 
                                            groupings = sce$animal,
                                            fun = "sum") 
  
  # Transpose the matrix
  sce_agg <- t(sce_agg)
  
  # Prepare the metadata
  sce_metadata <- data.frame(animal = as.numeric(colnames(sce_agg))) %>%
    left_join(read.xlsx("../samples.xlsx"), by = "animal") %>%
    column_to_rownames("animal")
  
  # Build the DESeq2 object
  dds <- DESeqDataSetFromMatrix(sce_agg, 
                                colData = sce_metadata, 
                                design = ~ inocula)
  
  # Transform counts for data visualization
  rld <- tryCatch({
    vst(dds)
  }, error = function(cond) varianceStabilizingTransformation(dds))
  
  # Plot PCA
  pca_plot <- plotPCA_custom(rld) + ggtitle(current_cluster)
  ggsave(paste0("PCA_plots/all_timepoints/", gsub("/", "-", current_cluster, fixed = T), ".png"), width = 6, height = 4)
}
 
sapply(all_clusters, generate_PCA_plot, seurat_obj = sr_merged)