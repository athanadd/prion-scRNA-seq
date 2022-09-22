# Load packages
library("Seurat")
library("openxlsx")
library("ggplot2")
library("ggrepel")
library("RColorBrewer")
library("tidyverse")
library("ensembldb")
library("AnnotationHub")
library("clusterProfiler")
 
# Load the Seurat object
TIMEPOINT <- "20dpi"
sr <- readRDS(paste0("./seurat_objects/", TIMEPOINT, "_sr_renamed_qc_integrated_annotated_filtered.rds"))
 
# Load the annotation resource.
ah <- AnnotationHub()
 
# fetch one of the databases
ahOrgDb <- ah[["AH92582"]]
 
# Prepare a vector of all all cluster names
all_clusters <- unique(as.character(sr$cluster_full_name))
 
## ORA - Over-Representation Analysis
 
# Create directory
dir.create("cluster_profiler/ORA", showWarnings = F, recursive = T)
 
run_ORA <- function(cluster_, degs.filtered, seurat_obj, ontology){
  genes <- subset(degs.filtered, subset = cluster == cluster_)$gene
  if (length(genes) == 0) {
    return()
  }
  print(paste0("working on cluster: ", cluster_))
  ego <- enrichGO(gene          = genes,
                  universe      = row.names(seurat_obj),
                  OrgDb         = ahOrgDb,
                  keyType       = "SYMBOL",
                  ont           = ontology,
                  pAdjustMethod = "BH")
  ego <- head(ego)
  
  if(!is.null(ego) && nrow(ego) > 0) {
    cbind(cluster = cluster_, ego)
  }
}
run_ORA_list <- function(cluster_, degs.filtered, seurat_obj, ontology, count_cutoff){
  genes <- subset(degs.filtered, subset = cluster == cluster_)$gene
  if (length(genes) == 0) {
    return()
  }
  print(paste0("working on cluster: ", cluster_))
  ego <- enrichGO(gene          = genes,
                  universe      = row.names(seurat_obj),
                  OrgDb         = ahOrgDb,
                  keyType       = "SYMBOL",
                  ont           = ontology,
                  pAdjustMethod = "BH")
  if(!is.null(ego)) {
    ego@result <- ego@result[ego@result$Count >= count_cutoff,]
    return(ego)
  }
}
 
ora.BP_list <- lapply(all_clusters,
                 run_ORA_list,
                degs.filtered = degs.filtered,
                seurat_obj = sr,
                ontology = "BP",
                count_cutoff = 3)
names(ora.BP_list) <- all_clusters
dotplot(merge_result(ora.BP_list), font.size = 12, title = paste0("ORA - BP - ", TIMEPOINT)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(paste0("cluster_profiler/ORA/", "ORA_BP_", TIMEPOINT, ".png"), width = 8, height = 6)
 
ora.CC_list <- lapply(all_clusters,
                      run_ORA_list,
                      degs.filtered = degs.filtered,
                      seurat_obj = sr,
                      ontology = "CC",
                      count_cutoff = 3)
names(ora.CC_list) <- all_clusters
dotplot(merge_result(ora.CC_list), font.size = 12, title = paste0("ORA - CC - ", TIMEPOINT)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(paste0("cluster_profiler/ORA/", "ORA_CC_", TIMEPOINT, ".png"), width = 8, height = 6)
 
ora.MF_list <- lapply(all_clusters,
                      run_ORA_list,
                      degs.filtered = degs.filtered,
                      seurat_obj = sr,
                      ontology = "MF",
                      count_cutoff = 3)
names(ora.MF_list) <- all_clusters
dotplot(merge_result(ora.MF_list), font.size = 12, title = paste0("ORA - MF - ", TIMEPOINT)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(paste0("cluster_profiler/ORA/", "ORA_MF_", TIMEPOINT, ".png"), width = 8, height = 6)
 
ora.BP <- map_dfr(all_clusters,
                    run_ORA,
                    degs.filtered = degs.filtered,
                    seurat_obj = sr,
                    ontology = "BP")
ora.CC <- map_dfr(all_clusters,
                    run_ORA,
                    degs.filtered = degs.filtered,
                    seurat_obj = sr,
                    ontology = "CC")
ora.MF <- map_dfr(all_clusters,
                    run_ORA,
                    degs.filtered = degs.filtered,
                    seurat_obj = sr,
                    ontology = "MF")
 
worksheets <- list(BP = ora.BP, CC = ora.CC, MF = ora.MF)
write.xlsx(worksheets, paste0("cluster_profiler/ORA/", "ORA_", TIMEPOINT, ".xlsx"), overwrite = T)
 
## GSEA - Gene Set Enrichment Analysis
 
# Create directory
dir.create("cluster_profiler/GSEA", showWarnings = F, recursive = T)
 
# Run GSEA for each cluster separately
 
run_GSEA_list <- function(cluster_, seurat_obj, ontology){
  print(paste0("working on cluster: ", cluster_))
 
  # Subset the Seurat object to keep cluster of interest and only RML and CD1 groups
  srTmp <- subset(seurat_obj, subset = cluster_full_name == cluster_ & inocula %in% c("RML", "CD1"))
 
  # Perform a fast Wilcoxon rank sum test using presto
  gsea.genes <- presto::wilcoxauc(srTmp, group_by = 'inocula')
  gsea.genes <- gsea.genes[which(gsea.genes$group == "RML"),]
 
  geneList <- gsea.genes$logFC
  names(geneList) <- gsea.genes$feature
  geneList <- sort(geneList, decreasing = TRUE)
 
  ego <- tryCatch({
    gseGO(geneList     = geneList,
          OrgDb        = ahOrgDb,
          ont          = ontology,
          keyType       = "SYMBOL",
          minGSSize    = 10,
          maxGSSize    = 500,
          pvalueCutoff = 0.05)
  },
  error = function(cond) NULL)
 
  if(!is.null(ego) && nrow(ego) > 0) {
    return(ego@result)
  }
}
 
gsea.BP <- lapply(all_clusters,
                  run_GSEA_list,
                  seurat_obj = sr,
                  ontology = "BP")
names(gsea.BP) <- all_clusters
 
gsea.MF <- lapply(all_clusters,
                  run_GSEA_list,
                  seurat_obj = sr,
                  ontology = "MF")
names(gsea.MF) <- all_clusters
 
gsea.CC <- lapply(all_clusters,
                  run_GSEA_list,
                  seurat_obj = sr,
                  ontology = "CC")
names(gsea.CC) <- all_clusters
 
# Save the results
worksheets <- list(BP = bind_rows(gsea.BP, .id = "cluster"),
                   CC = bind_rows(gsea.CC, .id = "cluster"),
                   MF = bind_rows(gsea.MF, .id = "cluster"))
write.xlsx(worksheets, paste0("cluster_profiler/GSEA/", "GSEA_", TIMEPOINT, "_per_cluster.xlsx"), overwrite = T)
 
# Run GSEA for all cells of all clusters
gsea.genes <- presto::wilcoxauc(subset(sr, subset = inocula %in% c("RML", "CD1")),
                              group_by = 'inocula')
gsea.genes <- gsea.genes[which(gsea.genes$group == "RML"),]
 
geneList <- gsea.genes$logFC
names(geneList) <- gsea.genes$feature
geneList <- sort(geneList, decreasing = TRUE)
 
run_GSEA_all_clusters <- function(ontology) {
  gsea <- gseGO(geneList = geneList,
                   OrgDb = ahOrgDb,
                   ont = ontology,
                   keyType = "SYMBOL",
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05)
  godata <- GOSemSim::godata('org.Mm.eg.db', ont = ontology)
  gsea <- enrichplot::pairwise_termsim(gsea, method="Wang", semData = godata)
  return(gsea)
}
 
gsea.BP <- run_GSEA_all_clusters(ontology = "BP")
gsea.CC <- run_GSEA_all_clusters(ontology = "CC")
gsea.MF <- run_GSEA_all_clusters(ontology = "MF")
 
# Save results
worksheets <- list(BP = gsea.BP, CC = gsea.CC, MF = gsea.MF)
write.xlsx(worksheets, paste0("cluster_profiler/GSEA/", "GSEA_", TIMEPOINT, "_all_clusters.xlsx"), overwrite = T)
 
# Save gseaResult objects for the generation of plots
saveRDS(gsea.BP, paste0("cluster_profiler/GSEA/", "gseaResult_BP_", TIMEPOINT, ".rds"))
saveRDS(gsea.CC, paste0("cluster_profiler/GSEA/", "gseaResult_CC_", TIMEPOINT, ".rds"))
saveRDS(gsea.MF, paste0("cluster_profiler/GSEA/", "gseaResult_MF_", TIMEPOINT, ".rds"))
 
# Plot
ridgeplot(gsea.BP) +
  labs(x = "enrichment distribution", y = "GO terms") +
  ggtitle(paste0("GSEA - BP - ", TIMEPOINT))
ggsave(paste0("cluster_profiler/GSEA/GSEA_BP_", TIMEPOINT, ".png"), width = 10, height = 10)
 
ridgeplot(gsea.CC) +
  labs(x = "enrichment distribution", y = "GO terms") +
  ggtitle(paste0("GSEA - CC - ", TIMEPOINT))
ggsave(paste0("cluster_profiler/GSEA/GSEA_CC_", TIMEPOINT, ".png"), width = 10, height = 14)
 
ridgeplot(gsea.MF) +
  labs(x = "enrichment distribution", y = "GO terms") +
  ggtitle(paste0("GSEA - MF - ", TIMEPOINT))
ggsave(paste0("cluster_profiler/GSEA/GSEA_MF_", TIMEPOINT, ".png"), width = 10, height = 10)
 
 
# Cleanup
rm(ah, gse.BP, gse.CC, gse.MF, ora.BP, ora.CC, ora.MF,
   worksheets, ahOrgDb, run_GSEA, run_ORA, all_clusters,
   godata_MF, godata_CC, godata_BP)