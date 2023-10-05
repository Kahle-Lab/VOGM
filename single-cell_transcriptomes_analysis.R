library(Seurat)
library(readxl)
library(tidyverse)
Disease <- readRDS(file='~/Desktop/VOGM/final_merged_vasc_cb.rds')
#VOGM <- readRDS(file='/gpfs/ycga/scratch60/kahle/ga352/VGAM/brain.VOGM.human.counts.rds')
#Meta <- read_excel('/gpfs/ycga/scratch60/kahle/ga352/VGAM/Meta.xlsx')


seurat.obj <- CreateSeuratObject(counts = VOGM)


Meta2 <- Meta[-1]
row.names(Meta2) <- seurat.obj@assays$RNA@counts@Dimnames[[2]]
seurat.obj <- AddMetaData(seurat.obj, Meta2)



seurat.obj <- NormalizeData(seurat.obj)
seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat.obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat.obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(seurat.obj)
seurat.obj <- ScaleData(seurat.obj, features = all.genes)
seurat.obj <- RunPCA(seurat.obj, features = all.genes)
ElbowPlot(seurat.obj)

seurat.obj <- FindNeighbors(seurat.obj, dims = 1:20)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)

seurat.obj <- RunUMAP(seurat.obj, dims = 1:10)

seurat.obj <- SetIdent(seurat.obj, value='celltype')

DimPlot(seurat.obj, reduction = "umap", raster = FALSE)

saveRDS(seurat.obj, file = '/gpfs/ycga/scratch60/kahle/ga352/VGAM/Shujuan_VGAM_VOGM.rds')
seurat <- readRDS(file = '~/Desktop/Shujuan_VGAM.rds')

rm(VOGM, Meta, Meta2, plot1, plot2)

### VOGM!
### Author: Garrett Allington
### Date of Last Update: Friday, December 2nd, 2022
### Built in R Version 4.2.0

### Monocle 3 is still in beta so I recommend a clean install each time to prevent deprecation
devtools::install_github("satijalab/seurat-wrappers", ref="feat/monocle3")

### Open new device
plot.new()
dev.new()

### Load packages
library(Seurat)
library(monocle3)
library(htmlwidgets)
library(DelayedArray)
library(DelayedMatrixStats)
library(BiocGenerics)
library(limma)
library(S4Vectors)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(reticulate)
library(dplyr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(SeuratWrappers)
library(readxl)

### Load pre-compiled seurat object
seurat <- Disease
rm(Disease, monocle_obj)

### Print Clusters
pdf(file = '~/Desktop/Scripts/VOGM/Cluster_Map.pdf', width = 16, height = 9)
DimPlot(seurat, reduction = "umap", raster = FALSE) #Change raster to true to reduce resolution and file size
dev.off()

### Find Markers
VOGM.markers <- FindAllMarkers(seurat)
saveRDS(VOGM.markers, file = "~/Desktop/Scripts/VOGM/VOGM.markers.rds")
#VOGM.markers <-readRDS(file = "~/Desktop/Scripts/VOGM/VOGM.markers.rds")
TopMarks <- VOGM.markers %>%
  group_by(cluster) 
write.table(TopMarks, file = "~/Desktop/Scripts/VOGM/TopMarks.txt")

VOGM_Top20_Markers <- VOGM.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

write.table(VOGM_Top20_Markers, file = "~/Desktop/Scripts/VOGM/VOGM_Top20_Markers.txt")

######################################
#### Convert Seurat object to CDS ####
######################################

### Set up components and parameters for conversion
gene_annotation <- as.data.frame(rownames(seurat@assays[["SCT"]]@meta.features),
                                 row.names = rownames(seurat@assays[["SCT"]]@meta.features))
colnames(gene_annotation) <- "gene_short_name"
cell_metadata <- as.data.frame(seurat@assays[["SCT"]]@counts@Dimnames[[2]],
                               row.names = seurat@assays[["SCT"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
New_matrix <- seurat@assays[["SCT"]]@counts
New_matrix <- New_matrix[rownames(seurat@assays[["SCT"]]@meta.features), ]
expression_matrix <- New_matrix

### Create CDS and recreate seurat clustering exactly as it previously appeared
library(monocle3)
cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)
recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
list_cluster <- seurat@active.ident
names(list_cluster) <- seurat@assays[["SCT"]]@data@Dimnames[[2]]
cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-seurat@reductions[["umap"]]@cell.embeddings
cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings



### Remove excess objects (can easily recreate later if needed). 
### This speeds up the learn graph and prevents memory issues
rm(seurat, cell_metadata, expression_matrix, New_matrix, gene_annotation, recreate.partition, list_cluster)

#################################
#### Monocle3-Based Analysis ####
#################################

### Learn graph, this step usually takes a significant period of time for larger samples (Depending on system and samples, expect well over an hour)
cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = T)

saveRDS(cds_from_seurat, file = "/home/ga352/scratch60/VGAM/VOGM_CDS_after_learn_graph.RDS")
cds_from_seurat <- readRDS(file = "~/Desktop/VOGM_CDS_after_learn_graph.RDS")

### Plot cluster info without trajectory
pdf(file = '~/Desktop/Scripts/VOGM/Graphs/Monocle3_Cluster_Map.pdf', width = 16, height = 9)
plot_cells(cds_from_seurat, 
           color_cells_by = 'cluster',
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph = FALSE)
dev.off()

### Plot cluster info with trajectory
pdf(file = '~/Desktop/Scripts/VOGM/Graphs/Monocle3_Cluster_Map_With_Trajectory.pdf', width = 16, height = 9)
plot_cells(cds_from_seurat, 
           color_cells_by = 'cluster',
           label_groups_by_cluster=TRUE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           show_trajectory_graph = TRUE)
dev.off()

### Before we order the cells in pseudotime, we need to identify which is the principal node/cells
### We can either use their default GUI (which unfortunately is not compatible with 3D), or, if there is a clear progenitor cell type, select them directly
### For the sake of completeness, here is how you select a specific cell type. The example here is neuroprogenitor cells (NPC)
root_cell_list <- grep("NPC", colData(cds_from_seurat))
root_cell_list <- counts(cds_from_seurat)@Dimnames[[2]][root_cell_list]
cds_from_seurat <- order_cells(cds_from_seurat, reduction_method = "UMAP", root_cells = root_cell_list)

### This is using the GUI (inferior to selecting all root cells directly) 
cds_from_seurat <- order_cells(cds_from_seurat, reduction_method = "UMAP")

### Change name to keep record of ordered cells and remove superfluous file. 
cds <- cds_from_seurat
rm(cds_from_seurat) #You may want to save the cds before ordering if you are unsure about your ordering

### Load packages necessary to make pretty graphs :)
library(RColorBrewer)
library(viridis)

### Plot cell clusters without trajectory
pdf(file = "~/Desktop/Scripts/VOGM/Graphs/Monocle3_Cluster_Map.pdf")
plot_cells(cds, color_cells_by="cluster", group_cells_by = "cluster", label_groups_by_cluster = TRUE, 
           graph_label_size = 1.5, group_label_size = 4, trajectory_graph_segment_size = 0.5, show_trajectory_graph = FALSE)
dev.off()

### Plot cell clusters with trajectory
pdf(file = "~/Desktop/Scripts/VOGM/Graphs/Monocle3_Cluster_Map_With_Trajectory.pdf")
plot_cells(cds, color_cells_by="cluster", group_cells_by = "cluster", label_groups_by_cluster = TRUE, 
           graph_label_size = 1.5, group_label_size = 4, trajectory_graph_segment_size = 0.5)
dev.off()

### Plot pseudotime trajectory
pdf(file = "~/Desktop/Scripts/VOGM/Graphs/Monocle3_Clusters_With_Pseudotime.pdf")
plot_cells(cds, color_cells_by="pseudotime", group_cells_by = "cluster", label_groups_by_cluster = TRUE, show_trajectory_graph = FALSE,
           graph_label_size = 1.5, group_label_size = 2, trajectory_graph_segment_size = 0.5)
dev.off()

### Find marker genes expressed by each cluster. Can take some time for larger data sets 
marker_test_res <- top_markers(cds, group_cells_by="cluster")

### Filter to find top 3 markers and create graph
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.05) %>%
  group_by(cell_group) %>%
  top_n(3, marker_score) #Change 3 to find any number of marker genes per cell type desired. Anything above 3-5 gets too busy for most data sets
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

### Make cell cluster marker expression graph
pdf(file = "~/Desktop/Scripts/VOGM/Graphs/Cluster_Top3_Marker_Gene_Epression.pdf", width = 8, height = 14)
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    axis_order = "group_marker",
                    max.size=10)
dev.off()

### Load in gene lists
AVM <- c("KRAS", "NRAS", "BRAF", "MAP2K1", "NOTCH4", "EPHB4", "RASA1", "ALK1", "SMAD4", "ENG", "GNAQ")
CCM <- c("CCM2", "PIK3CA", "MAP3K3", "KRIT1", "PDCD10")
VOGM <- c("RASA1", "EPHB4", "ENG", "ACVRL1", "PTPN11", "NOTCH1")
AVM_VOGM <- c("KRAS", "NRAS", "BRAF", "MAP2K1", "NOTCH4", "EPHB4", "RASA1", "ALK1", "SMAD4", "ENG", "GNAQ", "ACVRL1", "PTPN11", "NOTCH1")
CHD <- c('ABCC9', 	'ABCD3', 	'ACTB', 	'ACVR2B', 	'ADAMTS10', 	'ADNP', 	'ANKRD11', 	'ARHGAP31', 	'ARID1A', 	'ARL13B', 	'ARMC4', 	'ASXL1', 	'ATIC', 	'B3GALT6', 	'BBS1', 	'BBS10', 	'TRIM32', 	'BBS12', 	'BBS2', 	'ARL6', 	'BBS4', 	'BBS5', 	'MKKS', 	'BBS7', 	'TTC8', 	'BBS9', 	'BCOR', 	'BRAF', 	'CACNA1C', 	'CBL', 	'CCDC103', 	'CCDC114', 	'CCDC151', 	'CCDC39', 	'CCDC40', 	'CD96', 	'CDKN1C', 	'CEP290', 	'CEP41', 	'CEP57', 	'CFC1', 	'CHD4', 	'CHD7', 	'CHST3', 	'CITED2', 	'COL1A1', 	'COL1A2', 	'COL2A1', 	'COL3A1', 	'COL5A1', 	'COL5A2', 	'COX7B', 	'CREBBP', 	'CRELD1', 	'DDX11', 	'DGCR2', 	'DHCR7', 	'DNAAF1', 	'DNAAF2', 	'DNAAF3', 	'DNAH11', 	'DNAH5', 	'DNAI1', 	'DNAI2', 	'DNAL1', 	'DOCK6', 	'DYNC2H1', 	'DYX1C1', 	'ECE1', 	'EFTUD2', 	'EHMT1', 	'ELN', 	'EOGT', 	'ESCO2', 	'EVC', 	'EVC2', 	'FBN1', 	'FBN2', 	'FGF8', 	'FGFR1', 	'FIG4', 	'FKTN', 	'FLNA', 	'FOXC1', 	'FOXC2', 	'FTO', 	'GATA4', 	'GATA6', 	'GDF1', 	'GJA1', 	'GLI3', 	'GPC3', 	'GPC6', 	'HCCS', 	'HOXA1', 	'HRAS', 	'IFT122', 	'IFT140', 	'IFT80', 	'INVS', 	'IRX5', 	'JAG1', 	'C5ORF42', 	'KANSL1', 	'KAT6B', 	'KIF7', 	'KMT2A', 	'KMT2D', 	'KRAS', 	'LBR', 	'LEFTY2', 	'LRP2', 	'LTBP4', 	'MAP2K1', 	'MAP2K2', 	'MED12', 	'MED13L', 	'MEGF8', 	'MGP', 	'MID1', 	'MKS1', 	'MYH6', 	'NEK1', 	'NF1', 	'NIPBL', 	'NKX2-5', 	'NKX2-6', 	'NME8', 	'NODAL', 	'NOTCH1', 	'NOTCH2', 	'NPHP3', 	'NPHP4', 	'NEK8', 	'NR2F2', 	'NRAS', 	'NSD1', 	'NSDHL', 	'OFD1', 	'PEX1', 	'PEX13', 	'PHGDH', 	'PITX2', 	'PKD1', 	'PKD2', 	'PLOD1', 	'PKD1L1', 	'PQBP1', 	'PTEN', 	'PTPN11', 	'RAB23', 	'RAD21', 	'RAF1', 	'RAI1', 	'RBM10', 	'RBM8A', 	'RIT1', 	'RNU4ATAC', 	'ROR2', 	'RPGRIP1L', 	'RPL11', 	'RPL35A', 	'RPL5', 	'RPS10', 	'RPS17', 	'RPS19', 	'RPS24', 	'RPS26', 	'RPS7', 	'RSPH4A', 	'RSPH9', 	'SALL1', 	'SEMA3E', 	'SETBP1', 	'SF3B4', 	'SH3PXD2B', 	'SHH', 	'SHOC2', 	'SKI', 	'SMAD3', 	'SMAD4', 	'SMARCA4', 	'SMARCB1', 	'SMARCE1', 	'SMC3', 	'SMS', 	'SOS1', 	'SOX2', 	'SOX9', 	'STAMBP', 	'STRA6', 	'TAB2', 	'TBX1', 	'TBX20', 	'TBX3', 	'TBX5', 	'TCOF1', 	'TFAP2B', 	'TGFBR1', 	'TGFBR2', 	'TLL1', 	'TSC1', 	'TSC2', 	'TTC21B', 	'TWIST1', 	'UBR1', 	'WDR19', 	'WDR35', 	'WDR60', 	'ZEB2', 	'ZFPM2', 	'ZIC3', 	'ACAN', 	'ADAMTS6', 	'ANKS6', 	'AP1B1', 	'AP2B1', 	'BICC1', 	'CC2D2A', 	'CNTRL', 	'NAT8', 	'CXCR4', 	'DAW1', 	'DCTN5', 	'DNM2', 	'DRC1', 	'FOXJ1', 	'FREM2', 	'FUZ', 	'C1ORF127', 	'HECTD1', 	'IFT74', 	'LOX', 	'LRP1', 	'LTBP1', 	'MMP21', 	'MYH10', 	'NDST1', 	'PCSK5', 	'PDE2A', 	'PLXND1', 	'PRDM1', 	'PRICKLE1', 	'PSKH1', 	'PTK7', 	'ROBO1', 	'SMAD6', 	'SNX17', 	'SUFU', 	'TAB1', 	'TBC1D32', 	'TMEM67', 	'ZBTB14', 	'FLT4', 	'KDR', 	'IQGAP1', 	'SMAD2', 	'KDM5B', 	'FOXP1','DYRK1A')
Height <- c('ACAN',	'ADAM28',	'ADAMTS10',	'ADAMTS17',	'ADAMTS3',	'ADAMTSL3',	'AK092571',	'AKD1',	'AL117656',	'AL161980',	'ANKRD13B',	'ATP6V1E2',	'B3GNT8',	'BC030091',	'BMP2',	'BMP6',	'BNC2',	'BRCA1',	'C17orf42',	'C17orf82',	'C20orf199',	'C3orf63',	'C6orf1',	'C6orf173',	'CCDC100',	'CCDC108',	'CCDC66',	'CCDC91',	'CDK2AP1',	'CDK6',	'CLIC4',	'CLPS',	'CPN1',	'CWF19L1',	'CYP19A1',	'DDX27',	'DNM3',	'DTL',	'DYM',	'ECM2',	'EFEMP1',	'EIF2AK3',	'ESR1',	'ETV6',	'EXOSC5',	'FAM173A',	'FANCE',	'FARP2',	'FBLN2',	'FBXW11',	'FGFR4',	'FLI1',	'FNDC3B',	'FOLH1',	'FRS2',	'FUBP3',	'GDF5',	'GH1',	'GHSR',	'GNA12',	'GNPTAB',	'GPC5',	'H1FX',	'HAGHL',	'HDLBP',	'HHIP',	'HLA-__B',	'HMGA1',	'HMGA2',	'HSS00017874',	'HSS00085450',	'HSS00174467',	'ID4',	'IGF1R',	'IGF2BP2',	'IGF2BP3',	'IHH',	'INSR',	'INTS7',	'ITPR3',	'JMJD4',	'KCNQ1',	'L3MBTL3',	'LIN28',	'LIN28B',	'LPAR1',	'LRRC37B',	'LTBP1',	'LTBP2',	'LUZP1',	'MC4R',	'MEF2C',	'MFAP2',	'MICA',	'MKL2',	'MSTP9',	'MTMR11',	'MYO9B',	'N4BP2L2',	'NFATC4',	'NFIC',	'NME2',	'NOG',	'NPPC',	'NPR3',	'NSD1',	'NUCB2',	'PAPPA',	'PAPPA2',	'PCSK5',	'PDS5B',	'PEX2',	'PFAAP5',	'PIP4K2B',	'PIP5K2B',	'PITX1',	'PKN2',	'PML',	'PPA2',	'PPAP2A',	'PPARD',	'PPIF',	'PRKCZ',	'PRKG2',	'PSMB3',	'PTCH1',	'PTPRJ',	'QSCN6L1',	'QSOX2',	'REST',	'RMI1',	'RNF135',	'RPL5',	'RUNX2',	'RYBP',	'SCMH1',	'SLBP',	'SLC22A4',	'SLC22A5',	'SLC23A3',	'SLC39A13',	'SLIT3',	'SMPD2',	'SNAP47',	'SOCS2',	'SOCS5',	'STAT2',	'STK36',	'TBX2',	'TCF19',	'TEAD1',	'TGFB2',	'TMEM4',	'TNS1',	'TP53I13',	'TRIP11',	'TSEN15',	'USP52',	'UTP6',	'VGLL2',	'VPS13C',	'ZBTB24',	'ZNF142',	'ZNF311',	'ZNXF1')
pVOGM <- c('NOTCH1', 'RASA1', 'PTPN11', 'MUC5B', 'KMT2D', 'RAB11FIP4', 'SMARCA1', 'KAT6A', 'ANK2', 'LPGAT1', 'EPHB4', 'ENG', 'ACVRL1')
MMD <- c('RNF213', 'GUCY1A3', 'BRCC3', 'CMC4', 'MTPC1', 'DIAPH1', 'CNOT3', 'ACTA2', 'CHD4', 'SETD5', 'ACTA2', 'PCNT', 'YY1AP1', 'JAG1', 'SMARCAL1')
MasterList <- c(AVM, CCM, VOGM, pVOGM,MMD)
MasterList <- MasterList[!duplicated(MasterList)]

### Plot by-cluster enrichment for associated genes
pdf(file="~/Desktop/Scripts/VOGM/Graphs/Cluster_CCM_Gene_Expression.pdf", height = 5, width = 5)
plot_genes_by_group(cds,
                    CCM,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    axis_order = "group_marker",
                    max.size=10)
dev.off()

### Plot gene expression in clusters
pdf(file="~/Desktop/Scripts/VOGM/Graphs/Gene_Cluster_Map.pdf", height = 18, width = 32)
plot_cells(cds,
           genes=CCM,
           color_cells_by = "cluster",
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

##########################################
#### Differential Expression Analysis ####
##########################################

### Subset CDS by gene lists
cds_subset_ML <- cds[rowData(cds)$gene_short_name %in% CCM,]

### Call fit_models() to fit a generalized linear model for each gene in a cell_data_set
gene_fits <- fit_models(cds_subset_ML, model_formula_str = "~cluster")

### Print fit coefficient table
fit_coefs <- coefficient_table(gene_fits)

### Filter for q_value < 0.05
fit_coefs_filtered <- fit_coefs
fit_coefs_filtered %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

### Create violin plots displaying expression
pdf(file="~/Desktop/Scripts/VOGM/Graphs/VOGM_Gene_Expression_Violin.pdf", height = 8, width = 8)
plot_genes_violin(cds_subset_ML) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

##########################
#### Modular Analysis ####
##########################

### Load requisite packages
library(coda)
library(LearnBayes)
library(gmodels)
library(expm)

### Run Graph_Test -- This is very computationally-intensive (~8hr run time)
pr_graph_test_res <- graph_test(cds, neighbor_graph="knn", reduction_method = "UMAP",
                                k = 25, alternative = "greater", expression_family = "quasipoisson")

saveRDS(pr_graph_test_res, file = "~/Desktop/Scripts/VOGM/Files/pr_graph_test_res.RDS")
#pr_graph_test_res <-readRDS(file = "~/Desktop/Scripts/VOGM/Files/pr_graph_test_res.RDS")

### Subset by q_value
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

### Preprocess CDS
test <- seurat@assays$RNA@counts@Dimnames[[1]]
test <- as.data.frame(test)
cds@preprocess_aux$gene_loadings <- test
rownames(cds@preprocess_aux$gene_loadings) <- cds@preprocess_aux$gene_loadings[,1]
CDSP <- preprocess_cds(cds, method = 'PCA')

### Build gene modules
gene_module_df <- find_gene_modules(CDSP[pr_deg_ids,],
                                    reduction_method = "UMAP",
                                    max_components = 2,
                                    umap.metric = "cosine",
                                    umap.min_dist = 0.1,
                                    umap.n_neighbors = 15L,
                                    umap.fast_sgd = FALSE,
                                    umap.nn_method = "annoy",
                                    k = 20,
                                    leiden_iter = 1,
                                    partition_qval = 0.05,
                                    weight = FALSE,
                                    resolution = .01,
                                    random_seed = 0L,
                                    cores = 1,
                                    verbose = T)

### Export gene module data frame and analyze modular gene ontology profiles
write.table(gene_module_df, file = "~/Desktop/Scripts/VOGM/Files/gene_modules.txt")
#gene_module_df <- read_excel("~/Desktop/Scripts/VOGM/Files/gene_modules.xlsx", sheet = "gene_modules")

### Continuously recreate modules until you are happy with module resolution and specificity. 
### This will take many iterations

### Reshape gene_module_df so each module is a column
GMods <- gene_module_df %>% pivot_wider(names_from =module, values_from = id)
GMods <- GMods[, -c(1:3)]
GMods %>% select(order(colnames(GMods)))
original_cols <- colnames(GMods)
colnames(GMods) <- paste0("Module_" ,original_cols)
for (col in colnames(GMods)) {
  assign(col, GMods[[col]]) # create new object with column name
}
GMod_Cols <- as.vector(colnames(GMods))

# loop through each character vector
for (col in GMod_Cols) {
  # remove NA values
  assign(col, na.omit(get(col)))
}

# determine length of longest vector
max_length <- max(sapply(mget(GMod_Cols), length))

# loop through each character vector
for (col in GMod_Cols) {
  # get current vector
  current_vec <- get(col)
  # add NA values to end of vector until it is the same length as the longest vector
  current_vec <- c(current_vec, rep(NA, max_length - length(current_vec)))
  # update object in working environment with new vector
  assign(col, current_vec)
}

### Make table
# create empty list to store character vectors
char_vec_list <- list()

# loop through each module number
for (module_num in 1:23) {
  # get name of current module
  current_name <- paste0("Module_", module_num)
  # get current module vector and add it to list
  current_vec <- get(current_name)
  char_vec_list[[module_num]] <- current_vec
}

# convert list of character vectors to dataframe
GMods <- data.frame(char_vec_list)

# rename columns of dataframe
colnames(GMods) <- paste0("Module_", 1:23)


# Save table
write.table(GMods, file = "~/Desktop/Scripts/VOGM/Files/GMods.txt")
GMods <- read_excel("~/Desktop/Scripts/VOGM/Files/GMods.xlsx") 

### Plot cell-gene module maps
pdf(file="~/Desktop/Scripts/VOGM/Graphs/UMAP_Gene_Module_Expression.pdf", height = 20, width = 20)
plot_cells(CDSP, 
           genes=gene_module_df,
           group_cells_by="cluster",
           color_cells_by="cluster",
           show_trajectory_graph=FALSE)
dev.off()

### Create cell group dataframe for clusters 
cell_group_df <- tibble::tibble(cell=row.names(colData(CDSP)), 
                                cell_group=CDSP@clusters@listData$UMAP$clusters)

### Create aggregate gene expression matrix
agg_mat <- aggregate_gene_expression(CDSP, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

#Making agg_mat_NS from agg_mat
mods <- rep(0,nrow(agg_mat_df)) #Make list of module titles
for(i in 1:nrow(agg_mat_df)){
  mods[i] <- (paste0('Module_',i))
}



#colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
agg_mat_mtx <- as.matrix(agg_mat)
agg_mat_df <- as.data.frame(agg_mat_mtx)
write.table(agg_mat_df, file = "~/Desktop/Scripts/VOGM/Files/agg_mat.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
#agg_mat_df <- read_excel("~/Desktop/Scripts/VOGM/Files/agg_mat.xlsx", sheet = "agg_mat") 
#agg_mat_df <- as.data.frame(agg_mat_df) #Only needed if re-uploading from excel


### Create cluster and gene module heatmap
pdf(file = "~/Desktop/Scripts/VOGM/Graphs/Module_Cluster_Heatmap.pdf", width = 8, height = 10)
pheatmap::pheatmap(agg_mat_df, cluster_rows=FALSE, cluster_cols=FALSE,
                   scale="column", #clustering_method="ward.D2",
                   border_color = NA,
                   color = magma(50),
                   cellwidth = 50, 
                   cellheight = 20, 
                   main = "Cluster Gene Module Enrichment", 
                   #color = colorRampPalette(c("blue", "seagreen4", "gold"))(50),
                   fontsize=12)
dev.off()

### Create cluster and gene module heatmap with normalization and significance
library(readxl)
agg_mat_NS <- read_excel("~/Desktop/Scripts/VOGM/Files/agg_mat.xlsx", sheet = "Sheet1", col_names = TRUE)
agg_mat_NS$Module <- factor(agg_mat_NS$Module, levels=unique(agg_mat_NS$Module))
agg_mat_NS$Module <- factor(agg_mat_NS$Module, levels=c("Module_1", "Module_2", "Module_3", "Module_4", "Module_5", "Module_6",
                                                        "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12",
                                                        "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", 
                                                        "Module_19", "Module_20", "Module_21", "Module_22", "Module_23"))

pdf(file = '~/Desktop/Scripts/VOGM/Graphs/Modular_CellType_Enrichment.pdf', height=9, width=16)
ggplot(agg_mat_NS, aes(x=Module, y= Cell)) +
  geom_tile(aes(fill = Score), color = "white") +
  scale_fill_gradient2(low = 'black', mid = 'maroon', high = 'papayawhip') +
  theme(
    plot.title = element_text(face="bold", hjust = 0.5),
    axis.line        = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    axis.text.y      = element_text(angle = 0 , size = 12 ), 
    axis.title.x     = element_text(angle = 0, size = 12, face = "bold"),
    axis.title.y     = element_text(angle = 90, size = 12, face = "bold"),
    legend.position  = "top", 
    legend.title     = element_text(size=12, face = "bold"),
    legend.text      = element_text(size=10) 
  ) +
  geom_text(aes(label=ifelse(as.numeric(Score)>1.3, round(Score, digits=1), '')), color = "black") + 
  ggtitle("Modular Cell-Type Enrichment") +
  guides(fill=guide_legend(title="-log10(Adjusted P-Value)"))
dev.off()

############################################
#### Hypergeometric Enrichment Analysis ####
############################################

### Set up iterative analysis
library(readxl)
GMods <- read_excel("~/Desktop/Scripts/VOGM/Files/gene_modules.xlsx", sheet = "Sheet1") 
ModColSums <- colSums(!is.na(GMods))
gene_number <- 26812 #Change to number of transcripts observed in scRNA data

CONDITIONS.GENES <- NULL
CONDITIONS.GENES[['AVM']] <- AVM 
CONDITIONS.GENES[['AVM_VOGM']] <- AVM_VOGM
CONDITIONS.GENES[['CCM']] <- CCM
CONDITIONS.GENES[['MasterList']] <- MasterList
CONDITIONS.GENES[['CHD']] <- CHD
CONDITIONS.GENES[['pVOGM']] <- pVOGM
CONDITIONS.GENES[['VOGM']] <- pVOGM
CONDITIONS.GENES[['MMD']] <- MMD
CONDITIONS.GENES[['Height']] <- Height

MME <- NULL
ModColSums <- colSums(!is.na(GMods))
for (disease in 1:length(CONDITIONS.GENES) ){
  pVals <- data.frame(var1=c(4))
  for(i in 1:ncol(GMods)){
    pVals[nrow(pVals) + 1,] = c(phyper(length(intersect(CONDITIONS.GENES[[disease]], 
                                                        GMods[[i]])),length(CONDITIONS.GENES[[disease]]),
                                       (gene_number-length(CONDITIONS.GENES[[disease]])), ModColSums[i], lower.tail=FALSE))
  }
  pVals = pVals[-1, ,drop = FALSE]
  pVals$Module <- names(GMods)
  pVals$pVals <- pVals$var1
  pVals$LOG <- -log10(pVals$pVals)
  pVals$Set <- names(CONDITIONS.GENES[disease])
  pVals <- pVals[,-1]
  MME[[disease]] <- pVals
}
MME <- bind_rows(MME) #Collapse into DF

### AVM
AVM_pVals <- data.frame(var1=c(4))
for(i in 1:ncol(GMods)){
  AVM_pVals[nrow(AVM_pVals) + 1,] = c(phyper(length(intersect(AVM, GMods[[i]])),length(AVM),(gene_number-length(AVM)), ModColSums[i], lower.tail=FALSE))
}
AVM_pVals = AVM_pVals[-1,]
write.table(AVM_pVals, "~/Desktop/Scripts/VOGM/Files/AVM_pVals.txt", sep = " ")

### CCM
CCM_pVals <- data.frame(var1=c(4))
for(i in 1:ncol(GMods)){
  CCM_pVals[nrow(CCM_pVals) + 1,] = c(phyper(length(intersect(CCM, GMods[[i]])),length(CCM),(gene_number-length(CCM)), ModColSums[i], lower.tail=FALSE))
}
CCM_pVals = CCM_pVals[-1,]
write.table(CCM_pVals, "~/Desktop/Scripts/VOGM/Files/CCM_pVals.txt", sep = " ")

### VOGM
VOGM_pVals <- data.frame(var1=c(4))
for(i in 1:ncol(GMods)){
  VOGM_pVals[nrow(VOGM_pVals) + 1,] = c(phyper(length(intersect(VOGM, GMods[[i]])),length(VOGM),(gene_number-length(VOGM)), ModColSums[i], lower.tail=FALSE))
}
VOGM_pVals = VOGM_pVals[-1,]
write.table(VOGM_pVals, "~/Desktop/Scripts/VOGM/Files/VOGM_pVals.txt", sep = " ")

### pVOGM
pVOGM_pVals <- data.frame(var1=c(4))
for(i in 1:ncol(GMods)){
  pVOGM_pVals[nrow(pVOGM_pVals) + 1,] = c(phyper(length(intersect(pVOGM, GMods[[i]])),length(pVOGM),(gene_number-length(pVOGM)), ModColSums[i], lower.tail=FALSE))
}
pVOGM_pVals = pVOGM_pVals[-1,]
write.table(pVOGM_pVals, "~/Desktop/Scripts/VOGM/Files/pVOGM_pVals.txt", sep = " ")

### AVM_VOGM
AVM_VOGM_pVals <- data.frame(var1=c(4))
for(i in 1:ncol(GMods)){
  AVM_VOGM_pVals[nrow(AVM_VOGM_pVals) + 1,] = c(phyper(length(intersect(AVM_VOGM, GMods[[i]])),length(AVM_VOGM),(gene_number-length(AVM_VOGM)), ModColSums[i], lower.tail=FALSE))
}
AVM_VOGM_pVals = AVM_VOGM_pVals[-1,]
write.table(AVM_VOGM_pVals, "~/Desktop/Scripts/VOGM/Files/AVM_VOGM_pVals.txt", sep = " ")

### MMD
MMD_pVals <- data.frame(var1=c(4))
for(i in 1:ncol(GMods)){
  MMD_pVals[nrow(MMD_pVals) + 1,] = c(phyper(length(intersect(MMD, GMods[[i]])),length(MMD),(gene_number-length(MMD)), ModColSums[i], lower.tail=FALSE))
}
MMD_pVals = MMD_pVals[-1,]
write.table(MMD_pVals, "~/Desktop/Scripts/VOGM/Files/MMD_pVals.txt", sep = " ")

### CHD
CHD_pVals <- data.frame(var1=c(4))
for(i in 1:ncol(GMods)){
  CHD_pVals[nrow(CHD_pVals) + 1,] = c(phyper(length(intersect(CHD, GMods[[i]])),length(CHD),(25004-length(CHD)), ModColSums[i], lower.tail=FALSE))
}
CHD_pVals = CHD_pVals[-1,]
write.table(CHD_pVals, "~/Desktop/Scripts/VOGM/Files/CHD_pVals.txt", sep = " ")

### Height
Height_pVals <- data.frame(var1=c(4))
for(i in 1:ncol(GMods)){
  Height_pVals[nrow(Height_pVals) + 1,] = c(phyper(length(intersect(Height, GMods[[i]])),length(Height),(25004-length(Height)), ModColSums[i], lower.tail=FALSE))
}
Height_pVals = Height_pVals[-1,]
write.table(Height_pVals, "~/Desktop/Scripts/VOGM/Files/Height_pVals.txt", sep = " ")

### Import the concatenated MME table
library(readxl)
MME <- read_excel( "~/Desktop/Scripts/VOGM/Files/MME.xlsx")
MME$Module <- as.character(MME$Module)
MME$Module <- factor(MME$Module, levels=unique(MME$Module))
MME$Module <- factor(MME$Module, levels=c("Module_1", "Module_2", "Module_3", "Module_4", "Module_5", "Module_6",
                                          "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12",
                                          "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", 
                                          "Module_19", "Module_20", "Module_21", "Module_22", "Module_23"))
MME$Set <- factor(MME$Set, levels=c("AVM", "VOGM", "AVM_VOGM", "pVOGM", "CCM", "MasterList", "MMD", "CHD", "Height"))

### Create hypergeometric enrichment heatmap
library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape)

pdf(file = '~/Desktop/Scripts/VOGM/Graphs/Modular_Enrichment.pdf', height=9, width=16)
ggplot(MME, aes(x=Module, y= Set)) +
  geom_tile(aes(fill = LOG), color = "white") +
  scale_fill_gradient(low = "aliceblue", high = "navy") +
  theme(
    plot.title = element_text(face="bold", hjust = 0.5),
    axis.line        = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    axis.text.y      = element_text(angle = 0 , size = 12 ), 
    axis.title.x     = element_text(angle = 0, size = 12, face = "bold"),
    axis.title.y     = element_text(angle = 90, size = 12, face = "bold"),
    legend.position  = "top", 
    legend.title     = element_text(size=12, face = "bold"),
    legend.text      = element_text(size=10) 
  ) +
  geom_text(aes(label=ifelse(as.numeric(LOG)>2.4, round(LOG, digits=1), '')), color = "white") + 
  ggtitle("Module Enrichment") +
  guides(fill=guide_legend(title="-log10(p adj)"))
dev.off()

##################
#### EnrichR #####
##################

### This part of the script will query EnrichR for your gene modules and save results as tables
### Tables will then automatically be created to bar graphs for GO and WP terms 
### Other databases (e.g. transcription factors) and functionality (e.g. manhattan plots) are possible through EnrichR, but not included in this script.
### EnrichR requires R version 4.0 or higher

### Load requisite packages
library(rjson)
library(enrichR)
library(readxl)
library(forcats)

### Load in modules
GeneModules <- GMods
### Check cxn, set site, and verify databases 
listEnrichrSites()
setEnrichrSite("Enrichr") #Human or mouse gene input
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)


### Set EnrichR analysis databases to query
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", "WikiPathway_2021_Human")

### Separate each Gene Module into separate dataframes
for(i in 1:ncol(GeneModules)){
  tmp <- subset(GeneModules[i])
  tmp <- tmp[!is.na(tmp)]
  tmp <- as.data.frame(tmp)
  write.table(tmp,
              file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('Module',i,'.csv')),
              sep=',',
              row.names=FALSE)
}

### Run each module df through EnrichR
for(i in 1:ncol(GeneModules)){
  tmp <- read.csv(file=file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('Module',i,'.csv')))
  enriched <- enrichr(tmp[[1]], dbs)
  write.table(enriched$GO_Biological_Process_2021,
              file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('BP_Module',i,'.csv')),
              sep=',',
              row.names=FALSE)
  write.table(enriched$GO_Cellular_Component_2021,
              file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('CC_Module',i,'.csv')),
              sep=',',
              row.names=FALSE)
  write.table(enriched$GO_Molecular_Function_2021,
              file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('MF_Module',i,'.csv')),
              sep=',',
              row.names=FALSE)
  write.table(enriched$WikiPathway_2021_Human,
              file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('WP_Module',i,'.csv')),
              sep=',',
              row.names=FALSE)
}

### Calculate negative Log10 for BP
for (i in 1:ncol(GeneModules)) {
  tmp <- read_csv(file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('BP_Module',i,'.csv')))
  tmp$LOG <- with(tmp, (-log10(tmp$Adjusted.P.value)))
  tmp <- as_tibble(tmp)
  tmp %>% arrange(desc(LOG))
  tmp <- as.data.frame(tmp)
  assign(paste0('tmp',i), tmp)
  write.table(tmp,
              file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('BP_Module',i,'.csv')),
              sep=',',
              row.names=FALSE)
}

### Calculate negative Log10 for CC
for (i in 1:ncol(GeneModules)) {
  tmp <- read_csv(file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('CC_Module',i,'.csv')))
  tmp$LOG <- with(tmp, (-log10(tmp$Adjusted.P.value)))
  tmp <- as_tibble(tmp)
  tmp %>% arrange(desc(LOG))
  tmp <- as.data.frame(tmp)
  assign(paste0('tmp',i), tmp)
  write.table(tmp,
              file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('CC_Module',i,'.csv')),
              sep=',',
              row.names=FALSE)
}

### Calculate negative Log10 for MF
for (i in 1:ncol(GeneModules)) {
  tmp <- read_csv(file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('MF_Module',i,'.csv')))
  tmp$LOG <- with(tmp, (-log10(tmp$Adjusted.P.value)))
  tmp <- as_tibble(tmp)
  tmp %>% arrange(desc(LOG))
  tmp <- as.data.frame(tmp)
  assign(paste0('tmp',i), tmp)
  write.table(tmp,
              file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('MF_Module',i,'.csv')),
              sep=',',
              row.names=FALSE)
}

### Calculate negative Log10 for WP
for (i in 1:ncol(GeneModules)) {
  tmp <- read_csv(file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('WP_Module',i,'.csv')))
  tmp$LOG <- with(tmp, (-log10(tmp$Adjusted.P.value)))
  tmp <- as_tibble(tmp)
  tmp %>% arrange(desc(LOG))
  tmp <- as.data.frame(tmp)
  assign(paste0('tmp',i), tmp)
  write.table(tmp,
              file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('WP_Module',i,'.csv')),
              sep=',',
              row.names=FALSE)
}

### Make Graphs for ALL Modules 
#install.packages("egg")
library(egg)
library(tidyverse)
for (i in 1:ncol(GeneModules)) {
  BP <- read_csv(file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('BP_Module',i,'.csv')))
  CC <- read_csv(file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('CC_Module',i,'.csv')))
  MF <- read_csv(file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('MF_Module',i,'.csv')))
  WP <- read_csv(file= file.path('~/Desktop/Scripts/VOGM/EnrichR/', paste0('WP_Module',i,'.csv')))
  BPO <- BP[order(-BP$LOG),][1:10,]
  CCO <- CC[order(-CC$LOG),][1:10,]
  MFO <- MF[order(-MF$LOG),][1:10,]
  WPO <- WP[order(-WP$LOG),][1:10,]
  BPG <- ggplot(BPO, aes(reorder(Term, LOG), LOG)) +
    geom_bar(stat="identity", fill="navy", alpha=.95, width=0.8) +
    scale_y_continuous(expand = c(0,0.05)) + 
    scale_x_discrete(expand = c(0,0)) +
    coord_flip() +
    ylab("-Log(Adj. P-Value)") +
    xlab("") +
    theme_bw() +
    geom_hline(yintercept = 1.28, color = "green", size=.3) +
    theme(
      axis.text.x=element_text(),
      axis.text.y=element_text(face = "bold", size = 14),
      plot.title = element_text(color="Black", size=14, face="bold", hjust=0.5)) +
    ggtitle(paste0("2021 GO Biological Processes of Module ", i))
  CCG <- ggplot(CCO, aes(reorder(Term, LOG), LOG)) +
    geom_bar(stat="identity", fill="navy", alpha=.95, width=0.8) +
    scale_y_continuous(expand = c(0,0.05)) + 
    scale_x_discrete(expand = c(0,0)) +
    coord_flip() +
    ylab("-Log(Adj. P-Value)") +
    xlab("") +
    theme_bw() +
    geom_hline(yintercept = 1.28, color = "green", size=.3) +
    theme(
      axis.text.x=element_text(),
      axis.text.y=element_text(face = "bold", size = 14),
      plot.title = element_text(color="Black", size=14, face="bold", hjust=0.5)) +
    ggtitle(paste0("2021 GO Cellular Components of Module ", i))
  MFG <- ggplot(MFO, aes(reorder(Term, LOG), LOG)) +
    geom_bar(stat="identity", fill="navy", alpha=.95, width=0.8) +
    scale_y_continuous(expand = c(0,0.05)) + 
    scale_x_discrete(expand = c(0,0)) +
    coord_flip() +
    ylab("-Log(Adj. P-Value)") +
    xlab("") +
    theme_bw() +
    geom_hline(yintercept = 1.28, color = "green", size=.3) +
    theme(
      axis.text.x=element_text(),
      axis.text.y=element_text(face = "bold", size = 14),
      plot.title = element_text(color="Black", size=14, face="bold", hjust=0.5)) +
    ggtitle(paste0("2021 GO Molecular Functions of Module ", i))
  WPG <- ggplot(WPO, aes(reorder(Term, LOG), LOG)) +
    geom_bar(stat="identity", fill="navy", alpha=.95, width=0.8) +
    scale_y_continuous(expand = c(0,0.05)) + 
    scale_x_discrete(expand = c(0,0)) +
    coord_flip() +
    ylab("-Log(Adj. P-Value)") +
    xlab("") +
    theme_bw() +
    geom_hline(yintercept = 1.28, color = "green", size=.3) +
    theme(
      axis.text.x=element_text(),
      axis.text.y=element_text(face = "bold", size = 14),
      plot.title = element_text(color="Black", size=14, face="bold", hjust=0.5)) +
    ggtitle(paste0("2021 GO WikiPathways of Module ", i))
  pdf(file = file.path('~/Desktop/Scripts/VOGM/EnrichR/Graphs/', paste0('Module_', i, '.pdf')), width = 30, height = 14)
  ggarrange(BPG, CCG, MFG, WPG, ncol=2)
  dev.off()
}


