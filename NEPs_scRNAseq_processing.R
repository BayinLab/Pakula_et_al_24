# Load test data from STARsolo run, process with Seurat

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(stringr)
library(dplyr)
library(ggplot2)
library(gridExtra)

set.seed(123)

## FUNCTIONS
# Read in STARsolo data
read_solo <- function (s) {
    spliced_matrix <- ReadMtx(
      mtx = paste(s, "Solo.out/Velocyto/filtered/spliced.mtx", sep=""), 
          features = paste(s, "Solo.out/Velocyto/filtered/features.tsv", sep=""),
          cells = paste(s, "Solo.out/Velocyto/filtered/barcodes.tsv", sep="")
    )
    unspliced_matrix <- ReadMtx(
      mtx = paste(s, "Solo.out/Velocyto/filtered/unspliced.mtx", sep=""), 
      features = paste(s, "Solo.out/Velocyto/filtered/features.tsv", sep=""),
      cells = paste(s, "Solo.out/Velocyto/filtered/barcodes.tsv", sep="")
    )
    seurat_object <- CreateSeuratObject(counts = spliced_matrix, assay = "spliced", project=s)
    seurat_object[['unspliced']] <- CreateAssayObject(counts = unspliced_matrix)
    
    seurat_object[["RNA"]] <- seurat_object[["spliced"]]
    DefaultAssay(seurat_object) <- "RNA"
    
    return(seurat_object)
}

# list of seurat objects
so_list = list()

# Loop through sample names and make a list of Seurat objects
for (x in c("P1_1", "P1_2", "P2_IR_1", "P2_IR_2", "P2_nonIR_1", "P2_nonIR_2", "P3_IR_1", "P3_IR_2", 
            "P3_nonIR_2", "P5_IR_1", "P5_nonIR_1")) {
  s <- read_solo(x)
  so_list <- append(so_list, s)
}


# Run NormalizeData and FindVariableFeatures for each dataset
ifnb.list <- lapply(X = so_list, FUN = function(x) {
                            x <- NormalizeData(x)
                            x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
                          })

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

# find integration anchors
rna.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

# Integrate the datasets - this command creates an 'integrated' data assay
nep.rna.combined.solo <- IntegrateData(anchorset = rna.anchors)

###############
# ADD METADATA
###############
# Add meta data on treatment - IR vs. nonIR
treatment <- c()

for (i in 1:length(nep.rna.combined.solo@meta.data$orig.ident)) {
  if(grepl('nonIR', nep.rna.combined.solo@meta.data$orig.ident[i])) {
    treatment[i] <- 'nonIR'
  }
  else if (grepl('IR', nep.rna.combined.solo@meta.data$orig.ident[i])){
    treatment[i] <- 'IR'
  }
  else {
    treatment[i] <-'Pre'
  }
}
nep.rna.combined.solo[['treatment']] <- treatment

# Add meta data for timepoint
timepoint <- c()

for (i in 1:length(nep.rna.combined.solo@meta.data$orig.ident)) {
  y <- str_split(nep.rna.combined.solo@meta.data$orig.ident[i], '_')
  timepoint[i] <- y[[1]][1]
}
nep.rna.combined.solo[['timepoint']] <- timepoint


######################
# DIMENSION REDUCTION
######################
DefaultAssay(nep.rna.combined.solo) <- "integrated"

# SCTransform normalisation of UMI counts
nep.rna.combined.solo <- SCTransform(nep.rna.combined.solo)

# PCA
nep.rna.combined.solo <- RunPCA(nep.rna.combined.solo)
ElbowPlot(nep.rna.combined.solo)

# UMAP (n.neighbours = 5 vs. 50 not hugely different)
nep.rna.combined.solo <- RunUMAP(nep.rna.combined.solo, dims = 1:20, n.neighbours = 30)
DimPlot(nep.rna.combined.solo, reduction = "umap", group.by="orig.ident", pt.size=0.01) 
ggsave("nep.rna.combined.solo.sample.png", width=10, height=10)
DimPlot(nep.rna.combined.solo, reduction = "umap", group.by="timepoint", pt.size=0.01) 
ggsave("nep.rna.combined.solo.timepoint.png", width=10, height=10)
DimPlot(nep.rna.combined.solo, reduction = "umap", group.by="treatment", pt.size=0.01) 
ggsave("nep.rna.combined.solo.treatment.png", width=10, height=10)

# Plot markers of interest
FeaturePlot(
  object = nep.rna.combined.solo,
  features = c('Pax2', 'Slc6a11', 'Ascl1', 'Hopx', 'Foxj1', 'Atoh1'),
  pt.size = 0.01,
  max.cutoff = 'q95',
  ncol = 2
) 
ggsave("nep.rna.combined.solo.markers.png", width=10, height=15)

# Calculate the proportion of mitochondrial reads
nep.rna.combined.solo[["percent.mt"]] <- PercentageFeatureSet(nep.rna.combined.solo, pattern = "^mt-", assay="RNA")

# Plot QC metrics 
VlnPlot(nep.rna.combined.solo, pt.size=0.1, group.by = "orig.ident", features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("nep.rna.combined.solo.qc.png", width=10, height=5)

# Filter out unwanted cells
nep.rna.combined.solo <- subset(nep.rna.combined.solo, subset = nFeature_RNA > 1500 & nCount_RNA < 40000 & percent.mt < 5)

###############
# ADJUST CLUSTERING FOR CELL CYCLE
###############

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
m.cc.genes <- lapply(cc.genes, str_to_title)
s.genes <- m.cc.genes$s.genes
g2m.genes <- m.cc.genes$g2m.genes
DefaultAssay(nep.rna.combined.solo) <- 'integrated'

# Check genes associated with PCs for cell cycle (PCs 10 and 11 look guilty)
RunPCA(nep.rna.combined.solo, features = VariableFeatures(nep.rna.combined.solo), ndims.print = 1:20, nfeatures.print = 20)

# Cell cycle genes (n.b. histones are named differently in Gencode)
FeaturePlot(
  object = nep.rna.combined.solo,
  features = c('Hist1h1b', 'Hist1h2ap', 'Ube2c', 'Pcna', "Top2a"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 2
)

# Score each cell for cell cycle phase
nep.rna.combined.solo <- CellCycleScoring(nep.rna.combined.solo, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
# view cell cycle scores and phase assignments
head(nep.rna.combined.solo[[]])

# Plot phase on UMAP
DimPlot(nep.rna.combined.solo, reduction = "umap", group.by = "Phase")
ggsave("nep.rna.combined.solo.ccphase.png", width=10, height=10)

# Regress out the difference between S and G2M scores rather than the whole cell cycle signal 
# to preserve the distinction between quiescent and proliferating cells
nep.rna.combined.solo$CC.Difference <- nep.rna.combined.solo$S.Score - nep.rna.combined.solo$G2M.Score
nep.rna.combined.solo <- SCTransform(nep.rna.combined.solo, vars.to.regress = c("percent.mt", "CC.Difference"))

# Confirm cell cycle genes cause separation by cell cycle
nep.rna.combined.solo <- RunPCA(nep.rna.combined.solo)
nep.rna.combined.solo <- RunUMAP(nep.rna.combined.solo, dims = 1:20, n.neighbors = 30)
DimPlot(nep.rna.combined.solo, reduction="umap", group.by="Phase")
ggsave("nep.rna.combined.solo.ccphase.corrected.png", width=10, height=10)

# Cluster cells
nep.rna.combined.solo <- FindNeighbors(nep.rna.combined.solo, dims = 1:20)
nep.rna.combined.solo <- FindClusters(nep.rna.combined.solo, algorithm=4) # Leiden clustering

#######################
# ANNOTATE CELL CLUSTERS
#######################

# Identify cluster markers
DefaultAssay(nep.rna.combined.solo) <- 'RNA'
nep.rna.combined.solo.markers <- FindAllMarkers(nep.rna.combined.solo, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)

# Get top10 markers for each cluster
nep.rna.combined.solo.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

# Make a heatmap of the markers
DoHeatmap(nep.rna.combined.solo, assay='SCT', features = top10$gene) + 
  FontSize(x.title = 4, y.title = 2, main=5, x.text=4, y.text=4) + 
  NoLegend()
ggsave("nep.rna.combined.solo.marker.heatmap.png", width=10, height=10)

# Assign cluster names
new.cluster.ids <- c("GCP_1", "BgL_NEP_1", "BgL_NEP_2", "GCP_2", "Ascl1_NEP_1",
                     "BgL_NEP_3", "GCP_3", "Ascl1_NEP_2", "Astrocyte", "BgL_NEP_4",
                     "Interneuron", "GCP_4", "Ependymal_!", "GCP_5", "Oligodendrocyte",
                     "Meninges", "Microglia_1", "Ependymal_2", "Ependymal_3", "low_feature_count",
                     "Microglia_2")
names(new.cluster.ids) <- levels(nep.rna.combined.solo)
nep.rna.combined.solo <- RenameIdents(nep.rna.combined.solo, new.cluster.ids)

# View cluster names on UMAP
DimPlot(nep.rna.combined.solo, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("nep.rna.combined.solo.named_clusters.png", width=10, height=10)

# Save data
save.image("nep.rna.combined.solo.annotated.RData")
