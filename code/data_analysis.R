library(Seurat)
library(dplyr)
library(ggsci)
library(gridExtra)
library(tibble)
source("~/GSE139999-data-analysis/preprocessing.R")

# -------------------------------------------------------
# Read data
# -------------------------------------------------------
gem_female <- prepareData_10X(
    filepath = "path/to/GSM4151331/10X/matrices/",  
    pct.Cells = 0.2,
    minNumGenes = 300,
    maxNumGenes = 5000,
    minNumHk = 65,
    maxPct.Mito = 2
)
gem_male <- prepareData_10X(
    filepath = "path/to/GSM4151330/10X/matrices/",
    pct.Cells = 0.2,
    minNumGenes = 300,
    maxNumGenes = 6000,
    minNumHk = 70,
    maxPct.Mito = 2
)


sobj1 <- CreateSeuratObject(
    counts = gem_female, 
    project = "female", 
    min.cells = round(0.01 * ncol(gem_female))
)
sobj1$stim <- "female"
sobj1 <- subset(sobj1, subset = nFeature_RNA > 500) %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

sobj2 <- CreateSeuratObject(
    counts = gem_male, 
    project = "male", 
    min.cells = round(0.01 * ncol(gem_male))
)
sobj2$stim <- "male"
sobj2 <- subset(sobj2, subset = nFeature_RNA > 500) %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

# -------------------------------------------------------
# Integrate
# -------------------------------------------------------
anchors <- FindIntegrationAnchors(
  object.list = list(sobj1, sobj2), 
  dims = 1:20
)
sobj.combined <- IntegrateData(anchorset = anchors, dims = 1:20)

# -------------------------------------------------------
# Perform integrated analysis
# -------------------------------------------------------
DefaultAssay(sobj.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
sobj.combined <- ScaleData(
  sobj.combined, 
  verbose = FALSE, 
  features = rownames(sobj.combined)
) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

sobj.combined <- subset(sobj.combined, idents = c(0,1,2,3,4,5,6,7,8,9,10))

sobj.combined <- BuildClusterTree(sobj.combined)
PlotClusterTree(sobj.combined)

# Merge clusters 0+1+2+5 and 3+10
Idents(sobj.combined)[WhichCells(sobj.combined, idents = 0)] <- 1
Idents(sobj.combined)[WhichCells(sobj.combined, idents = 2)] <- 1
Idents(sobj.combined)[WhichCells(sobj.combined, idents = 5)] <- 1
Idents(sobj.combined)[WhichCells(sobj.combined, idents = 10)] <- 3
Idents(sobj.combined)[WhichCells(sobj.combined, idents = 9)] <- 3

# Rename idents
new.cluster.ids <- c(5, 4, 3, 6, 2, 1)
names(new.cluster.ids) <- levels(sobj.combined)
sobj.combined <- RenameIdents(sobj.combined, new.cluster.ids)
Idents(sobj.combined) <- factor(Idents(sobj.combined), levels = c(1,2,3,4,5,6))

p1 <- DimPlot(sobj.combined, label = T) + NoAxes() + NoLegend()
p2 <- DimPlot(sobj.combined, group.by = "stim", split.by = 'stim') + NoAxes() + NoLegend()
gridExtra::grid.arrange(p1, p2)

DimPlot(sobj.combined, label = T) + NoAxes() + NoLegend()
DimPlot(sobj.combined, group.by = "stim", split.by = 'stim') + NoAxes() + NoLegend()

sobj.combined <- BuildClusterTree(sobj.combined)
PlotClusterTree(sobj.combined)

# -------------------------------------------------------
# Cluster Frequency per sample
# -------------------------------------------------------
stats <- data.frame(
  ClusterID = unique(Idents(sobj.combined)), 
  pct_Female = NA, 
  pct_Male = NA
)
for (i in 1:nrow(stats)) {
  ID <- stats[i, "ClusterID"]
  cells <- WhichCells(sobj.combined, idents = ID)
  stats[i, "pct_Female"] <- round(length(grep("_1", cells))*100/length(grep("_1", colnames(sobj.combined))), 2)
  stats[i, "pct_Male"] <- round(length(grep("_2", cells))*100/length(grep("_2", colnames(sobj.combined))), 2)
}
df.m <- reshape::melt(stats)

ggplot(df.m, aes(x=variable, y=value, fill=ClusterID)) + 
    geom_bar(stat="identity") + 
    labs(title = "", y="pct of cells from dataset", x="Dataset") + 
    theme_minimal() +
    theme(
        axis.text=element_text(size=10), 
        plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black")
    )

# -------------------------------------------------------
# Plot expression of known markers
# -------------------------------------------------------
genes <- c(
    "ADGRE1", "ICAM2", "CD93", "NAPSA", "CD226", "CIITA", "CCR2", "RETNLA",
    "MRC1", "TIMD4", "MKI67", "TUBB5", "TOP2A", "BIRC5", "UBE2C", "CD209A",
    "H2.AA", "H2.EB1", "FCRLS"
  )

VlnPlot(sobj.combined, genes, pt.size = 0, ncol = 5)

p <- FeaturePlot(
    sobj.combined, 
    features = genes, 
    pt.size = 0.01, 
    combine = FALSE, 
    min.cutoff = 0, 
    max.cutoff = 3
)

# -------------------------------------------------------
# Find Consensus markers
# -------------------------------------------------------
markers <- FindConservedMarkers(
    sobj.combined, 
    ident.1 = 1, 
    test.use="MAST", 
    only.pos=TRUE, 
    grouping.var = "stim"
) %>%
    rownames_to_column("gene") %>%
    mutate(cluster = 1)

for (ID in c(2, 3, 4, 5, 6)) {
    m <- FindConservedMarkers(
        sobj.combined,
        ident.1 = ID,
        test.use="MAST",
        only.pos=TRUE,
        grouping.var = "stim"
    ) %>%
        rownames_to_column("gene") %>%
        mutate(cluster = ID)
    
    markers <- rbind(markers, m)
}

# Show markers that are in sobj.combined
markers <- markers[which(markers$gene %in% rownames(sobj.combined)), ]
# Select top 10 markers per cluster by p_val_adj
top10 <- markers %>% 
    group_by(cluster) %>% 
    top_n(n = -10, wt = female_p_val_adj)

DoHeatmap(sobj.combined, features = top10$gene, raster = F) + 
    NoLegend() +
    FontSize(y.text = 5)

genes <- unique(top10$gene)
plots <- VlnPlot(sobj.combined, genes, pt.size = 0, combine = F)
for (i in 1:length(genes)) {
    plots[[i]] + 
        NoLegend() + 
        theme(
            axis.text=element_text(size=7),
            plot.margin = unit(c(0,0,0,0), "cm"),
            axis.title = element_text(size = 7),
            title = element_text(size = 7),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour = "black")
        )
}

