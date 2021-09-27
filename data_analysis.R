library(Seurat)
library(dplyr)
library(ggsci)
library(gridExtra)
source("~/Google Drive/bin/batadalab_scrnaseq_utils.R")

# -------------------------------------------------------
# Read data
# -------------------------------------------------------
female_data = get_scrnaseq_data("p01e21/jenkins_female", from_RDS = F)
male_data = get_scrnaseq_data("p01e20/jenkins_male", from_RDS = F)

gem_female = female_data$gem
gem_male = male_data$gem

sobj1 <- CreateSeuratObject(counts = gem_female, project = "female", min.cells = round(0.01*ncol(gem_female)))
sobj1$stim <- "female"
sobj1 <- subset(sobj1, subset = nFeature_RNA > 500)
sobj1 <- NormalizeData(sobj1, verbose = FALSE)
sobj1 <- FindVariableFeatures(sobj1, selection.method = "vst", nfeatures = 2000)

sobj2 <- CreateSeuratObject(counts = gem_male, project = "male", min.cells = round(0.01*ncol(gem_male)))
sobj2$stim <- "male"
sobj2 <- subset(sobj2, subset = nFeature_RNA > 500)
sobj2 <- NormalizeData(sobj2, verbose = FALSE)
sobj2 <- FindVariableFeatures(sobj2, selection.method = "vst", nfeatures = 2000)

# -------------------------------------------------------
# Integrate
# -------------------------------------------------------
anchors <- FindIntegrationAnchors(object.list = list(sobj1, sobj2), dims = 1:20)
sobj.combined <- IntegrateData(anchorset = anchors, dims = 1:20)

sobj.combined@meta.data$ref.ident = refLabs[rownames(sobj.combined@meta.data)]

# -------------------------------------------------------
# Perform integrated analysis
# -------------------------------------------------------
DefaultAssay(sobj.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
sobj.combined <- ScaleData(sobj.combined, verbose = FALSE, features = rownames(sobj.combined))
sobj.combined <- RunPCA(sobj.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
sobj.combined <- RunUMAP(sobj.combined, reduction = "pca", dims = 1:20)
sobj.combined <- FindNeighbors(sobj.combined, reduction = "pca", dims = 1:20)
sobj.combined <- FindClusters(sobj.combined, resolution = 0.5)

sobj.combined = subset(sobj.combined, idents = c(0,1,2,3,4,5,6,7,8,9,10))

sobj.combined = BuildClusterTree(sobj.combined)
PlotClusterTree(sobj.combined)

# Merge clusters 0+1+2+5 and 3+10
Idents(sobj.combined)[WhichCells(sobj.combined, idents = 0)] = 1
Idents(sobj.combined)[WhichCells(sobj.combined, idents = 2)] = 1
Idents(sobj.combined)[WhichCells(sobj.combined, idents = 5)] = 1
Idents(sobj.combined)[WhichCells(sobj.combined, idents = 10)] = 3
Idents(sobj.combined)[WhichCells(sobj.combined, idents = 9)] = 3

# Rename idents
new.cluster.ids <- c(5, 4, 3, 6, 2, 1)
names(new.cluster.ids) <- levels(sobj.combined)
sobj.combined <- RenameIdents(sobj.combined, new.cluster.ids)
Idents(sobj.combined) = factor(Idents(sobj.combined), levels = c(1,2,3,4,5,6))

p1 = DimPlot(sobj.combined, label = T) + NoAxes() + NoLegend()
p2 = DimPlot(sobj.combined, group.by = "stim", split.by = 'stim') + NoAxes() + NoLegend()
gridExtra::grid.arrange(p1, p2)

pdf("~/Google Drive/Data/Analysis/collaborations/jenkins/Sept2019/plots/24Sept/UMAP.pdf", width = 3, height = 3)
DimPlot(sobj.combined, label = T) + NoAxes() + NoLegend()
dev.off()

pdf("~/Google Drive/Data/Analysis/collaborations/jenkins/Sept2019/plots/24Sept/UMAP_individuals.pdf", width = 6, height = 3)
DimPlot(sobj.combined, group.by = "stim", split.by = 'stim') + NoAxes() + NoLegend()
dev.off()

sobj.combined = BuildClusterTree(sobj.combined)
#pdf("~/Google Drive/Data/Analysis/collaborations/jenkins/Sept2019/plots/24Sept/ClusterTree.pdf")
PlotClusterTree(sobj.combined)
#dev.off()

# -------------------------------------------------------
# Cluster Frequency per sample
# -------------------------------------------------------
stats = data.frame(ClusterID = unique(Idents(sobj.combined)), pct_Female = NA, pct_Male = NA)
for (i in 1:nrow(stats)) {
  ID = stats[i, "ClusterID"]
  cells = WhichCells(sobj.combined, idents = ID)
  stats[i, "pct_Female"] = round(length(grep("_1", cells))*100/length(grep("_1", colnames(sobj.combined))), 2)
  stats[i, "pct_Male"] = round(length(grep("_2", cells))*100/length(grep("_2", colnames(sobj.combined))), 2)
}
df.m = reshape::melt(stats)

pdf("~/Google Drive/Data/Analysis/collaborations/jenkins/Sept2019/plots/24Sept/Cluster_frequencies.pdf", width = 2, height = 3)
ggplot(df.m, aes(x=variable, y=value, fill=ClusterID)) + geom_bar(stat="identity") + 
  labs(title = "", y="pct of cells from dataset", x="Dataset") + theme_minimal() +
  theme(axis.text=element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black"))
dev.off()

# -------------------------------------------------------
# Plot expression of known markers
# -------------------------------------------------------
# Annotate clusters
#pdf("~/Google Drive/Data/Analysis/collaborations/jenkins/Sept2019/plots/24Sept/Prolif_markers.pdf")
#VlnPlot(sobj.combined, c("MKI67", "TUBB5", "TOP2A", "BIRC5", "UBE2C"), pt.size = 0)
#dev.off()

#pdf("~/Google Drive/Data/Analysis/collaborations/jenkins/Sept2019/plots/24Sept/pan_markers.pdf")
#VlnPlot(sobj.combined, c("ADGRE1", "ICAM2", "CD93"), pt.size = 0)
#dev.off()

#pdf("~/Google Drive/Data/Analysis/collaborations/jenkins/Sept2019/plots/24Sept/mono_markers.pdf")
#VlnPlot(sobj.combined, c("MRC1", "CCR2", "RETNLA"), pt.size = 0)
#dev.off()

# -------------------------------------------------------
# Plot expression of known markers
# -------------------------------------------------------
genes = c("ADGRE1", "ICAM2", "CD93", "NAPSA", "CD226", "CIITA", "CCR2", "RETNLA", "MRC1", "TIMD4", "MKI67", "TUBB5", 
          "TOP2A", "BIRC5", "UBE2C", "CD209A", "H2.AA", "H2.EB1", "FCRLS")

pdf("~/Google Drive/Data/Analysis/collaborations/jenkins/Sept2019/plots/24Sept/lineage_markers_vln.pdf", width = 8, height = 8)
VlnPlot(sobj.combined, genes, pt.size = 0, ncol = 5)
dev.off()

p = FeaturePlot(sobj.combined, features = genes, pt.size = 0.01, combine = FALSE, min.cutoff = 0, max.cutoff = 3)
#DoHeatmap(sobj.combined, genes)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + NoLegend()
  pdf(paste("~/Google Drive/Data/Analysis/collaborations/jenkins/Sept2019/plots/24Sept/FeaturePlots/", genes[i], ".pdf", sep = ""), width = 1, height = 1)
  plot(p[[i]])
  dev.off()
}

# -------------------------------------------------------
# Find Consensus markers
# -------------------------------------------------------
markers = FindConservedMarkers(sobj.combined, ident.1 = 1, test.use="MAST", only.pos=TRUE, grouping.var = "stim")
markers$gene = rownames(markers)
markers$cluster = 1
for (ID in c(2, 3, 4, 5, 6)) {
  m = FindConservedMarkers(sobj.combined, ident.1 = ID, test.use="MAST", only.pos=TRUE, grouping.var = "stim")
  m$gene = rownames(m)
  m$cluster = ID
  markers = rbind(markers, m)
}

write.table(markers, file = "~/Google Drive/Data/Analysis/collaborations/jenkins/Sept2019/data/24Sept/Consensus_markers.txt", quote = F, sep = "\t")

# Show markers that are in sobj.combined
markers = markers[which(markers$gene %in% rownames(sobj.combined)), ]
# Select top 10 markers per cluster by p_val_adj
top30 = markers %>% group_by(cluster) %>% top_n(n = -30, wt = female_p_val_adj)

tiff("~/Google Drive/Data/Analysis/collaborations/jenkins/Sept2019/plots/24Sept/Heatmap.tiff")
DoHeatmap(sobj.combined, features = top30$gene, raster = F) + NoLegend() + FontSize(y.text = 5)
dev.off()

genes = unique(top30$gene)
plots = VlnPlot(sobj.combined, genes, pt.size = 0, combine = F)
for (i in 1:length(genes)) {
  pdf(paste("~/Google Drive/Data/Analysis/collaborations/jenkins/Sept2019/plots/24Sept/VlnPlots/", genes[i], ".pdf", sep = ""), width = 1, height = 1)
  print(plots[[i]] + NoLegend() + theme(axis.text=element_text(size=7), plot.margin = unit(c(0,0,0,0), "cm"),
                                        axis.title = element_text(size = 7), title = element_text(size = 7),
                                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       panel.background = element_blank(), axis.line.x = element_line(colour = "black"),
                                       axis.line.y = element_line(colour = "black"))) 
  dev.off()
}

saveRDS(sobj.combined, file="~/Google Drive/Data/Analysis/collaborations/jenkins/Sept2019/data/24Sept/sobj.combined.rds")


