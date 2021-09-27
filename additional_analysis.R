library(Seurat)
library(dplyr)
library(ggsci)
library(gridExtra)
library(scID)
source("~/Google Drive/bin/batadalab_scrnaseq_utils.R")

# -------------------------------------------------------
# Read data
# -------------------------------------------------------
sobj.combined <- readRDS("~/Google Drive/Data/Analysis/collaborations/jenkins/Sept2019/data/24Sept/sobj.combined.rds")

# -------------------------------------------------------
# Cluster Frequency per sample - no clusters 1,2,6
# -------------------------------------------------------
stats <- data.frame(
  ClusterID = factor(c(3,4,5)), 
  pct_Female = NA, 
  pct_Male = NA
)
all_cells <- WhichCells(sobj.combined, idents = c(3,4,5))
for (i in 1:nrow(stats)) {
  ID <- stats[i, "ClusterID"]
  cells <- WhichCells(sobj.combined, idents = ID)
  stats[i, "pct_Female"] <- round(length(grep("_1", cells))*100/length(grep("_1", all_cells)), 2)
  stats[i, "pct_Male"] <- round(length(grep("_2", cells))*100/length(grep("_2", all_cells)), 2)
}
df.m <- reshape::melt(stats)

pdf("~/Google Drive/Data/Analysis/collaborations/jenkins/Feb2020/plots/Cluster_frequencies_345.pdf", width = 3, height = 3)
ggplot(df.m, aes(x=variable, y=value, fill=ClusterID)) + geom_bar(stat="identity") + 
  labs(title = "", y="pct of cells from dataset", x="Dataset") + theme_minimal() +
  theme(axis.text=element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black"))
dev.off()

# -------------------------------------------------------
# Cluster Frequency per sample - no clusters 1,2
# -------------------------------------------------------
stats <- data.frame(ClusterID = factor(c(3,4,5,6)), pct_Female = NA, pct_Male = NA)
all_cells = WhichCells(sobj.combined, idents = c(3,4,5,6))
for (i in 1:nrow(stats)) {
  ID = stats[i, "ClusterID"]
  cells = WhichCells(sobj.combined, idents = ID)
  stats[i, "pct_Female"] = round(length(grep("_1", cells))*100/length(grep("_1", all_cells)), 2)
  stats[i, "pct_Male"] = round(length(grep("_2", cells))*100/length(grep("_2", all_cells)), 2)
}
df.m = reshape::melt(stats)

pdf("~/Google Drive/Data/Analysis/collaborations/jenkins/Feb2020/plots/Cluster_frequencies_3456.pdf", width = 3, height = 3)
ggplot(df.m, aes(x=variable, y=value, fill=ClusterID)) + geom_bar(stat="identity") + 
  labs(title = "", y="pct of cells from dataset", x="Dataset") + theme_minimal() +
  theme(axis.text=element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black"))
dev.off()

# -------------------------------------------------------
# Map proliferating cells to combined clusters
# -------------------------------------------------------
labels = Idents(sobj.combined)

female_labels = labels[grep("_1", names(labels))]
names(female_labels) = unlist(lapply(names(female_labels), function(x) strsplit(x, "_")[[1]][1]))

male_labels = labels[grep("_2", names(labels))]
names(male_labels) = unlist(lapply(names(male_labels), function(x) strsplit(x, "_")[[1]][1]))

# run scID for female cells
data_f = get_scrnaseq_data('p01e21/jenkins_female', from_RDS = F)
gem_f = data_f$gem

reference_fem = gem_f[, names(female_labels)[which(female_labels != 6)]]
target_fem = gem_f[, names(female_labels)[which(female_labels == 6)]]
target_fem_cpm = counts_to_cpm(target_fem)

female_res = scid_multiclass(target_gem = target_fem_cpm, reference_gem = reference_fem, 
                             reference_clusters = female_labels[colnames(reference_fem)], 
                             only_pos = F, estimate_weights_from_target = F)
table(female_res$labels)
table(female_res$markers$cluster)

make_heatmap(target_fem_cpm, female_res$labels, female_res$markers)
make_heatmap(counts_to_cpm(reference_fem), female_labels[colnames(reference_fem)], female_res$markers)

# run scID for male cells
data_m = get_scrnaseq_data('p01e20/jenkins_male', from_RDS = F)
gem_m = data_m$gem

reference_m = gem_m[, names(male_labels)[which(male_labels != 6)]]
target_m = gem_m[, names(male_labels)[which(male_labels == 6)]]
target_m_cpm = counts_to_cpm(target_m)

male_res = scid_multiclass(target_gem = target_m_cpm, reference_gem = reference_m, 
                           reference_clusters = male_labels[colnames(reference_m)], 
                           only_pos = F, estimate_weights_from_target = F)
table(male_res$labels)
table(male_res$markers$cluster)

make_heatmap(target_m_cpm, male_res$labels, male_res$markers)
make_heatmap(counts_to_cpm(reference_m), male_labels[colnames(reference_m)], male_res$markers)

# -------------------------------------------------------
# Repeat Cluster Frequency per sample (no clusters 1,2,6) after reassignment of Cluster 6 cells
# -------------------------------------------------------
sobj.combined.new = sobj.combined
Idents(sobj.combined.new)[paste(names(female_res$labels[which(female_res$labels == 3)]), "1", sep = "_")] = 3
Idents(sobj.combined.new)[paste(names(female_res$labels[which(female_res$labels == 4)]), "1", sep = "_")] = 4
Idents(sobj.combined.new)[paste(names(female_res$labels[which(female_res$labels == 5)]), "1", sep = "_")] = 5

Idents(sobj.combined.new)[paste(names(male_res$labels[which(male_res$labels == 3)]), "2", sep = "_")] = 3
Idents(sobj.combined.new)[paste(names(male_res$labels[which(male_res$labels == 4)]), "2", sep = "_")] = 4
Idents(sobj.combined.new)[paste(names(male_res$labels[which(male_res$labels == 5)]), "2", sep = "_")] = 5

table(Idents(sobj.combined))
table(Idents(sobj.combined.new))

stats = data.frame(ClusterID = factor(c(3,4,5)), pct_Female = NA, pct_Male = NA)
all_cells = WhichCells(sobj.combined, idents = c(3,4,5))
for (i in 1:nrow(stats)) {
  ID = stats[i, "ClusterID"]
  cells = WhichCells(sobj.combined, idents = ID)
  stats[i, "pct_Female"] = round(length(grep("_1", cells))*100/length(grep("_1", all_cells)), 2)
  stats[i, "pct_Male"] = round(length(grep("_2", cells))*100/length(grep("_2", all_cells)), 2)
}
df.m = reshape::melt(stats)

pdf("~/Google Drive/Data/Analysis/collaborations/jenkins/Feb2020/plots/Cluster_frequencies_345_after_reassignment.pdf", width = 3, height = 3)
ggplot(df.m, aes(x=variable, y=value, fill=ClusterID)) + geom_bar(stat="identity") + 
  labs(title = "", y="pct of cells from dataset", x="Dataset") + theme_minimal() +
  theme(axis.text=element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black"))
dev.off()

# find what proportion of cells in each cluster are in cycle
stats = data.frame(ClusterID = factor(c(3,4,5)), pct_Female = NA, pct_Male = NA)
for (i in 1:nrow(stats)) {
  ID = stats[i, "ClusterID"]
  cells = WhichCells(sobj.combined, idents = ID)
  stats[i, "pct_Female"] = round(length(female_res$labels[which(female_res$labels == ID)])*100/length(grep("_1", cells)), 2)
  stats[i, "pct_Male"] = round(length(male_res$labels[which(male_res$labels == ID)])*100/length(grep("_2", cells)), 2)
}
df.m = reshape::melt(stats)

pdf("~/Google Drive/Data/Analysis/collaborations/jenkins/Feb2020/plots/percent_cycling_after_reassignment.pdf", width = 3, height = 3)
ggplot(df.m, aes(x=ClusterID, y=value, fill=variable)) + geom_bar(position="dodge", stat="identity") + 
  labs(title = "", y="pct of cells from dataset", x="Dataset") + theme_minimal() +
  theme(axis.text=element_text(size=10), plot.margin = unit(c(0,0,0,0), "cm"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black"))
dev.off()

saveRDS(sobj.combined.new, file = "~/Google Drive/Data/Analysis/collaborations/jenkins/Feb2020/sobj_combined_reassigned.rds")

