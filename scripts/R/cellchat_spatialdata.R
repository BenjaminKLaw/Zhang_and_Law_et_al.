library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(Matrix)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(tidyr)
library(raster)
library(grid)
library(CellChat)


dir <- "path/to/spatial_data"
setwd(dir)

## save individual spatial replicate from integrated data that have cell type annotation and gene names (instead of Ensembl ID)
lung.combined = readRDS("path/to/lung_integrated.RDS")

## save annotated spatial data
lung_spatial <- lung.combined[, lung.combined$tech == 'DBiT_seq']
lung_spatial_rep1 <- lung_spatial[,lung_spatial$orig.ident == "sample_rep1"]
lung_spatial_rep2 <- lung_spatial[,lung_spatial$orig.ident == "sample_rep2"]
lung_spatial_rep3 <- lung_spatial[,lung_spatial$orig.ident == "sample_rep3"]
lung_spatial_rep4 <- lung_spatial[,lung_spatial$orig.ident == "sample_rep4"]
saveRDS(lung_spatial, file="lung_spatial.RDS")
saveRDS(lung_spatial_rep1, file="lung_spatial_rep1.RDS")
saveRDS(lung_spatial_rep2, file="lung_spatial_rep2.RDS")
saveRDS(lung_spatial_rep3, file="lung_spatial_rep3.RDS")
saveRDS(lung_spatial_rep4, file="lung_spatial_rep4.RDS")

# create cellchat object
# preparing replicate data
data.input1 <- GetAssayData(lung_spatial_rep1, layer = "data") # normalized data matrix

# to convert gene name
labels1 <- Idents(lung_spatial_rep1)
meta1 <- data.frame(labels = labels1, slices = "A1") # create a dataframe of the cell labels
meta1$labels = droplevels(meta1$labels, exclude = setdiff(levels(meta1$labels),unique(meta1$labels)))


## load lung_spatial_rep2
data.input2 <- GetAssayData(lung_spatial_rep2, layer = "data") # normalized data matrix
labels2 <- Idents(lung_spatial_rep2)
meta2 <- data.frame(labels = labels2, slices = "A2") # create a dataframe of the cell labels
meta2$labels = droplevels(meta2$labels, exclude = setdiff(levels(meta2$labels),unique(meta2$labels)))

## load lung_spatial_rep3
data.input3 <- GetAssayData(lung_spatial_rep3, layer = "data") # normalized data matrix
labels3 <- Idents(lung_spatial_rep3)
meta3 <- data.frame(labels = labels3, slices = "A3") # create a dataframe of the cell labels
meta3$labels = droplevels(meta3$labels, exclude = setdiff(levels(meta3$labels),unique(meta3$labels)))

## load lung_spatial_rep4
data.input4 <- GetAssayData(lung_spatial_rep4, layer = "data") # normalized data matrix
labels4 <- Idents(lung_spatial_rep4)
meta4 <- data.frame(labels = labels4, slices = "A4") # create a dataframe of the cell labels
meta4$labels = droplevels(meta4$labels, exclude = setdiff(levels(meta4$labels),unique(meta4$labels)))

genes.common <- intersect(intersect(intersect(rownames(data.input1), rownames(data.input2)),rownames(data.input3)),rownames(data.input4))
data.input <- cbind(data.input1[genes.common, ], data.input2[genes.common, ], data.input3[genes.common, ], data.input4[genes.common, ])
meta <- rbind(meta1, meta2, meta3,meta4)
rownames(meta) <- colnames(data.input)
meta$labels <- droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
meta$slices <- factor(meta$slices, levels = c("A1", "A2", "A3","A4"))
unique(meta$labels) # check the cell labels

# load spatial imaging information
id1 <- Idents(lung_spatial_rep1)
id2 <- Idents(lung_spatial_rep2)
id3 <- Idents(lung_spatial_rep3)
id4 <- Idents(lung_spatial_rep4)
df1 <- data.frame(V1 = names(id1),V2=id1 )
df2 <- data.frame(V1 = names(id2),V2=id2 )
df3 <- data.frame(V1 = names(id3),V2=id3 )
df4 <- data.frame(V1 = names(id4),V2=id4 )
test1 <- df1 %>% separate(V1, c("A", "B"),  sep = "x")
test1 <- test1 %>% separate(B,c("B","C"), sep = "_")
test2 <- df2 %>% separate(V1, c("A", "B"),  sep = "x")
test2 <- test2 %>% separate(B,c("B","C"), sep = "_")
test3 <- df3 %>% separate(V1, c("A", "B"),  sep = "x")
test3 <- test3 %>% separate(B,c("B","C"), sep = "_")
test4 <- df4 %>% separate(V1, c("A", "B"),  sep = "x")
test4 <- test4 %>% separate(B,c("B","C"), sep = "_")
spatial.locs1 <- test1[,c(1,2)]
spatial.locs2 <- test2[,c(1,2)]
spatial.locs3 <- test3[,c(1,2)]
spatial.locs4 <- test4[,c(1,2)]
spatial.locs <- rbind(spatial.locs1, spatial.locs2, spatial.locs3,spatial.locs4)
rownames(spatial.locs) <- colnames(data.input)

multiplier <- 25
spatial.locs$A <- as.numeric(spatial.locs$A)
spatial.locs$B <- as.numeric(spatial.locs$B)
spatial.locs_multiplied <- spatial.locs %>% mutate(across(where(is.numeric), ~ . * multiplier))
conversion.factor1 = 0.8; 
spot.size1 = 10
spatial.factors1 = data.frame(ratio = conversion.factor1, tol = spot.size1/2)
spatial.factors <- rbind(spatial.factors1,spatial.factors1,spatial.factors1,spatial.factors1)
rownames(spatial.factors) <- c("A1","A2","A3","A4")
# finally create cellchat object
cellchat4 <- createCellChat(object = data.input, meta = meta, group.by = "labels", datatype = "spatial", spatial.factors = spatial.factors, coordinates = spatial.locs_multiplied)
cellchat4

# set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat4@DB <- CellChatDB.use

# Identify over-expressed ligands or receptors.
cellchat4 <- subsetData(cellchat4) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat4 <- identifyOverExpressedGenes(cellchat4)
cellchat4 <- identifyOverExpressedInteractions(cellchat4)

# Infer cell-cell communication at a ligand-receptor pair level
cellchat4 <- computeCommunProb(cellchat4, type = "truncatedMean", trim = 0.1,distance.use = TRUE, interaction.range = 250, scale.distance = 3.9, contact.dependent = TRUE, contact.range = 20)
# Filter the cell-cell communication
cellchat4 <- filterCommunication(cellchat4, min.cells = 10)

# Infer cell-cell communication at a signaling pathway level.
cellchat4 <- computeCommunProbPathway(cellchat4)

# Calculate aggregated cell-cell communication network
cellchat4 <- aggregateNet(cellchat4)
# Compute the network centrality scores
cellchat4 <- netAnalysis_computeCentrality(cellchat4, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#  export the CellChat object together with the inferred cell- cell communication networks 
saveRDS(cellchat4, file = "cellchat_output.rds")

## visulaization
## Visualize differential interactions between different celltypes
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ptm = Sys.time()
netVisual_heatmap(cellchat4,targets.use = "Neuron", measure = "weight")
dev.copy(pdf, "Neuron interaction strength with other cell types.pdf", width = 8, height = 6)
dev.off()  # Close the device to save the file

netVisual_heatmap(cellchat4, measure = "weight")
dev.copy(pdf, "interaction strength with different cell types.pdf", width = 8, height = 6)
dev.off()  # Close the device to save the file

# check the significant interaction between each cell type with others
celllist <- unique(meta$labels)
for (i in 1:length(celllist)) {
  # bubble plot
  netVisual_bubble(cellchat4, sources.use = celllist[i], targets.use = c(1:20), remove.isolate = FALSE)
  # save
  ggsave(filename=paste0(celllist[i], "_SignificantInteraction.pdf"), width = 5, height = 20, units = 'in', dpi = 300)
}
pathways.show.all <- cellchat4@netP$pathways

# Hierarchy plot for each significant interactions
for (i in 1:length(pathways.show.all)) {
  # Hierarchy plot
  netVisual_aggregate(cellchat4, signaling = pathways.show.all[i])
  # save
  filename=paste0("sample_",pathways.show.all[i], "_Hierarchy.pdf")
  dev.copy(pdf, filename)
  dev.off()
}
# Centrality score plot for each significant interactions
for (i in 1:length(pathways.show.all)) {
  # Hierarchy plot
  netAnalysis_signalingRole_network(cellchat4, signaling = pathways.show.all[i], width = 8, height = 2.5, font.size = 10)
  # save
  filename=paste0("sample_",pathways.show.all[i], "_Centrality score.pdf")
  dev.copy(pdf, filename)
  dev.off()
}
# Heatmap plot for each significant interactions
for (i in 1:length(pathways.show.all)) {
  # Heatmap
  myplot <- netVisual_heatmap(cellchat4, signaling = pathways.show.all[i], color.heatmap = "Reds", width = 10, height = 5, font.size = 10)
  print(myplot)
  #netAnalysis_signalingRole_network(cellchat4, signaling = pathways.show.all[i], width = 8, height = 2.5, font.size = 10)
  # save
  filename=paste0("sample_",pathways.show.all[i], "_Heatmap.pdf")
  dev.copy(pdf, filename)
  dev.off()
}
# Contribution of each signaling pathway
for (i in 1:length(pathways.show.all)) {
  # contribution
  myplot <- netAnalysis_contribution(cellchat4, signaling = pathways.show.all[i])
  print(myplot)
  #netAnalysis_signalingRole_network(cellchat4, signaling = pathways.show.all[i], width = 8, height = 2.5, font.size = 10)
  # save
  filename=paste0("sample_",pathways.show.all[i], "_Contribution.pdf")
  dev.copy(pdf, filename)
  dev.off()
}
