## CARD deconvolution
## Benjamin Law
## 20240128
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
#library(rhdf5)
library(Matrix)
library(sctransform)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(raster)
library(OpenImageR)
library(ggpubr)
library(grid)
library(wesanderson)
library(CARD)
library(MuSiC)
library(SingleCellExperiment)

setwd('~/ChanLab/Spatial_Lung/sandbox/scRNA/CARD/')

# load spatial data
data <- readRDS("~/ChanLab/Spatial_Lung/raw_data/scRNA/integrated/Lung_dbit_srf_CellReports_Science_integrated.RDS")

# Run the analysis on each dataset invididually
################################################################################
# Analyze the E11 dataset via card
Idents(data) <- data$seurat_clusters

datasets <- c('E11lungR1', 'E11lungR2', 'E11lungR3', 'ElllungR4', 'E11lungR5', 'E12lungR1', 'E12lungR2', 'E12lungR3', 'E12lungR4', 'E12lungR5', 'E13lungR1', 'E13lungR2', 'E13lungR3', 'E13lungR4')

# load the reference datasets
scRNA_seq <- readRDS("~/ChanLab/Spatial_Lung/raw_data/scRNA/integrated/Lung_srf_CellReports_Science_integrated_reference_SCT.RDS")
DimPlot(scRNA_seq, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# change to gene expression matrix
scRNA_seqdata <- scRNA_seq@assays$SCT$counts
scRNA_seqdata[1:4,1:4]

# keep the named idents rather than the seurat cluster numbers
#Idents(scRNA_seq) <- scRNA_seq$seurat_clusters

# add cell type information to metadata & load metadata
cell_type.info <- data.frame(cell_type = scRNA_seq@active.ident, row.names= colnames(scRNA_seq))
scRNA_seq <- AddMetaData(object = scRNA_seq, metadata = cell_type.info)
scRNA_seqmeta <- scRNA_seq@meta.data
scRNA_seqmeta$sampleInfo = "sample1"
scRNA_seqmeta[1:10,]
options(future.globals.maxSize= 10891289600)

# create an environment of objects for each of the datasets
for (exp in datasets){
  exp = 'E13lungR4'
  temp_data <- subset(x = data, subset = orig.ident == exp)

  # Make the spatial transcriptomic dataset. This should be a matrix with rows as genes and columns as a spatial location. These are the cellBC from the integrated dataset
  spatial_count = temp_data@assays$SCT@counts
  spatial_count[1:4, 1:4]

  # make a dataframe that has the spatial x and y coordinates
  locdata <- data.frame(row.names = colnames(temp_data))
  temp <- data.frame(X = colnames(temp_data))
  temp <- temp %>% separate(X, c("A", "B"), sep = "x")
  temp <- temp %>% separate(B, c("B", "C", "D"), sep = "_")
  locdata$x=as.numeric(as.character(temp$A))
  locdata$y=as.numeric(as.character(temp$B))

  # Create the CARD object
  CARD_obj = createCARDObject(sc_count = scRNA_seqdata, sc_meta = scRNA_seqmeta, spatial_count = spatial_count, spatial_location = locdata, ct.varname = "cell_type", ct.select = unique(scRNA_seqmeta$cell_type), sample.varname = 'sampleInfo', minCountGene = 100, minCountSpot = 5)
  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
  
  proportions = CARD_obj@Proportion_CARD
  write.csv(proportions, file = paste0(exp, '_CARD_proportions.csv'))
  
  # open a pdf for the results
  pdf(paste0(exp, '_CARD_plots.pdf'))
  ## set the colors. Here, I just use the colors in the manuscript, if the color is not provided, the function will use default color in the package. 
  colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
             "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
             "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
  
  p1 <- CARD.visualize.pie(proportion = CARD_obj@Proportion_CARD,
                           spatial_location = CARD_obj@spatial_location, 
                           colors = colors,
                           radius = 0.5)
  print(p1)

  ## select the cell type that we are interested
  ct.visualize = c('Epithelial', 'MesenchymeProgenitor', 'Mesenchyme1', 'SmoothMuscle', 'Mesothelium', 'VascularEndothelial', 'Immune', 'Precartilage', 'Neuron')
  ## visualize the spatial distribution of the cell type proportion
  p2 <- CARD.visualize.prop(
    proportion = CARD_obj@Proportion_CARD,        
    spatial_location = CARD_obj@spatial_location, 
    ct.visualize = ct.visualize,                 ### selected cell types to visualize
    colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
    NumCols = 6,
    pointSize = 1)                                 ### number of columns in the figure panel
  print(p2)
  
  # visualize the cell type proportion correlation
  p3 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
  print(p3)
  
  ## refined spatial map
  # imputation on the newly grided spatial locations
  CARD_obj = CARD.imputation(CARD_obj, NumGrids = 2500, ineibor = 10, exclude = NULL)
  
  ## Visualize the newly grided spatial locations to see if the shape is correctly detected. If not, the user can provide the row names of the excluded spatial location data into the CARD.imputation function
  location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
                                         y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
  rownames(location_imputation) = rownames(CARD_obj@refined_prop)
  
  p4<- ggplot(location_imputation, 
         aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          legend.position="bottom",
          panel.background = element_blank(),
          plot.background = element_blank(),
          panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))
  print(p4)
  
  # Visualize the cell type proportion at an enhanced resolution
  p5 <- CARD.visualize.prop(
    proportion = CARD_obj@refined_prop,                         
    spatial_location = location_imputation,            
    ct.visualize = ct.visualize,                    
    colors = c("lightblue","lightyellow","red"),    
    NumCols = 6,
    pointSize = 1)                                  
  print(p5)
  
  dev.off()
  
  #Save the environment for each dbit sample
  save(CARD_obj, file = paste0("Lung_CARD_", exp, ".RData"))
}
