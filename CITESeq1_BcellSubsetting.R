## ##Re-clustering of all B cells in the CITE-Seq (1) dataset####
####Setup####
#Set working directory to folder containing files

##Packages##
library(Seurat)
library(Polychrome)
library(SeuratDisk)
library(ggplot2)
library(viridis)
library(clustree)
library(patchwork)
library(stringr)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(HGNChelper)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(sc2marker)
sessionInfo()

##Potential Packages##
library(clustifyr)
library(ArrayExpress)
library(GEOquery)
library(clustifyrdatahub)
library(devtools)
library(rjags)
library(infercnv)
library(clustree)
library(enrichR)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AnnotationHub)
library(scRepertoire)
library(decoupleR)
library(OmnipathR)
library(grid)
library(gridExtra)
library(msigdbr)
library(fgsea)
library(presto)
library(SCpubr)
library(EnhancedVolcano)
library(ggraph)
library(RColorBrewer)
library(openxlsx)
library(ensembldb)
library(DOSE)
library(Ibex)
library(HGNChelper)
library(openxlsx)
library(VennDiagram)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(sc2marker)
sessionInfo()

##Colour paneles##
col_con1 <- createPalette(50,  c("#2A9D8F", "#E9C46A", "#E76F51"))
col_con1 <-as.character(col_con1)
col_con2 <- c("#C94055", "#EFE3A2", "#4FA4C2", "#E80DFF", "#26FF0D", "#0D65FB", "#0DF6BD", "#F8000D", "#FE0DBD", "#FDA30D",
              "#FCCFEB", "#B4F2A7", "#C5F00D", "#FC6A00", "#F8AA92", "#94B2F4", "#FE8BBD", "#F24798", "#22FAFB", "#B898FB") #Adjusted for 20x bright colours
col_con3 <- createPalette(50,  c("#CE4257", "#FFF2B0", "#4BA3C3"))
col_con3 <-as.character(col_con2)
col_con.genotypes <- c("#666666", "#f2dd38", "#4ba3de", "#14905f") #For Genotypes
col_con.mouse <- c("grey1", "azure4", "gold1", "chocolate2", "blue", "forestgreen", "darkolivegreen1")#For individual mice
col_con.sex <- c("darkred", "darkturquoise") #For Sex
col_con_clonotype <- rev(c("grey39", "#4FA4C2", "chartreuse4", "gold2", "darkturquoise","firebrick3")) #Clonotype all
col_con.clonotype.no.hyper <- c( "darkturquoise", "gold2", "chartreuse4", "grey39") #Clonotype no hyperexpansion

##Vector orders##
clo_order <- c("NA", "Rare (0 < X <= 1e-04)", "Small (1e-04 < X <= 0.001)", "Medium (0.001 < X <= 0.01)", "Large (0.01 < X <= 0.1)", "Hyperexpanded (0.1 < X <= 1)")
clo_order.reverse <- c("Hyperexpanded (0.1 < X <= 1)", "Large (0.01 < X <= 0.1)", "Medium (0.001 < X <= 0.01)", "Small (1e-04 < X <= 0.001)","Rare (0 < X <= 1e-04)", "NA")
gen_order <- c("WT", "BCL6", "E1020K", "E1020K_BCL6")
mouse_order <- c("WT_1", "WT_2", "BCL6_1", "BCL6_2", "E1020K", "E1020K_BCL6_1", "E1020K_BCL6_2")

##Functions##
tfidf = function(data,target,universe){
  if(!all(target %in% universe))
    stop('Target must be a subset of universe')
  nObs = Matrix::rowSums(data[,target,drop=FALSE]>0)
  nTot = Matrix::rowSums(data[,universe,drop=FALSE]>0)
  tf = nObs/length(target)
  idf = log(length(universe)/nTot)
  score = tf*idf
  #Calculate p-value for significance based on using a hypergeometric distribution to simulate the results of infinite random sampling
  pvals = phyper(nObs-1,nTot,length(universe)-nTot,length(target),lower.tail=FALSE)
  qvals = p.adjust(pvals,method='BH')
  ntf = (exp(-idf)*length(universe)-tf*length(target))/(length(universe)-length(target))
  return(data.frame(geneFrequency=tf,
                    geneFrequencyOutsideCluster=ntf,
                    geneFrequencyGlobal=exp(-idf),
                    geneExpression=Matrix::rowMeans(data[,target,drop=FALSE]),
                    geneExpressionOutsideCluster = Matrix::rowMeans(data[,universe[!(universe%in%target)],drop=FALSE]),
                    geneExpressionGlobal = Matrix::rowMeans(data),
                    idf=idf,
                    tfidf=score,
                    qval=qvals)[order(score,decreasing=TRUE),])
}


####Loading Seurat Objects#####
#Load Seurat Object containing all cells with clonotype data lacking NAs - "Allcells"
Allcells <- LoadH5Seurat(file.choose("CITESeq1_AllCells_nNA.h5seurat")) #With clonotype data, excluding NAs
Allcells$seurat_clusters <- Allcells$wsnn_res.0.8
Idents(Allcells) <- Allcells$seurat_clusters

#DimPlot Allcells
names(col_con1) <- levels(Allcells$seurat_clusters) #One colour associated with one cluster
UMAP1 <- DimPlot(Allcells, reduction = "wnn.umap", cols = col_con1, pt.size = 1, label = TRUE, label.box = TRUE, label.size = 8) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("All cells") +
  theme(plot.title = element_text(color="black", size= 30, face="bold"),
        legend.text = element_text(size = 20),
        text = element_text(size = 20),
        legend.position = "bottom") +
  NoLegend() #Can remove if legend is needed
print(UMAP1)
ggsave("UMAP1.tiff", width = 30, height = 30, units = "cm", UMAP1, compression = "lzw")


####Subsetting into B, T and Other cells####
DefaultAssay(Allcells) <- "ADT"
Unmapped <- FeaturePlot(Allcells, features = "Unmapped", reduction = "wnn.umap", pt.size = 1, order = FALSE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("Unmapped binding") +
  scale_colour_gradientn(colours = magma(10)) +
  theme(plot.title = element_text(color="black", size=30, face="bold"),
        legend.text = element_text(size = 20),
        text = element_text(size = 30),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
print(Unmapped)
ggsave("Unmapped.tiff", width = 30, height = 30, units = "cm", Unmapped, compression = "lzw")

CD19 <- FeaturePlot(Allcells, features = "Cd19", reduction = "wnn.umap", pt.size = 1, order = FALSE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD19 cell surface expression") +
  scale_colour_gradientn(colours = magma(10)) +
  theme(plot.title = element_text(color="black", size=30, face="bold"),
        legend.text = element_text(size = 20),
        text = element_text(size = 30),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
print(CD19)
ggsave("CD19.tiff", width = 30, height = 30, units = "cm", CD19, compression = "lzw")

CD4 <- FeaturePlot(Allcells, features = "Cd4", reduction = "wnn.umap", pt.size = 1, order = FALSE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD4 cell surface expression") +
  scale_colour_gradientn(colours = magma(10)) +
  theme(plot.title = element_text(color="black", size=30, face="bold"),
        legend.text = element_text(size = 20),
        text = element_text(size = 30),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
print(CD4)
ggsave("CD4.tiff", width = 30, height = 30, units = "cm", CD4, compression = "lzw")

CD8a <- FeaturePlot(Allcells, features = "Cd8a", reduction = "wnn.umap", pt.size = 1, order = FALSE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD8a cell surface expression") +
  scale_colour_gradientn(colours = magma(10)) +
  theme(plot.title = element_text(color="black", size=30, face="bold"),
        legend.text = element_text(size = 20),
        text = element_text(size = 30),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
print(CD8a)
ggsave("CD8a.tiff", width = 30, height = 30, units = "cm", CD8a, compression = "lzw")

DefaultAssay(Allcells) <- "RNA"
Ncr1 <- FeaturePlot(Allcells, features = "Ncr1", reduction = "wnn.umap", pt.size = 1, order = FALSE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("Ncr1 gene expression") +
  scale_colour_gradientn(colours = mako(10)) +
  theme(plot.title = element_text(color="black", size=30, face="bold"),
        legend.text = element_text(size = 20),
        text = element_text(size = 30),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
print(Ncr1)
ggsave("Ncr1.tiff", width = 30, height = 30, units = "cm", Ncr1, compression = "lzw")


Prdm1 <- FeaturePlot(Allcells, features = "Prdm1", reduction = "wnn.umap", pt.size = 1, order = FALSE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("Prdm1 gene expression") +
  scale_colour_gradientn(colours = mako(10)) +
  theme(plot.title = element_text(color="black", size=30, face="bold"),
        legend.text = element_text(size = 20),
        text = element_text(size = 30),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
print(Prdm1)
ggsave("Prdm1.tiff", width = 30, height = 30, units = "cm", Prdm1, compression = "lzw")

patch.Bcells <- UMAP1 | CD19 | Prdm1
ggsave("patch.Bcells.tiff", width = 90, height = 30, units = "cm", patch.Bcells, compression = "lzw")
Bcell.clusters <- c("0", "1", "2", "10","11", "17", "18", "24", "25", "26", "29", "31") #12 cludsters

patch.Tcells <- CD4 | CD8a | Ncr1
ggsave("patch.Tcells.tiff", width = 90, height = 30, units = "cm", patch.Tcells, compression = "lzw")
Tcell.clusters <- c("3", "4", "5", "6", "7", "9", "14", "16", "22", "28", "30") #11 clusters

Other.cells <- c("8", "12", "13", "15", "19", "20", "21", "23", "27", "32", "33", "34") #12 clusters

#Final plots
patch.broad.clustering <- patch.Bcells / patch.Tcells
ggsave("patch.broad.clustering.tiff", width = 90, height = 60, units = "cm", patch.broad.clustering, compression = "lzw")


####Re-clustering B cells####
DefaultAssay(Allcells) <- "RNA"
Bcells <- subset(Allcells, idents = Bcell.clusters)

DefaultAssay(Bcells) <- "ADT"
rownames(Bcells)

DefaultAssay(Bcells) <- "RNA"
rownames(Bcells)

#Remove previous ssn_res columns#
head(Bcells[[]])
columns.to.remove <- c("RNA_snn_res.0.2", "ADT_snn_res.0.2", "wsnn_res.0.8")
for(i in columns.to.remove) {
  Bcells[[i]] <- NULL
}
head(Bcells[[]])

#RNA normalisation#
DefaultAssay(Bcells) <- "RNA" #For log normalisation
Bcells <- NormalizeData(Bcells, verbose = TRUE)
Bcells <- FindVariableFeatures(Bcells, nfeatures = 3000)
all.genes <- rownames(Bcells)
Bcells <- ScaleData(Bcells, features = all.genes) #Scaling based on genes

#Visualisation#
top20 <-  head(VariableFeatures(Bcells), 20)
plot1.1 <-  VariableFeaturePlot(Bcells)
top20_plot <-  LabelPoints(plot = plot1.1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)
print(top20_plot)
ggsave("top20_plot.tiff", width = 30, height = 20, units = "cm", top20_plot, compression = "lzw")


#RNA PCA#
Bcells <- RunPCA(Bcells, verbose = FALSE, features = VariableFeatures(object = Bcells))
pca_variance <- Bcells@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component",
     xlim = c(0, 50))
abline(h = 0.01) #37

#RNA clustering#
DefaultAssay(Bcells) <- "RNA" #For log normalisation

Bcells <- FindNeighbors(Bcells, dims = 1:37)
Bcells <- FindClusters(Bcells, resolution = 0.2, verbose = FALSE) #0.2 for the resolution
clustree.RNA <- clustree(Bcells, prefix = "RNA_snn_res.") + theme(legend.position="bottom")
ggsave("Bcells.clustree.RNA.tiff", width = 60, height = 40, units = "cm", clustree.RNA, compression = "lzw")

Bcells <- RunUMAP(Bcells, dims = 1:37, assay = 'RNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
Bcells_p1 <- DimPlot(Bcells, label = TRUE, reduction = "rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, cols = col_con2) +  ggtitle("RNA Clustering") + theme_bw() + NoLegend()
Bcells_p1 <- Bcells_p1 + theme(plot.title = element_text(color="black", size=25, face="bold"))
ggsave("Bcells_p1.tiff", width = 30, height = 30, units = "cm", Bcells_p1, compression = "lzw")


#ADT normalisation
DefaultAssay(Bcells) <- "ADT"
VariableFeatures(Bcells) <- rownames(Bcells[["ADT"]])
Bcells <- NormalizeData(Bcells, normalization.method = "CLR", margin = 2)
all.abs <- rownames(Bcells)
Bcells <- ScaleData(Bcells, features = all.abs) #Scaling for all abs, unmapped as well

#ADT PCA
Bcells <- RunPCA(Bcells, reduction.name = 'apca', approx = FALSE)
apca_variance <- Bcells@reductions$apca@stdev^2
plot(apca_variance/sum(apca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #27

#ADT clusterin
Bcells <- FindNeighbors(Bcells, dims = 1:27, reduction = "apca")
Bcells <- FindClusters(Bcells, resolution = 0.3, verbose = FALSE) #0.3 for the resolution
clustree.ADT <- clustree(Bcells, prefix = "ADT_snn_res.") + theme(legend.position = "bottom")
ggsave("Bcells.clustree.ADT.tiff", width = 60, height = 40, units = "cm", clustree.ADT, compression = "lzw")

Bcells <- RunUMAP(Bcells, reduction = 'apca', dims = 1:27, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
Bcells_p2 <- DimPlot(Bcells, label = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, cols = col_con2) +  ggtitle("ADT Clustering") + theme_bw() + NoLegend()
Bcells_p2 <- Bcells_p2 + theme(plot.title = element_text(color="black", size=25, face="bold"))
ggsave("Bcells_p2.tiff", width = 30, height = 30, units = "cm", Bcells_p2, compression = "lzw")

#WNN
DefaultAssay(Bcells) <- "RNA" #For log normalisation

#Combine into wnn plot
Bcells <- FindMultiModalNeighbors(
  Bcells, reduction.list = list("pca", "apca"), 
  dims.list = list(1:37, 1:27), modality.weight.name = "RNA.weight")

#WNN clustering
Bcells <- FindClusters(Bcells, graph.name = "wsnn", algorithm = 3, resolution = 0.2, verbose = TRUE) #0.2 for the resolution
clustree.wnn <- clustree(Bcells, prefix = "wsnn_res.") + theme(legend.position="bottom")
ggsave("Bcells.clustree.wnn.tiff", width = 60, height = 40, units = "cm", clustree.wnn, compression = "lzw")

Bcells <- RunUMAP(Bcells, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Bcells_p3 <- DimPlot(Bcells, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, cols = col_con2) +  ggtitle("WNN Clusters") + theme_bw() + NoLegend()
Bcells_p3 <- Bcells_p3 + theme(plot.title = element_text(color="black", size = 25, face = "bold")) + xlab("UMAP1") + ylab("UMAP2")
ggsave("Bcells_p3.tiff", width = 30, height = 30, units = "cm", Bcells_p3, compression = "lzw")

#Final UMAP plots
patch.Bcell.dimplots <- Bcells_p1 | Bcells_p2 | Bcells_p3
ggsave("patch.Bcell.dimplots.tiff", width = 90, height = 30, units = "cm", patch.Bcell.dimplots, compression = "lzw")

#Save SeuratObejct
SaveH5Seurat(Bcells, "CITESeq1_Bcells_nBCRnNA.h5seurat", overwrite = TRUE)


####Load B cells Seurat object####
Bcells <- LoadH5Seurat(file.choose("CITESeq1_Bcells_nBCRnNA.h5seurat"))#With clonotype data, excluding NAs, no BCR genes (except for class inforamtion)
Bcells$seurat_clusters <- Bcells$wsnn_res.0.2
Idents(Bcells) <- Bcells$seurat_clusters

my_levels <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
               "11", "12")
Idents(Bcells) <- factor(Idents(Bcells), levels = my_levels)
Bcells$seurat_clusters <- factor(Idents(Bcells), levels = my_levels)

