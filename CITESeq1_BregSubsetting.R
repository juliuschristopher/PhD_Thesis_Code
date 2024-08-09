## ##Sub-clustering of Cluster 2 (B regulatory cells) CITE-Seq (1) dataset####
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
library(msigdbr)
library(presto)
library(fgsea)
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
sessionInfo()

##Colour paneles##
col_con1 <- createPalette(50,  c("#2A9D8F", "#E9C46A", "#E76F51"))
col_con1 <-as.character(col_con1)
col_con2 <- c("#C94055", "#EFE3A2", "#4FA4C2", "#E80DFF", "#26FF0D", "#0D65FB", "#0DF6BD", "#F8000D", "#FE0DBD", "#FDA30D",
              "#FCCFEB", "#B4F2A7", "#C5F00D", "#FC6A00", "#F8AA92", "#94B2F4", "#FE8BBD", "#F24798", "#22FAFB", "#B898FB") #Adjusted for 20x bright colours
col_con3 <- createPalette(50,  c("#CE4257", "#FFF2B0", "#4BA3C3"))
col_con3 <-as.character(col_con2)
col_con6 <- c("#B4F2A7", "#E6C266", "#BFEC00", "#FF00FC", "#E76F51", "skyblue")
col_con.genotypes <- c("#666666", "#fdcc66", "#0f80ff", "#91bd95") #For Genotypes
col_con.mouse <- c("grey1", "azure4", "gold1", "chocolate2", "blue", "forestgreen", "darkolivegreen1")#For individual mice
col_con.sex <- c("darkred", "darkturquoise") #For Sex
col_con_clonotype <- rev(c("grey39", "#4FA4C2", "chartreuse4", "gold2", "darkturquoise","firebrick3")) #Clonotype all
col_con.clonotype.no.hyper <- c( "darkturquoise", "gold2", "chartreuse4", "grey39") #Clonotype no hyperexpansion


##Vector orders##
clo_order <- c("NA", "Small (1e-04 < X <= 0.001)", "Medium (0.001 < X <= 0.01)", "Large (0.01 < X <= 0.1)", "Hyperexpanded (0.1 < X <= 1)")
clo_order.reverse <- c("Hyperexpanded (0.1 < X <= 1)", "Large (0.01 < X <= 0.1)", "Medium (0.001 < X <= 0.01)", "Small (1e-04 < X <= 0.001)","NA")
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

#Load Pathway analysis datasets
msigdbr_species() #Show species
Hallmark <- msigdbr(species="Mus musculus", category="H") #Hallmark genes
C1 <- msigdbr(species="Mus musculus", category="C1") #Positional gene sets - corresponding to human chromosome cytogenetic bands
C2 <- msigdbr(species="Mus musculus", category="C2") #Curated gene sets, including chemical and genetic perturbations (CGP) and canonical pathways (CP) (BioCarta, KEGG, PID, Reactome and WikiPathways)
C2_CP_WIKIPATHWAYS <- msigdbr(species="Mus musculus", category="C2", subcategory = "CP:WIKIPATHWAYS")
C2_CP_REACTOME <- msigdbr(species="Mus musculus", category="C2", subcategory = "CP:REACTOME")
C2_CP_PID <- msigdbr(species="Mus musculus", category="C2", subcategory = "CP:PID")
C2_CP_KEGG <- msigdbr(species="Mus musculus", category="C2", subcategory = "CP:KEGG")
C2_CP_BIOCARTA <- msigdbr(species="Mus musculus", category="C2", subcategory = "CP:BIOCARTA")
C3 <- msigdbr(species="Mus musculus", category="C3") #Regulatory target gene sets
C3_TFT_GTRD <- msigdbr(species="Mus musculus", category="C3", subcategory = "TFT:GTRD") #Regulatory target gene sets Transcription Factor Targets - GTRD
C3_TFT_LEGACY <- msigdbr(species="Mus musculus", category="C3", subcategory = "TFT:TFT_Legacy") #Regulatory target gene sets Transcription Factor Targets - Legacy
C4 <- msigdbr(species="Mus musculus", category="C4") #Computational gene sets
C5 <- msigdbr(species="Mus musculus", category="C5") #Ontology gene sets
C5_BP <- msigdbr(species="Mus musculus", category="C5", subcategory = "GO:BP") #GO Biological process
C5_CC <- msigdbr(species="Mus musculus", category="C5", subcategory = "GO:CC") #GO Cellular component
C5_MF <- msigdbr(species="Mus musculus", category="C5", subcategory = "GO:MF") #GO Molecular function
C6 <- msigdbr(species="Mus musculus", category="C6") #Oncogenicsignature gene sets
C7_IMMUNESIGDB <- msigdbr(species="Mus musculus", category="C7", subcategory = "IMMUNESIGDB") #Immunologic signature gene sets - ImmuneSigDB gene sets - 4872 sets
C8 <- msigdbr(species="Mus musculus", category="C8") #Cell type signature gene sets


####Load B cell Suerat Object####
Bcells <- LoadH5Seurat(file.choose("CITESeq1_Bcells_nBCRnNA.h5seurat"))#With clonotype data, excluding NAs, no BCR genes (except for class inforamtion)
Bcells$seurat_clusters <- Bcells$wsnn_res.0.2
Idents(Bcells) <- Bcells$seurat_clusters

my_levels <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
               "11", "12")
Idents(Bcells) <- factor(Idents(Bcells), levels = my_levels)
Bcells$seurat_clusters <- factor(Idents(Bcells), levels = my_levels)

names(col_con2) <- levels(Bcells$seurat_clusters) #One colour associated with one cluster
UMAP2 <- DimPlot(Bcells, reduction = "wnn.umap", cols = col_con2, pt.size = 2, label = TRUE, label.box = TRUE, label.size = 8) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("B cells") +
  theme(plot.title = element_text(color="black", size= 30, face="bold"),
        legend.text = element_text(size = 20),
        text = element_text(size = 20)) +
  NoLegend()
print(UMAP2)
ggsave("UMAP2.tiff", width = 30, height = 30, units = "cm", UMAP2, compression = "lzw")


####Labeling of B cell Clusters####
Idents(Bcells) <- Bcells$seurat_clusters
Bcells <- RenameIdents(Bcells,
                       `0` = "Follicular B cells", `1` = "Marginal Zone B cells", `2` = "B regulatory cells", `3` = "High Mito. genes",
                       `4` = "Germinal Center Light Zone B cells", `5` = "Atypical B cells", `6` = "Transitional B cells", `7` = "High Ribo. genes",
                       `8` = "CTLA4+ B cells", `9` = "Plasma Cells/Plasmablasts", `10` = "IFN-responsive B cells", `11` = "Outlier (1)", `12` = "Outlier (2)")
Bcells[["Label"]] <- Idents(Bcells)
Idents(Bcells) <- Bcells$Label

names(col_con2) <- levels(Bcells$Label) #One colour associated with one cluster
UMAP3 <- DimPlot(Bcells, reduction = "wnn.umap", cols = col_con2, pt.size = 2, label = TRUE, label.box = TRUE, label.size = 6 , repel = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("B cells") +
  theme(plot.title = element_text(color="black", size= 30, face="bold"),
        legend.text = element_text(size = 20),
        text = element_text(size = 20)) +
  NoLegend()
print(UMAP3)
ggsave("UMAP3.tiff", width = 30, height = 30, units = "cm", UMAP3, compression = "lzw")


####Remove Outlier Clusters####
Bcells_all_clusters <- Bcells
Bcells <- subset(Bcells, idents = c("Outlier (1)", "Outlier (2)"), invert = TRUE)
Bcells$Label

UMAP4 <- DimPlot(Bcells, reduction = "wnn.umap", cols = col_con2, pt.size = 2, label = TRUE, label.box = TRUE, label.size = 6 , repel = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("B cells") +
  theme(plot.title = element_text(color="black", size= 30, face="bold"),
        legend.text = element_text(size = 20),
        text = element_text(size = 20)) +
  NoLegend()
print(UMAP4)
ggsave("UMAP4.tiff", width = 30, height = 30, units = "cm", UMAP4, compression = "lzw")


####Re-clustering B cells####
DefaultAssay(Bcells) <- "RNA"
Idents(Bcells) <- Bcells$seurat_clusters
UMAP2
Cluster2 <- subset(Bcells, idents = "2")

DefaultAssay(Cluster2) <- "ADT"
rownames(Cluster2)

DefaultAssay(Cluster2) <- "RNA"
rownames(Cluster2)

#Remove previous ssn_res columns#
head(Cluster2[[]])
columns.to.remove <- c("RNA_snn_res.0.2", "ADT_snn_res.0.3", "wsnn_res.0.2")
for(i in columns.to.remove) {
  Cluster2[[i]] <- NULL
}
head(Cluster2[[]])

#RNA normalisation#
DefaultAssay(Cluster2) <- "RNA" #For log normalisation
Cluster2 <- NormalizeData(Cluster2, verbose = TRUE)
Cluster2 <- FindVariableFeatures(Cluster2, nfeatures = 3000)
all.genes <- rownames(Cluster2)
Cluster2 <- ScaleData(Cluster2, features = all.genes) #Scaling based on genes

#Visualisation#
top20 <-  head(VariableFeatures(Cluster2), 20)
plot1.1 <-  VariableFeaturePlot(Cluster2)
top20_plot <-  LabelPoints(plot = plot1.1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)
print(top20_plot)
ggsave("top20_plot.tiff", width = 30, height = 20, units = "cm", top20_plot, compression = "lzw")


#RNA PCA#
Cluster2 <- RunPCA(Cluster2, verbose = FALSE, features = VariableFeatures(object = Cluster2))
pca_variance <- Cluster2@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component",
     xlim = c(0, 100))
abline(h = 0.01) #50

#RNA clustering#
DefaultAssay(Cluster2) <- "RNA" #For log normalisation

Cluster2 <- FindNeighbors(Cluster2, dims = 1:50)
Cluster2 <- FindClusters(Cluster2, resolution = 1.2, verbose = FALSE) #1.2 for the resolution
clustree.RNA <- clustree(Cluster2, prefix = "RNA_snn_res.") + theme(legend.position="bottom")
ggsave("Cluster2.Bcells.clustree.RNA.tiff", width = 60, height = 40, units = "cm", clustree.RNA, compression = "lzw")

Cluster2 <- RunUMAP(Cluster2, dims = 1:50, assay = 'RNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
Cluster2_p1 <- DimPlot(Cluster2, label = TRUE, reduction = "rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, cols = col_con2) +  ggtitle("RNA Clustering") + theme_bw() + NoLegend()
Cluster2_p1 <- Cluster2_p1 + theme(plot.title = element_text(color="black", size=25, face="bold"))
ggsave("Cluster2_p1.tiff", width = 30, height = 30, units = "cm", Cluster2_p1, compression = "lzw")


#ADT normalisation
DefaultAssay(Cluster2) <- "ADT"
VariableFeatures(Cluster2) <- rownames(Cluster2[["ADT"]])
Cluster2 <- NormalizeData(Cluster2, normalization.method = "CLR", margin = 2)
all.abs <- rownames(Cluster2)
Cluster2 <- ScaleData(Cluster2, features = all.abs) #Scaling for all abs, unmapped as well

#ADT PCA
Cluster2 <- RunPCA(Cluster2, reduction.name = 'apca', approx = FALSE)
apca_variance <- Cluster2@reductions$apca@stdev^2
plot(apca_variance/sum(apca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #27

#ADT clusterin
Cluster2 <- FindNeighbors(Cluster2, dims = 1:27, reduction = "apca")
Cluster2 <- FindClusters(Cluster2, resolution = 1.0, verbose = FALSE) #1.0 for the resolution
clustree.ADT <- clustree(Cluster2, prefix = "ADT_snn_res.") + theme(legend.position = "bottom")
ggsave("Cluster2.Bcells.clustree.ADT.tiff", width = 60, height = 40, units = "cm", clustree.ADT, compression = "lzw")

Cluster2 <- RunUMAP(Cluster2, reduction = 'apca', dims = 1:27, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
Cluster2_p2 <- DimPlot(Cluster2, label = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, cols = col_con2) +  ggtitle("ADT Clustering") + theme_bw() + NoLegend()
Cluster2_p2 <- Cluster2_p2 + theme(plot.title = element_text(color="black", size=25, face="bold"))
ggsave("Cluster2s_p2.tiff", width = 30, height = 30, units = "cm", Cluster2_p2, compression = "lzw")

#WNN
DefaultAssay(Cluster2) <- "RNA" #For log normalisation

#Combine into wnn plot
Cluster2 <- FindMultiModalNeighbors(
  Cluster2, reduction.list = list("pca", "apca"), 
  dims.list = list(1:50, 1:27), modality.weight.name = "RNA.weight")

#WNN clustering
Cluster2 <- FindClusters(Cluster2, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = TRUE) #0.5 for the resolution
clustree.wnn <- clustree(Cluster2, prefix = "wsnn_res.") + theme(legend.position="bottom")
ggsave("Cluster2.Bcells.clustree.wnn.tiff", width = 60, height = 40, units = "cm", clustree.wnn, compression = "lzw")

Cluster2 <- RunUMAP(Cluster2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Cluster2_p3 <- DimPlot(Cluster2, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, cols = col_con2) +  ggtitle("WNN Clusters") + theme_bw() + NoLegend()
Cluster2_p3 <- Cluster2_p3 + theme(plot.title = element_text(color="black", size = 25, face = "bold")) + xlab("UMAP1") + ylab("UMAP2")
ggsave("Cluster2_p3.tiff", width = 30, height = 30, units = "cm", Cluster2_p3, compression = "lzw")

#Final UMAP plots
patch.Cluster2.dimplots <- Cluster2_p1 | Cluster2_p2 | Cluster2_p3
ggsave("patch.Cluster2.dimplots.tiff", width = 90, height = 30, units = "cm", patch.Cluster2.dimplots, compression = "lzw")

#Save SeuratObejct
SaveH5Seurat(Cluster2, "CITESeq1_Cluster2_Bcells_nBCRnNA.h5seurat", overwrite = TRUE)

####Load Cluster 2 Seurat Object####
Cluster2 <- LoadH5Seurat(file.choose("CITESeq1_Cluster2_Bcells_nBCRnNA.h5seurat")) #With clonotype data, excluding NAs, no BCR genes (except for class inforamtion)
Cluster2$seurat_clusters <- Cluster2$wsnn_res.0.5
Idents(Cluster2) <- Cluster2$seurat_clusters

