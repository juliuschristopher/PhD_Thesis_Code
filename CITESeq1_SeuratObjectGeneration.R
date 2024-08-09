## ##CITE-Seq (1) Script - Seurat Object Generations####
####Setup####
#Set working directory to folder containing files
#Load required packages
library(Polychrome)
library(Seurat)
library(future)
library(RColorBrewer)
library(viridis)
library(stringr)
library(Matrix)
library(ggplot2)
library(SeuratDisk)
library(clustree)
library(scRepertoire)
library(patchwork)

#Colour palettes
col_con1 <- createPalette(50,  c("#2A9D8F", "#E9C46A", "#E76F51"))
col_con1 <-as.character(col_con1)

#Alter working capacity
options(future.globals.maxSize= 10097152000) # 10Gb

####Load the 10X Cell Ranger output####
#Note: CITE-Seq (1) files are within the CITE-Seq2 folder, due to different ordering compared to the chapters in the PhD Thesis
#Read the 10x Cell Ranger Output
a_ge.data <- Read10X(data.dir = "~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/GE/A_WT_GE/outs/filtered_feature_bc_matrix")
b_ge.data <- Read10X(data.dir = "~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/GE/B_WT_GE/outs/filtered_feature_bc_matrix")
c_ge.data <- Read10X(data.dir = "~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/GE/C_BCL6_GE/outs/filtered_feature_bc_matrix")
d_ge.data <- Read10X(data.dir = "~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/GE/D_BCL6_GE/outs/filtered_feature_bc_matrix")
f_ge.data <- Read10X(data.dir = "~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/GE/F_E1020K_GE/outs/filtered_feature_bc_matrix")
g_ge.data <- Read10X(data.dir = "~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/GE/G_E1020K_BCL6_GE/outs/filtered_feature_bc_matrix")
h_ge.data <- Read10X(data.dir = "~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/GE/H_E1020K_BCL6_GE/outs/filtered_feature_bc_matrix")

#Add the sample to the cell names, consistent with antibody data below
colnames(a_ge.data)=gsub("-1","_a",colnames(a_ge.data))
colnames(b_ge.data)=gsub("-1","_b",colnames(b_ge.data))
colnames(c_ge.data)=gsub("-1","_c",colnames(c_ge.data))
colnames(d_ge.data)=gsub("-1","_d",colnames(d_ge.data))
colnames(f_ge.data)=gsub("-1","_f",colnames(f_ge.data))
colnames(g_ge.data)=gsub("-1","_g",colnames(g_ge.data))
colnames(h_ge.data)=gsub("-1","_h",colnames(h_ge.data))

#Capitalise with lower case letters (genes)
rownames(a_ge.data)=str_to_title(rownames(a_ge.data))
rownames(b_ge.data)=str_to_title(rownames(b_ge.data))
rownames(c_ge.data)=str_to_title(rownames(c_ge.data))
rownames(d_ge.data)=str_to_title(rownames(d_ge.data))
rownames(f_ge.data)=str_to_title(rownames(f_ge.data))
rownames(g_ge.data)=str_to_title(rownames(g_ge.data))
rownames(h_ge.data)=str_to_title(rownames(h_ge.data))
head(f_ge.data)

####Load 10X Antibody data####
#Read the 10x Antibody output
a_ab.data <- Read10X(data.dir = "~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/A_WT/umi_count",gene.column=1)
b_ab.data <- Read10X(data.dir = "~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/B_WT/umi_count",gene.column=1)
c_ab.data <- Read10X(data.dir = "~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/C_BCL6/umi_count",gene.column=1)
d_ab.data <- Read10X(data.dir = "~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/D_BCL6/umi_count",gene.column=1)
f_ab.data <- Read10X(data.dir = "~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/F_E1020K/umi_count",gene.column=1)
g_ab.data <- Read10X(data.dir = "~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/G_E1020K_BCL6/umi_count",gene.column=1)
h_ab.data <- Read10X(data.dir = "~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/H_E1020K_BCL6/umi_count",gene.column=1)


#Tidy up the rownames from the data
rownames(a_ab.data)=gsub("-[^-]+$","",rownames(a_ab.data),perl=TRUE)
rownames(b_ab.data)=gsub("-[^-]+$","",rownames(b_ab.data),perl=TRUE)
rownames(c_ab.data)=gsub("-[^-]+$","",rownames(c_ab.data),perl=TRUE)
rownames(d_ab.data)=gsub("-[^-]+$","",rownames(d_ab.data),perl=TRUE)
rownames(f_ab.data)=gsub("-[^-]+$","",rownames(f_ab.data),perl=TRUE)
rownames(g_ab.data)=gsub("-[^-]+$","",rownames(g_ab.data),perl=TRUE)
rownames(h_ab.data)=gsub("-[^-]+$","",rownames(h_ab.data),perl=TRUE)

#Capitalise with lower case letters (antobodies)
rownames(a_ab.data)=str_to_title(rownames(a_ab.data))
rownames(b_ab.data)=str_to_title(rownames(b_ab.data))
rownames(c_ab.data)=str_to_title(rownames(c_ab.data))
rownames(d_ab.data)=str_to_title(rownames(d_ab.data))
rownames(f_ab.data)=str_to_title(rownames(f_ab.data))
rownames(g_ab.data)=str_to_title(rownames(g_ab.data))
rownames(h_ab.data)=str_to_title(rownames(h_ab.data))
head(f_ab.data)

#Add the Sample to the cell names in each sample
colnames(a_ab.data)=paste(colnames(a_ab.data),"_a",sep="")
colnames(b_ab.data)=paste(colnames(b_ab.data),"_b",sep="")
colnames(c_ab.data)=paste(colnames(c_ab.data),"_c",sep="")
colnames(d_ab.data)=paste(colnames(d_ab.data),"_d",sep="")
colnames(f_ab.data)=paste(colnames(f_ab.data),"_f",sep="")
colnames(g_ab.data)=paste(colnames(g_ab.data),"_g",sep="")
colnames(h_ab.data)=paste(colnames(h_ab.data),"_h",sep="")
head(f_ab.data)

####Combine 10X Cell Ranger and Antibody Data into a Suerat Object####
m <- Matrix(nrow = nrow(a_ab.data), ncol = ncol(a_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(a_ab.data)
colnames(m)=colnames(a_ge.data)
common=intersect(colnames(a_ge.data),colnames(a_ab.data))
m[,common]=a_ab.data[,common]
a = CreateSeuratObject(counts = a_ge.data,project="a", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
a[["ADT"]] <- adt_assay

ncol(a_ge.data) #4609
ncol(a_ab.data) #4607
ncol(m) #4609


m <- Matrix(nrow = nrow(b_ab.data), ncol = ncol(b_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(b_ab.data)
colnames(m)=colnames(b_ge.data)
common=intersect(colnames(b_ge.data),colnames(b_ab.data))
m[,common]=b_ab.data[,common]
b = CreateSeuratObject(counts = b_ge.data,project="b", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
b[["ADT"]] <- adt_assay

ncol(b_ge.data) #5117
ncol(b_ab.data) #5115
ncol(m) #5117


m <- Matrix(nrow = nrow(c_ab.data), ncol = ncol(c_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(c_ab.data)
colnames(m)=colnames(c_ge.data)
common=intersect(colnames(c_ge.data),colnames(c_ab.data))
m[,common]=c_ab.data[,common]
c = CreateSeuratObject(counts = c_ge.data,project="c", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
c[["ADT"]] <- adt_assay

ncol(c_ge.data) #4818
ncol(c_ab.data) #4817
ncol(m) #4818


m <- Matrix(nrow = nrow(d_ab.data), ncol = ncol(d_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(d_ab.data)
colnames(m)=colnames(d_ge.data)
common=intersect(colnames(d_ge.data),colnames(d_ab.data))
m[,common]=d_ab.data[,common]
d = CreateSeuratObject(counts = d_ge.data,project="d", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
d[["ADT"]] <- adt_assay

ncol(d_ge.data) #5459
ncol(d_ab.data) #5458
ncol(m) #5459


m <- Matrix(nrow = nrow(f_ab.data), ncol = ncol(f_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(f_ab.data)
colnames(m)=colnames(f_ge.data)
common=intersect(colnames(f_ge.data),colnames(f_ab.data))
m[,common]=f_ab.data[,common]
f = CreateSeuratObject(counts = f_ge.data,project="f", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
f[["ADT"]] <- adt_assay

ncol(f_ge.data) #11609
ncol(f_ab.data) #11604
ncol(m) #11609


m <- Matrix(nrow = nrow(g_ab.data), ncol = ncol(g_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(g_ab.data)
colnames(m)=colnames(g_ge.data)
common=intersect(colnames(g_ge.data),colnames(g_ab.data))
m[,common]=g_ab.data[,common]
g = CreateSeuratObject(counts = g_ge.data,project="g", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
g[["ADT"]] <- adt_assay

ncol(g_ge.data) #6524
ncol(g_ab.data) #6521
ncol(m) #6524


m <- Matrix(nrow = nrow(h_ab.data), ncol = ncol(h_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(h_ab.data)
colnames(m)=colnames(h_ge.data)
common=intersect(colnames(h_ge.data),colnames(h_ab.data))
m[,common]=h_ab.data[,common]
h = CreateSeuratObject(counts = h_ge.data,project="h", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
h[["ADT"]] <- adt_assay

ncol(h_ge.data) #5879
ncol(h_ab.data) #5876
ncol(m) #5879


####Merge suerat objects####
experiments=c(a,b,c,d,f,g,h)
experiment_names=c("a","b","c","d","f","g","h")
experiment<-merge(x= a, y=c(b,c,d,f,g,h))
head(experiment[[]]) #21760 features across 44015 samples within 2 assays

####Quality control: Filtering####
#Mitochondrial QC metrics
experiment[["percent.Mt"]] <- PercentageFeatureSet(experiment, pattern = "^Mt-")

#Ribosomal QC metric
experiment[["percent.Ribo"]] <- PercentageFeatureSet(experiment, pattern = "^Rp[sl]")

#Remove where nCount_ADT = 0
DefaultAssay(experiment) <- "ADT"
experiment_noADT <- subset(experiment, nCount_ADT == 0) #9244 cells with no antibodies
experiment <- subset(experiment, nCount_ADT > 0) #Removing cells with no antibodies bound - 21760 features across 34771 after
DefaultAssay(experiment) <- "RNA"

#Visualize QC metrics as violin plot
RNA_QC <- VlnPlot(experiment, features = c("nFeature_RNA", "nCount_RNA", "percent.Mt", "percent.Ribo"))
ADT_QC_1 <- VlnPlot(experiment, features = "nFeature_ADT")
ADT_QC_2 <- VlnPlot(experiment, features = "nCount_ADT", y.max = 10000)

#FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
feature_count.RNA_vs_percent.Mt = FeatureScatter(experiment, feature1 = "nCount_RNA", feature2 = "percent.Mt") + NoLegend()  +
  ylab("% of mitochondrial genes") +
  xlab("UMI counts") + 
  geom_hline(yintercept = 5) 
feature_count.RNA_vs_feature.RNA = FeatureScatter(experiment, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() +
  ylab("Number of genes") +
  xlab("UMI counts") + 
  geom_hline(yintercept = 200) 
feature_count.ADT_vs_feature.ADT = FeatureScatter(experiment, feature1 = "nCount_ADT", feature2 = "nFeature_ADT", jitter = TRUE) + NoLegend() +
  ylab("Number of different antibodies") + ylim(0, 30) +
  xlab("ADT counts") + xlim(0, 40000)

#Generally aim to filter out unique feature counts over 2,500 and less than 200, and cells with percent mitochondrial genes above 5
filter_seurat = function(seurat_object){
  
  message("Performing filter by number of genes.")
  seurat_object = subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
  message("Now the object has ", dim(seurat_object)[1], " genes and ", dim(seurat_object)[2], " cells.")
  
  
  return(seurat_object)
}

experiment <-  filter_seurat(experiment)
head(experiment[[]]) #Now the object has 21733 genes and 32813 cells.

####Removeing BCR and TCR genes####
#Removing BCR genes
nrow(experiment) #21733 genes

DefaultAssay(experiment) <- "ADT"
ADT.vector <- rownames(experiment) #27 antibodies
DefaultAssay(experiment) <- "RNA"

noBCRcounts <- GetAssayData(experiment, assay = "RNA") #21733 x 32813
noBCRcounts <- noBCRcounts[!grepl("^Ig[hkl]v", rownames(noBCRcounts)), ] #21495 x 32813
noBCRcounts <- noBCRcounts[!grepl("^Ig[hkl]j", rownames(noBCRcounts)), ] #21495 x 32813
noBCRcounts <- noBCRcounts[!grepl("^Ig[kl]c", rownames(noBCRcounts)), ] #21490 x 32813
#Optional - removing Antibody class information: noBCRcounts <- noBCRcounts[!grepl("^Igh[adegm]", rownames(noBCRcounts)), ]

experiment <- subset(experiment, features = c(rownames(noBCRcounts), ADT.vector))

DefaultAssay(experiment) <- "RNA"
nrow(experiment) #21490 genes

DefaultAssay(experiment) <- "ADT"
rownames(experiment) #27 antibodies
DefaultAssay(experiment) <- "RNA"

#Removing TCR genes
noTCRcounts <- GetAssayData(experiment, assay = "RNA") #21490 x 32813
noTCRcounts <- noTCRcounts[!grepl("^Tr[ab][vjc]", rownames(noTCRcounts)), ] #21385 x 32813
#Optional - removing gdT cell information: noTCRcounts <- noTCRcounts[!grepl("^Tr[dg][vjc]", rownames(noTCRcounts)), ]

experiment <- subset(experiment, features = c(rownames(noTCRcounts), ADT.vector))

DefaultAssay(experiment) <- "RNA"
nrow(experiment) #21385 genes

DefaultAssay(experiment) <- "ADT"
rownames(experiment) #27 antibodies
DefaultAssay(experiment) <- "RNA"


####Add Metadata to Seurat Object####
head(experiment[[]])

#Add mouse and genotype information
Idents(experiment) <- experiment$orig.ident
experiment<- RenameIdents(experiment, `a` = "WT_1", `b` ="WT_2", `c` ="BCL6_1", `d` = "BCL6_2", `f` = "E1020K_1",`g` ="E1020K_BCL6_1", `h` = "E1020K_BCL6_2")
experiment[["Mouse"]] <- Idents(experiment)

Idents(experiment) <- experiment$orig.ident
experiment <- RenameIdents(experiment, `a` = "WT", `b` ="WT", `c` ="BCL6", `d` = "BCL6", `f` = "E1020K",`g` ="E1020K_BCL6", `h` = "E1020K_BCL6")
experiment[["Genotype"]] <- Idents(experiment)

Idents(experiment) <- experiment$orig.ident
experiment <- RenameIdents(experiment, `a` = "Female", `b` ="Male", `c` ="Female", `d` = "Male", `f` = "Male",`g` ="Male", `h` = "Male")
experiment[["Sex"]] <- Idents(experiment)

Idents(experiment) <- experiment$Mouse

####Normalisation, Scaling, Dimensionality Reduction and Clustering####
#RNA normalisation
DefaultAssay(experiment) <- "RNA" #For log normalisation

experiment <- NormalizeData(experiment, verbose = TRUE)
experiment <- FindVariableFeatures(experiment, nfeatures = 3000)
all.genes <- rownames(experiment)
experiment <- ScaleData(experiment, features = all.genes)

#Visualisation of top 20 most differentially expressed genes
top20 <-  head(VariableFeatures(experiment), 20)
plot1.1 <-  VariableFeaturePlot(experiment)
top20_plot <-  LabelPoints(plot = plot1.1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)
ggsave("top20_plot.RNA.tiff", width = 30, height = 20, units = "cm", top20_plot, compression = "lzw")

#RNA PCA
experiment <- RunPCA(experiment, verbose = FALSE, features = VariableFeatures(object = experiment))
pca_variance <- experiment@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #27

#RNA clustering
DefaultAssay(experiment) <- "RNA" #Ensure correct assay is selected

experiment <- FindNeighbors(experiment, dims = 1:27)
experiment <- FindClusters(experiment, resolution = 0.2, verbose = FALSE) #0.2 for the resolution
clustree.RNA <- clustree(experiment, prefix = "RNA_snn_res.") + theme(legend.position="bottom")
ggsave("clustree.RNA.tiff", width = 60, height = 40, units = "cm", clustree.RNA, compression = "lzw")

experiment <- RunUMAP(experiment, dims = 1:27, assay = 'RNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
experiment_p1 <- DimPlot(experiment, label = TRUE, reduction = "rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, cols = col_con1) +  ggtitle("RNA Clustering") + theme_bw() + NoLegend()
experiment_p1 <- experiment_p1 + theme(plot.title = element_text(color="black", size=25, face="bold"))
ggsave("experiment_p1.tiff", width = 30, height = 30, units = "cm", experiment_p1, compression = "lzw")


#ADT normalisation#
DefaultAssay(experiment) <- "ADT"
VariableFeatures(experiment) <- rownames(experiment[["ADT"]])
experiment <- NormalizeData(experiment, normalization.method = "CLR", margin = 2)
all.abs <- rownames(experiment)
experiment <- ScaleData(experiment, features = all.abs) #Scaling for all abs, unmapped as well

#ADT PCA#
experiment <- RunPCA(experiment, reduction.name = 'apca', approx = FALSE)
apca_variance <- experiment@reductions$apca@stdev^2
plot(apca_variance/sum(apca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #25

#ADT clustering
experiment <- FindNeighbors(experiment, dims = 1:25, reduction = "apca")
experiment <- FindClusters(experiment, resolution = 0.2, verbose = FALSE) #0.2 for the resolution
clustree.ADT <- clustree(experiment, prefix = "ADT_snn_res.") + theme(legend.position="bottom")
ggsave("clustree.ADT.tiff", width = 60, height = 40, units = "cm", clustree.ADT, compression = "lzw")

experiment <- RunUMAP(experiment, reduction = 'apca', dims = 1:25, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
experiment_p2 <- DimPlot(experiment, label = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, cols = col_con1) +  ggtitle("ADT Clustering") + theme_bw() + NoLegend()
experiment_p2 <- experiment_p2 + theme(plot.title = element_text(color="black", size=25, face="bold"))
ggsave("experiment_p2.tiff", width = 30, height = 30, units = "cm", experiment_p2, compression = "lzw")


#WNN#
DefaultAssay(experiment) <- "RNA" #For log normalisation

#Combine into wnn plot#
experiment <- FindMultiModalNeighbors(
  experiment, reduction.list = list("pca", "apca"), 
  dims.list = list(1:27, 1:25), modality.weight.name = "RNA.weight")

#WNN clustering#
experiment <- FindClusters(experiment, graph.name = "wsnn", algorithm = 3, resolution = 0.8, verbose = TRUE) #0.8 for the resolution
clustree.wnn <- clustree(experiment, prefix = "wsnn_res.") + theme(legend.position="bottom")
ggsave("clustree.wnn.tiff", width = 60, height = 40, units = "cm", clustree.wnn, compression = "lzw")

experiment <- RunUMAP(experiment, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
experiment_p3 <- DimPlot(experiment, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, cols = col_con1) +  ggtitle("Seurat Clusters") + theme_bw() + NoLegend()
experiment_p3 <- experiment_p3 + theme(plot.title = element_text(color="black", size = 25, face = "bold")) + xlab("UMAP1") + ylab("UMAP2")

patch.dimplots <- experiment_p1 | experiment_p2 | experiment_p3
ggsave("patch.dimplots.tiff", width = 90, height = 30, units = "cm", patch.dimplots, compression = "lzw")


####Loading VDJ data seperately####
##Load contig file##
a_cl.data <- read.csv("~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/A_WT_VDJ/outs/filtered_contig_annotations.csv")
b_cl.data <- read.csv("~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/B_WT_VDJ/outs/filtered_contig_annotations.csv")
c_cl.data <- read.csv("~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/C_BCL6_VDJ/outs/filtered_contig_annotations.csv")
d_cl.data <- read.csv("~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/D_BCL6_VDJ/outs/filtered_contig_annotations.csv")
f_cl.data <- read.csv("~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/F_E1020K_VDJ/outs/filtered_contig_annotations.csv")
g_cl.data <- read.csv("~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/G_E1020K_BCL6_VDJ/outs/filtered_contig_annotations.csv")
h_cl.data <- read.csv("~/Desktop/Datasets/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/H_E1020K_BCL6_VDJ/outs/filtered_contig_annotations.csv")

##Match barcode names with GE and ADT data##
a_cl.data$barcode=gsub("-1","_a",a_cl.data$barcode)
b_cl.data$barcode=gsub("-1","_b",b_cl.data$barcode)
c_cl.data$barcode=gsub("-1","_c",c_cl.data$barcode)
d_cl.data$barcode=gsub("-1","_d",d_cl.data$barcode)
f_cl.data$barcode=gsub("-1","_f",f_cl.data$barcode)
g_cl.data$barcode=gsub("-1","_g",g_cl.data$barcode)
h_cl.data$barcode=gsub("-1","_h",h_cl.data$barcode)
colnames(h_cl.data)

contig_list <- list(a_cl.data, b_cl.data, c_cl.data, d_cl.data, f_cl.data, g_cl.data, h_cl.data)
head(contig_list)

#Generate combined object
combined <- combineBCR(contig_list, samples = c("a", "b", "c", "d", "f", "g", "h"))
combined2 <- combineBCR(contig_list, samples = c("a", "b", "c", "d", "f", "g", "h"), removeNA = TRUE)

#Make sure barcodes are identical to GE and ADT data
combined$a$barcode=gsub("a_","",combined$a$barcode)
combined$b$barcode=gsub("b_","",combined$b$barcode)
combined$c$barcode=gsub("c_","",combined$c$barcode)
combined$d$barcode=gsub("d_","",combined$d$barcode)
combined$f$barcode=gsub("f_","",combined$f$barcode)
combined$g$barcode=gsub("g_","",combined$g$barcode)
combined$h$barcode=gsub("h_","",combined$h$barcode)
head(combined$h)

combined2$a$barcode=gsub("a_","",combined2$a$barcode)
combined2$b$barcode=gsub("b_","",combined2$b$barcode)
combined2$c$barcode=gsub("c_","",combined2$c$barcode)
combined2$d$barcode=gsub("d_","",combined2$d$barcode)
combined2$f$barcode=gsub("f_","",combined2$f$barcode)
combined2$g$barcode=gsub("g_","",combined2$g$barcode)
combined2$h$barcode=gsub("h_","",combined2$h$barcode)
head(combined2$h)

#Replace letters (indicating samples) with mouse names#
combined$a$sample <- "WT_1"
combined$b$sample <- "WT_2"
combined$c$sample <- "BCL6_1"
combined$d$sample <- "BCL6_2"
combined$f$sample <- "E1020K_1"
combined$g$sample <- "E1020K_BCL6_1"
combined$h$sample <- "E1020K_BCL6_2"

combined2$a$sample <- "WT_1"
combined2$b$sample <- "WT_2"
combined2$c$sample <- "BCL6_1"
combined2$d$sample <- "BCL6_2"
combined2$f$sample <- "E1020K_1"
combined2$g$sample <- "E1020K_BCL6_1"
combined2$h$sample <- "E1020K_BCL6_2"

#Rename objects in combined list
names(combined) <- c("WT_1", "WT_2", "BCL6_1", "BCL6_2", "E1020K_1", "E1020K_BCL6_1", "E1020K_BCL6_2")
names(combined2) <- c("WT_1", "WT_2", "BCL6_1", "BCL6_2", "E1020K_1", "E1020K_BCL6_1", "E1020K_BCL6_2")


####Merge clonotype data with Seurat object####
experiment_wNA <- combineExpression(combined, 
                                    experiment,
                                    group.by = "sample",
                                    cloneCall= "gene",
                                    proportion = TRUE)
head(experiment_wNA[[]])

experiment_nNA <- combineExpression(combined2, 
                                    experiment,
                                    group.by = "sample",
                                    cloneCall= "gene",
                                    proportion = TRUE)
head(experiment_nNA[[]])

SaveH5Seurat(experiment_wNA, "CITESeq2_AllCells_wNA.h5seurat", overwrite = TRUE) #With all cells, excluding cells with no ADTs. NAs from clonotype data kept
SaveH5Seurat(experiment_nNA, "CITESeq2_AllCells_nNA.h5seurat", overwrite = TRUE) #With all cells, excluding cells with no ADTs. NAs from clonotype data removed

####Load Seurat Object####
experiment <- LoadH5Seurat(file.choose("CITESeq2_AllCells_wNA.h5seurat")) #With clonotype data, including NAs, no BCR genes (except for class inforamtion)
experiment$seurat_clusters <- experiment$wsnn_res.0.1
Idents(experiment) <- experiment$seurat_clusters

#Or

experiment <- LoadH5Seurat(file.choose("CITESeq2_AllCells_nNA.h5seurat"))#With clonotype data, excluding NAs, no BCR genes (except for class inforamtion)
experiment$seurat_clusters <- experiment$wsnn_res.0.1
Idents(experiment) <- experiment$seurat_clusters

#DimPlot all cells
names(col_con1) <- levels(experiment$seurat_clusters) #One colour associated with one cluster
UMAP1 <- DimPlot(experiment, reduction = "wnn.umap", cols = col_con1, pt.size = 1, label = TRUE, label.box = TRUE, label.size = 8) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("All cells") +
  theme(plot.title = element_text(color="black", size= 30, face="bold"),
        legend.text = element_text(size = 20),
        text = element_text(size = 20),
        legend.position = "bottom") +
  NoLegend() #Can remove if legend is needed
print(UMAP1)
ggsave("UMAP1.tiff", width = 30, height = 30, units = "cm", UMAP1, compression = "lzw")


####Subsetting into B, T and Other cells####
DefaultAssay(experiment) <- "ADT"
Unmapped <- FeaturePlot(experiment, features = "Unmapped", reduction = "wnn.umap", pt.size = 1, order = FALSE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD19 cell surface expression") +
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

CD19 <- FeaturePlot(experiment, features = "Cd19", reduction = "wnn.umap", pt.size = 1, order = FALSE) +
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

CD4 <- FeaturePlot(experiment, features = "Cd4", reduction = "wnn.umap", pt.size = 1, order = FALSE) +
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

CD8a <- FeaturePlot(experiment, features = "Cd8a", reduction = "wnn.umap", pt.size = 1, order = FALSE) +
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

DefaultAssay(experiment) <- "RNA"
Ncr1 <- FeaturePlot(experiment, features = "Ncr1", reduction = "wnn.umap", pt.size = 1, order = FALSE) +
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


Prdm1 <- FeaturePlot(experiment, features = "Prdm1", reduction = "wnn.umap", pt.size = 1, order = FALSE) +
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


