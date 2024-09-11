## ##CITE-Seq (2): Seurat Object Generation####
#Julius Christopher Baeck

####Note: The Seurat object(s) have been generated with the Seurat version 4 pipeline and associated dependencies.
####      Seurat version 5 is now implemented with updated dependencies.
####      UMAP clustering is dependent on the versions of individual packages used for the data analysis.
####      Consequently, the UMAP plot(s) will look marginally different for each individual, depending on the versions of each package used for analysis.
####      Cell numbers, gene expression and clustering however will always stay the same, hence there is no difference in the biology, just in the representation on the UMAP.
####      The cell ranger output files can be provided to run the code below.
####      Furthermore, the Seurat object used for the data analysis within the PhD thesis can also be provided if results want to be reproduced.


####Setup####
#Set working directory folder of containing files
#Load required packages
library(Seurat)
library(future)
library(stringr)
library(Matrix)
library(ggplot2)
library(SeuratDisk)
library(clustree)
library(scRepertoire)
library(harmony)
sessionInfo()

#Alter working capacity
options(future.globals.maxSize= 10097152000) # 10Gb

####Load the 10X Cell Ranger output####
#Read the 10x Cell Ranger Output
#c1_ge.data <- Read10X(data.dir = "INPUT FILE PATH HERE") #Mouse KOCB10.3d
#c2_ge.data <- Read10X(data.dir = "INPUT FILE PATH HERE") #Mouse KOCB10.3g
#d1_ge.data <- Read10X(data.dir = "INPUT FILE PATH HERE") #Mouse KOBC10.3f
#d2_ge.data <- Read10X(data.dir = "INPUT FILE PATH HERE") #Mouse KOBC13.2b
#a_ge.data <- Read10X(data.dir = "INPUT FILE PATH HERE") #Mouse KOBC16.3g
#b_ge.data <- Read10X(data.dir = "INPUT FILE PATH HERE") #Mouse KOBC16.3c
#f_ge.data <- Read10X(data.dir = "INPUT FILE PATH HERE") #Mouse KOBC16.3a


#Add the sample to the cell names, consistent with antibody data below
colnames(c1_ge.data)=gsub("-1","_c1",colnames(c1_ge.data))
colnames(c2_ge.data)=gsub("-1","_c2",colnames(c2_ge.data))
colnames(d1_ge.data)=gsub("-1","_d1",colnames(d1_ge.data))
colnames(d2_ge.data)=gsub("-1","_d2",colnames(d2_ge.data))
colnames(a_ge.data)=gsub("-1","_a",colnames(a_ge.data))
colnames(b_ge.data)=gsub("-1","_b",colnames(b_ge.data))
colnames(f_ge.data)=gsub("-1","_f",colnames(f_ge.data))

#Capitalise with lower case letters (genes)
rownames(c1_ge.data)=str_to_title(rownames(c1_ge.data))
rownames(c2_ge.data)=str_to_title(rownames(c2_ge.data))
rownames(d1_ge.data)=str_to_title(rownames(d1_ge.data))
rownames(d2_ge.data)=str_to_title(rownames(d2_ge.data))
rownames(a_ge.data)=str_to_title(rownames(a_ge.data))
rownames(b_ge.data)=str_to_title(rownames(b_ge.data))
rownames(f_ge.data)=str_to_title(rownames(f_ge.data))
head(f_ge.data)

####Load 10X Antibody data####
#Read the 10x Antibody output
#c1_ab.data <- Read10X(data.dir = "INPUT FILE PATH HERE") #Mouse KOCB10.3d
#c2_ab.data <- Read10X(data.dir = "INPUT FILE PATH HERE") #Mouse KOCB10.3g
#d1_ab.data <- Read10X(data.dir = "INPUT FILE PATH HERE") #Mouse KOBC10.3f
#d2_ab.data <- Read10X(data.dir = "INPUT FILE PATH HERE") #Mouse KOBC13.2b
#a_ab.data <- Read10X(data.dir = "INPUT FILE PATH HERE") #Mouse KOBC16.3g
#b_ab.data <- Read10X(data.dir = "INPUT FILE PATH HERE") #Mouse KOBC16.3c
#f_ab.data <- Read10X(data.dir = "INPUT FILE PATH HERE") #Mouse KOBC16.3a

#Tidy up the rownames from the data
rownames(c1_ab.data)=gsub("-[^-]+$","",rownames(c1_ab.data),perl=TRUE)
rownames(c2_ab.data)=gsub("-[^-]+$","",rownames(c2_ab.data),perl=TRUE)
rownames(d1_ab.data)=gsub("-[^-]+$","",rownames(d1_ab.data),perl=TRUE)
rownames(d2_ab.data)=gsub("-[^-]+$","",rownames(d2_ab.data),perl=TRUE)
rownames(a_ab.data)=gsub("-[^-]+$","",rownames(a_ab.data),perl=TRUE)
rownames(b_ab.data)=gsub("-[^-]+$","",rownames(b_ab.data),perl=TRUE)
rownames(f_ab.data)=gsub("-[^-]+$","",rownames(f_ab.data),perl=TRUE)

#Capitalise with lower case letters (antobodies)
rownames(c1_ab.data)=str_to_title(rownames(c1_ab.data))
rownames(c2_ab.data)=str_to_title(rownames(c2_ab.data))
rownames(d1_ab.data)=str_to_title(rownames(d1_ab.data))
rownames(d2_ab.data)=str_to_title(rownames(d2_ab.data))
rownames(a_ab.data)=str_to_title(rownames(a_ab.data))
rownames(b_ab.data)=str_to_title(rownames(b_ab.data))
rownames(f_ab.data)=str_to_title(rownames(f_ab.data))
head(f_ab.data)

#Add the Sample to the cell names in each sample
colnames(c1_ab.data)=paste(colnames(c1_ab.data),"_c1",sep="")
colnames(c2_ab.data)=paste(colnames(c2_ab.data),"_c2",sep="")
colnames(d1_ab.data)=paste(colnames(d1_ab.data),"_d1",sep="")
colnames(d2_ab.data)=paste(colnames(d2_ab.data),"_d2",sep="")
colnames(a_ab.data)=paste(colnames(a_ab.data),"_a",sep="")
colnames(b_ab.data)=paste(colnames(b_ab.data),"_b",sep="")
colnames(f_ab.data)=paste(colnames(f_ab.data),"_f",sep="")
head(f_ab.data)

####Combine 10X Cell Ranger and Antibody Data into a Suerat Object####
m <- Matrix(nrow = nrow(c1_ab.data), ncol = ncol(c1_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(c1_ab.data)
colnames(m)=colnames(c1_ge.data)
common=intersect(colnames(c1_ge.data),colnames(c1_ab.data))
m[,common]=c1_ab.data[,common]
c1 = CreateSeuratObject(counts = c1_ge.data,project="c1", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
c1[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(c2_ab.data), ncol = ncol(c2_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(c2_ab.data)
colnames(m)=colnames(c2_ge.data)
common=intersect(colnames(c2_ge.data),colnames(c2_ab.data))
m[,common]=c2_ab.data[,common]
c2 = CreateSeuratObject(counts = c2_ge.data,project="c2", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
c2[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(d1_ab.data), ncol = ncol(d1_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(d1_ab.data)
colnames(m)=colnames(d1_ge.data)
common=intersect(colnames(d1_ge.data),colnames(d1_ab.data))
m[,common]=d1_ab.data[,common]
d1 = CreateSeuratObject(counts = d1_ge.data,project="d1", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
d1[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(d2_ab.data), ncol = ncol(d2_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(d2_ab.data)
colnames(m)=colnames(d2_ge.data)
common=intersect(colnames(d2_ge.data),colnames(d2_ab.data))
m[,common]=d2_ab.data[,common]
d2 = CreateSeuratObject(counts = d2_ge.data,project="d2", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
d2[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(a_ab.data), ncol = ncol(a_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(a_ab.data)
colnames(m)=colnames(a_ge.data)
common=intersect(colnames(a_ge.data),colnames(a_ab.data))
m[,common]=a_ab.data[,common]
a = CreateSeuratObject(counts = a_ge.data,project="a", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
a[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(b_ab.data), ncol = ncol(b_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(b_ab.data)
colnames(m)=colnames(b_ge.data)
common=intersect(colnames(b_ge.data),colnames(b_ab.data))
m[,common]=b_ab.data[,common]
b = CreateSeuratObject(counts = b_ge.data,project="b", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
b[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(f_ab.data), ncol = ncol(f_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(f_ab.data)
colnames(m)=colnames(f_ge.data)
common=intersect(colnames(f_ge.data),colnames(f_ab.data))
m[,common]=f_ab.data[,common]
f = CreateSeuratObject(counts = f_ge.data,project="f", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
f[["ADT"]] <- adt_assay
head(f[[]])

####Incoperate VDJ data####
#Load contig file
#c1_cl.data <- read.csv("INPUT FILE PATH HERE") #Mouse KOCB10.3d
#c2_cl.data <- read.csv("INPUT FILE PATH HERE") #Mouse KOCB10.3g
#d1_cl.data <- read.csv("INPUT FILE PATH HERE") #Mouse KOBC10.3f
#d2_cl.data <- read.csv("INPUT FILE PATH HERE") #Mouse KOBC13.2b
#a_cl.data <- read.csv("INPUT FILE PATH HERE") #Mouse KOBC16.3g
#b_cl.data <- read.csv("INPUT FILE PATH HERE") #Mouse KOBC16.3c
#f_cl.data <- read.csv("INPUT FILE PATH HERE") #Mouse KOBC16.3a

#match barcode names with GE and ADT data
c1_cl.data$barcode=gsub("-1","_c1",c1_cl.data$barcode)
c2_cl.data$barcode=gsub("-1","_c2",c2_cl.data$barcode)
d1_cl.data$barcode=gsub("-1","_d1",d1_cl.data$barcode)
d2_cl.data$barcode=gsub("-1","_d2",d2_cl.data$barcode)
a_cl.data$barcode=gsub("-1","_a",a_cl.data$barcode)
b_cl.data$barcode=gsub("-1","_b",b_cl.data$barcode)
f_cl.data$barcode=gsub("-1","_f",f_cl.data$barcode)

contig_list <- list(c1_cl.data, c2_cl.data, d1_cl.data, d2_cl.data, a_cl.data, b_cl.data, f_cl.data)
head(contig_list[[1]])

#Generate combined object
combined <- combineBCR(contig_list, samples = c("c1", "c2", "d1", "d2", "a", "b", "f"))
combined[[7]]
str(combined)
head(combined[[7]])

#Make sure barcodes are identical to GE and ADT data
combined$c1$barcode=gsub("c1_","",combined$c1$barcode)
combined$c2$barcode=gsub("c2_","",combined$c2$barcode)
combined$d1$barcode=gsub("d1_","",combined$d1$barcode)
combined$d2$barcode=gsub("d2_","",combined$d2$barcode)
combined$a$barcode=gsub("a_","",combined$a$barcode)
combined$b$barcode=gsub("b_","",combined$b$barcode)
combined$f$barcode=gsub("f_","",combined$f$barcode)
head(combined$f)

####Merge suerat objects####
experiments=c(c1,c2,d1,d2,a,b,f)
experiment_names=c("c1","c2","d1","d2","a","b","f")
experiment<-merge(x= c1, y=c(c2,d1,d2,a,b,f))
head(experiment[[]])

#Merge Seurat object with VDJ data
experiment <- combineExpression(combined, experiment, cloneCall="gene", group.by =  "sample")
head(experiment[[]])

####Quality control, filtering and normalisation####
#Mitochondrial QC metrics
experiment[["percent.Mt"]] <- PercentageFeatureSet(experiment, pattern = "^Mt-")

#Remove where nCount_ADT = 0
DefaultAssay(experiment) <- "ADT"
experiment <- subset(experiment, nCount_ADT > 0)
DefaultAssay(experiment) <- "RNA"

#Visualize QC metrics as violin plot
RNA_QC <- VlnPlot(experiment, features = c("nFeature_RNA", "nCount_RNA", "percent.Mt"))
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

#Generally aim to filter out unique feature counts over 2,500 and less than 200; and percent.mt over 5%
filter_seurat = function(seurat_object){
  
  message("Performing filter by number of genes and mitochondrial percentage.")
  seurat_object = subset(seurat_object, subset = nFeature_RNA > 200  & percent.Mt < 5 & nFeature_RNA < 2500)
  message("Now the object has ", dim(seurat_object)[1], " genes and ", dim(seurat_object)[2], " cells.")
  
  
  return(seurat_object)
}

head(experiment[[]])
experiment <-  filter_seurat(experiment) #22521 genes across 28913 cells

####Saving Seurat Object####
SaveH5Seurat(experiment, filename = "experiment.plain", overwrite = TRUE) #Non-batch corrected Seurat obejct

#Or

experiment <- LoadH5Seurat(file.choose("experiment.plain.h5seurat"))

