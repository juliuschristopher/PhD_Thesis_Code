## ##Re-clustering of all B cells in the CITE-Seq (2) dataset####
#Julius Christopher Baeck

####Note: The Seurat object(s) have been generated with the Seurat version 4 pipeline and associated dependencies.
####      Seurat version 5 is now implemented with updated dependencies.
####      UMAP clustering is dependent on the versions of individual packages used for the data analysis.
####      Consequently, the UMAP plot(s) will look marginally different for each individual, depending on the versions of each package used for analysis.
####      Cell numbers, gene expression and clustering however will always stay the same, hence there is no difference in the biology, just in the representation on the UMAP.
####      The cell ranger output files can be provided to run the code below.
####      Furthermore, the Seurat object used for the data analysis within the PhD thesis can also be provided if results want to be reproduced.


####Setup####
#Set working directory to folder containing files
##Packages##
library(Seurat)
library(stringr)
library(Matrix)
library(ggplot2)
library(SeuratDisk)
library(clustree)
library(patchwork)
library(harmony)
sessionInfo()

##Define functions##
###TFIDF
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
} #For downstreal analysis (not essential)

##Define colour palettes##
col_con = viridis(50)
nb.cols <- 13
mycolors <- colorRampPalette(brewer.pal(11, "RdYlBu"))(nb.cols)


####Load Seurat object - experiment.plain####
experiment.plain <- LoadH5Seurat(file.choose("experiment.plain.h5seurat"))
head(experiment.plain[[]])

test <- experiment.plain

##Selection of B cells clonotype + CD19 ADT > 2
DefaultAssay(test) <-  "RNA"
Bcells <- subset(test, subset = cloneType != "NA") #12794
DefaultAssay(Bcells) <-  "ADT"
Bcells <- subset(Bcells, subset = Cd19 > 2) #12636

B_cells <- Bcells
head(B_cells[[]])

###RNA####
###Normalise subset###
DefaultAssay(B_cells) <- "RNA" #For log normalisation

##RNA normalisation
B_cells <- NormalizeData(B_cells, verbose = TRUE)
B_cells <- FindVariableFeatures(B_cells, nfeatures = 3000)
B_cells <- ScaleData(B_cells)

##Visualisation
top20 <-  head(VariableFeatures(B_cells), 20)
plot1.1 <-  VariableFeaturePlot(B_cells)
top20_plot <-  LabelPoints(plot = plot1.1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)

##RNA PCA
B_cells <- RunPCA(B_cells, verbose = FALSE, features = VariableFeatures(object = B_cells))
pca_variance <- B_cells@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #35

##RNA clustering
DefaultAssay(B_cells) <- "RNA"

B_cells <- FindNeighbors(B_cells, dims = 1:35)
B_cells <- FindClusters(B_cells, resolution = 1.0, verbose = FALSE) #1.0 for the resolution
clustree(B_cells, prefix = "RNA_snn_res.") + theme(legend.position="bottom")
B_cells <-RunUMAP(B_cells, dims = 1:35, assay = 'RNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
B_cells_p1 <- DimPlot(B_cells, label = TRUE, reduction = "rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("RNA Clustering") + theme_void() + NoLegend()
B_cells_p1 <- B_cells_p1 + theme(plot.title = element_text(color="black", size=25, face="bold"))

####ADT####
DefaultAssay(B_cells) <- "ADT"

##ADT normalisation
VariableFeatures(B_cells) <- rownames(B_cells[["ADT"]])
B_cells <- NormalizeData(B_cells, normalization.method = "CLR", margin = 2)
B_cells <- ScaleData(B_cells)
B_cells <- RunPCA(B_cells, reduction.name = 'apca', approx = FALSE)

##ADT PCA
apca_variance <- B_cells@reductions$apca@stdev^2
plot(apca_variance/sum(apca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #27

##ADT clustering
B_cells <- FindNeighbors(B_cells, dims = 1:27, reduction = "apca")
B_cells <- FindClusters(B_cells, resolution = 1.0, verbose = FALSE) #1.0 for the resolution
clustree(B_cells, prefix = "ADT_snn_res.") + theme(legend.position="bottom")
B_cells <- RunUMAP(B_cells, reduction = 'apca', dims = 1:27, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
B_cells_p2 <- DimPlot(B_cells, label = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("ADT Clustering") + theme_void() + NoLegend()
B_cells_p2 <- B_cells_p2 + theme(plot.title = element_text(color="black", size=25, face="bold"))

####WNN####
DefaultAssay(B_cells) <- "RNA"

##Combine into wnn plot
B_cells <- FindMultiModalNeighbors(
  B_cells, reduction.list = list("pca", "apca"), 
  dims.list = list(1:35, 1:27), modality.weight.name = "RNA.weight")

##WNN clustering
B_cells <- FindClusters(B_cells, graph.name = "wsnn", algorithm = 3, resolution = 1.0, verbose = TRUE) #1.0 for the resolution
clustree(B_cells, prefix = "wsnn_res.") + theme(legend.position="bottom")
B_cells <- RunUMAP(B_cells, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
B_cells_p3 <- DimPlot(B_cells, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("WNN Clustering") + theme_void() + NoLegend()
B_cells_p3 <- B_cells_p3 + theme(plot.title = element_text(color="black", size=25, face="bold"))

##Change Idents
Idents(B_cells)
B_cells[["old.ident"]] <- Idents(B_cells)
Idents(B_cells) <- B_cells$orig.ident
B_cells<- RenameIdents(B_cells, `c1` = "Mb1 Cyp11a1 KO 1", `c2` = "Mb1 Cyp11a1 KO 2", `d1` = "Mb1 E1020K Cyp11a1 KO 1", `d2` = "Mb1 E1020K Cyp11a1 KO 2", `a` = "WT 1", `b` = "WT 2", `f` = "E1020K")
B_cells[["orig.ident"]] <- Idents(B_cells)
Idents(B_cells) <- B_cells$old.ident
Idents(B_cells)

####Harmony reductions####
Idents(B_cells) <- B_cells$orig.ident
B_cells[["Genotype"]] <- Idents(B_cells)
B_cells <- RenameIdents(B_cells, `Mb1 Cyp11a1 KO 1` = "Batch 1", `Mb1 Cyp11a1 KO 2` = "Batch 1", `Mb1 E1020K Cyp11a1 KO 1` = "Batch 1", `Mb1 E1020K Cyp11a1 KO 2` = "Batch 1", `WT 1` = "Batch 2", `WT 2` = "Batch 2", `E1020K` = "Batch 2")
B_cells[["Batch"]] <- Idents(B_cells)
Idents(B_cells) <- B_cells$old.ident
head(B_cells[[]])

##Batch effect for RNA
DefaultAssay(B_cells) <- "RNA"
B_cells = RunHarmony(B_cells, "Batch", plot_convergence = TRUE, reduction = "pca", assay.use = "RNA", reduction.save = "harmony.rna")
harmony_rna_embeddings = Embeddings(B_cells, 'harmony.rna')
B_cells = B_cells %>%
  RunUMAP(reduction = "harmony.rna", dims = 1:35, assay = 'RNA', reduction.name = 'harmony.rna.umap', reduction.key = 'harmony.rnaUMAP_') %>%
  FindNeighbors(reduction = "harmony.rna", dims = 1:35) %>%
  FindClusters(resolution = 1.0, verbose = FALSE)
B_cells_p4 <- DimPlot(B_cells, label = TRUE, reduction = "harmony.rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("Harmony RNA Clustering") + theme_void() + NoLegend()
B_cells_p4 <- B_cells_p4 + theme(plot.title = element_text(color="black", size=25, face="bold"))


##Batch effect for ADT
DefaultAssay(B_cells) <- "ADT"
B_cells = RunHarmony(B_cells, "Batch", plot_convergence = TRUE, reduction = "apca", assay.use = "ADT", reduction.save = "harmony.adt")
harmony_adt_embeddings = Embeddings(B_cells, 'harmony.adt')
B_cells = B_cells %>%
  RunUMAP(reduction = "harmony.adt", dims = 1:27, assay = 'ADT', reduction.name = 'harmony.adt.umap', reduction.key = 'harmony.adtUMAP_') %>%
  FindNeighbors(reduction = "harmony.adt", dims = 1:27) %>%
  FindClusters(resolution = 1.0, verbose = FALSE)
B_cells_p5 <- DimPlot(B_cells, label = TRUE, reduction = "harmony.adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("Harmony ADT Clustering") + theme_bw() + NoLegend()
B_cells_p5 <- B_cells_p5 + theme(plot.title = element_text(color="black", size=15, face="bold")) + xlab("UMAP1")


##WNN clustering with harmony reductions
DefaultAssay(B_cells) <- "RNA"
B_cells <- FindMultiModalNeighbors(
  B_cells, reduction.list = list("harmony.rna", "harmony.adt"), 
  dims.list = list(1:35, 1:27), modality.weight.name = "harmony.weight", weighted.nn.name = "harmony.weighted.nn", snn.graph.name = "harmony.wsnn")
B_cells <- FindClusters(B_cells, graph.name = "harmony.wsnn", algorithm = 3, resolution = 0.6, verbose = TRUE) #0.6 for the resolution
B_cells <- RunUMAP(B_cells, nn.name = "harmony.weighted.nn", reduction.name = "harmony.wnn.umap", reduction.key = "harmony.wnnUMAP_")
B_cells_p6 <- DimPlot(B_cells, label = TRUE, reduction = "harmony.wnn.umap", pt.size = 1.2, label.size = 3, label.box = TRUE) +  ggtitle("B cell clustering") + theme_bw() + NoLegend()
B_cells_p6 <- B_cells_p6 + theme(plot.title = element_text(color="black", size=15, face="bold")) + xlab("UMAP1") + ylab("UMAP2")

Batch_rna.umap <- DimPlot(B_cells, label = TRUE, reduction = "rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("RNA UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))
Batch_adt.umap <- DimPlot(B_cells, label = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("ADT UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))
Batch_wnn.umap <- DimPlot(B_cells, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("WNN UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))

Batch_harmony.rna.umap <- DimPlot(B_cells, label = TRUE, reduction = "harmony.rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("Harmony RNA UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))
Batch_harmony.adt.umap <- DimPlot(B_cells, label = TRUE, reduction = "harmony.adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("Harmony ADT UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))
Batch_harmony.wnn.umap <- DimPlot(B_cells, label = TRUE, reduction = "harmony.wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("Harmony WNN UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))

Batch_plot <- Batch_rna.umap + Batch_adt.umap + Batch_wnn.umap + Batch_harmony.rna.umap + Batch_harmony.adt.umap + Batch_harmony.wnn.umap

####Save B cells Seurat Object####
SaveH5Seurat(B_cells, filename = "CITESeq1_Bcells", overwrite = TRUE)

#Or
B_cells <- LoadH5Seurat(file.choose("CITESeq1_Bcells.h5seurat"))
