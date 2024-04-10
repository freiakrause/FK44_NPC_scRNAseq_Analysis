#This is the script to analyse FK44.1 scRNAseq COUNT data provided by BSF
#install.packages("renv")
renv::init()    #creates project library,  lockfile and Rprofile
#renv::snapshot() # updates lockfile abou current-used packages-> can be shared and reproduced by others when restore() is used
#renv::restore() # installs the exact same versions of packages determined in lockfile
#renv::update()) #
#renv::history() #
#renv::install("usethis")
#usethis::create_github_token()
#gitcreds::gitcreds_set()
renv::install("remotes")
renv::install("llvm")
renv::install("ggplot2", "dplyr" ,"RColorBrewer") # package "SingleR" required by tutorial but not found when trying isntallation
renv::install("igraph") # für das scheiß igraph (benötigt für seurat) nach tausend jahren troublehsooting gefunden: ich brauche: sudo apt install build-essential gfortran UND sudo apt install build-essential gfortran, dann gehts
renv::install("Seurat")

#fehlend für tutorial SingleR,celldex, SingleCellExperiment, müssen wahrscheinlich über bioconducter installiert weredn. ahb ich jetzt noch nicht geschafft
library(Seurat)
library(ggplot2)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)

#Loading and Combining output (non normalized) fromCellRanger
NPC_87.data <- Read10X(data.dir = "") #initializes Seurat object with the non normalized data, in the data dir, mtx file, barcodes.tsv files and genes.tsv files have to be
NPC_87 <- CreateSeuratObject(count = NPC_87.data, project = "FK44_NPC_87", min.cells = 3, min.features = 200)
NPC_88.data <- Read10X(data.dir = "") #initializes Seurat object with the non normalized data, in the data dir, mtx file, barcodes.tsv files and genes.tsv files have to be
NPC_88 <- CreateSeuratObject(count = NPC_88.data, project = "FK44_NPC_88", min.cells = 3, min.features = 200)
NPC_91.data <- Read10X(data.dir = "") #initializes Seurat object with the non normalized data, in the data dir, mtx file, barcodes.tsv files and genes.tsv files have to be
NPC_91 <- CreateSeuratObject(count = NPC_91.data, project = "FK44_NPC_91", min.cells = 3, min.features = 200)
NPC_92.data <- Read10X(data.dir = "") #initializes Seurat object with the non normalized data, in the data dir, mtx file, barcodes.tsv files and genes.tsv files have to be
NPC_92 <- CreateSeuratObject(count = NPC_92.data, project = "FK44_NPC:92", min.cells = 3, min.features = 200)
NPC_87.data<- NULL #erase from memory to save RAM
NPC_88.data<- NULL
NPC_91.data<- NULL
NPC_92.data<- NULL
NPC_combined <- merge(NPC_87, y = c(NPC_88, NPC_91, NPC_92), add.cell.ids = c("87","88", "91", "92"), project = "NPC_ALL")
head(colnames(NPC_ALL))
tail(colnames(NPC_ALL))
unique(sapply(X = strsplit(colnames(NPC_ALL), split = "_"), FUN = "[", 1))
table(NPC_ALL$orig.ident)
str(NPC_ALL)
meta <- NPC_ALL@meta.data
dim(meta)
head(meta)
summary(meta$nCount_RNA)
summary(meta$nFeature_RNA)

#Setting things up for qualitycontrol mitochondrial genes, ribosomal content, doublets
NPC_ALL[["percent.mt"]] <- PercentageFeatureSet(NPC_ALL, pattern = "^MT-") #
NPC_ALL[["percent.rb"]] <- PercentageFeatureSet(NPC_ALL, pattern = "^RP[SL]")
doublets <- read.table("data/update/scrublet_calls.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
NPC_ALL <- AddMetaData(NPC_ALL,doublets)
head(NPC_ALL[[]])

#Plot Meta Data
VlnPlot(NPC_ALL, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
FeatureScatter(NPC_ALL, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(NPC_ALL, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(NPC_ALL, feature1 = "nCount_RNA", feature2 = "percent.rb")
FeatureScatter(NPC_ALL, feature1 = "percent.rb", feature2 = "percent.mt")
FeatureScatter(NPC_ALL, feature1 = "nFeature_RNA", feature2 = "Doublet_score")

#Add Conclusions from meta data as QC columns
NPC_ALL[['QC']] <- ifelse(NPC_ALL@meta.data$Is_doublet == 'True','Doublet','Pass')
NPC_ALL[['QC']] <- ifelse(NPC_ALL@meta.data$nFeature_RNA < 500 & NPC_ALL@meta.data$QC == 'Pass','Low_nFeature',NPC_ALL@meta.data$QC)
NPC_ALL[['QC']] <- ifelse(NPC_ALL@meta.data$nFeature_RNA < 500 & NPC_ALL@meta.data$QC != 'Pass' & NPC_ALL@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',NPC_ALL@meta.data$QC,sep = ','),NPC_ALL@meta.data$QC)
NPC_ALL[['QC']] <- ifelse(NPC_ALL@meta.data$percent.mt > 15 & NPC_ALL@meta.data$QC == 'Pass','High_MT',NPC_ALL@meta.data$QC)
NPC_ALL[['QC']] <- ifelse(NPC_ALL@meta.data$nFeature_RNA < 500 & NPC_ALL@meta.data$QC != 'Pass' & NPC_ALL@meta.data$QC != 'High_MT',paste('High_MT',NPC_ALL@meta.data$QC,sep = ','),NPC_ALL@meta.data$QC)
table(NPC_ALL[['QC']])

#Plot only cells that pass QC
VlnPlot(subset(NPC_ALL, subset = QC == 'Pass'), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
#normalize data set to account for sequencing depth, default scale to 10 000 and log2-transform
NPC_ALL_Norm <- NormalizeData(NPC_ALL)
NPC_ALL_Norm <- FindVariableFeatures(NPC_ALL_Norm, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(NPC_ALL_Norm), 10)
plot1 <- VariableFeaturePlot(NPC_ALL_Norm)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

#Scale data
all.genes <- rownames(NPC_ALL_Norm)
NPC_ALL_Scaled <- ScaleData(NPC_ALL_Norm, features = all.genes)
#PCA 
NPC_ALL_Scaled <- RunPCA(NPC_ALL_Scaled, features = VariableFeatures(object = NPC_ALL_Scaled))

VizDimLoadings(NPC_ALL_Scaled, dims = 1:9, reduction = "pca") & 
  theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))
DimHeatmap(NPC_ALL_Scaled, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
DimPlot(NPC_ALL_Scaled, reduction = "pca")
ElbowPlot(NPC_ALL_Scaled) #It’s often good to find how many PCs can be used without much information loss. Look were the big drop happens

#now do clustering, resolution matters, default is 0.8 but play around with it to get best resolution
NPC_ALL_Scaled <- FindNeighbors(NPC_ALL_Scaled, dims = 1:10)
NPC_ALL_Scaled <- FindClusters(NPC_ALL_Scaled, resolution = 0.5)
NPC_ALL_Scaled <- RunUMAP(NPC_ALL_Scaled, dims = 1:10, verbose = F)
table(NPC_ALL_Scaled@meta.data$seurat_clusters)
DimPlot(NPC_ALL_Scaled,label.size = 4,repel = T,label = T)


NPC_ALL_Scaled_QC <- subset(NPC_ALL_Scaled, subset = QC == 'Pass')

DimPlot(NPC_ALL_Scaled_QC,label.size = 4,repel = T,label = T)
FeaturePlot(NPC_ALL_Scaled_QC,features = "percent.mt",label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))

#Also want to get cell cycle to nedd to convert cell cycle gene list which is for human into mous
library(gprofiler2)
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
NPC_ALL_Scaled_QC <- CellCycleScoring(NPC_ALL_Scaled_QC, s.features = mmus_s, g2m.features = mmus_g2m)
table(NPC_ALL_Scaled_QC[[]]$Phase)

VlnPlot(NPC_ALL_Scaled_QC,features = "percent.mt") & theme(plot.title = element_text(size=10))
FeaturePlot(NPC_ALL_Scaled_QC,features = "percent.rb",label.size = 4,repel = T,label = T) & theme(plot.title = element_text(size=10))
VlnPlot(NPC_ALL_Scaled_QC,features = "percent.rb") & theme(plot.title = element_text(size=10))

VlnPlot(NPC_ALL_Scaled_QC,features = c("nCount_RNA","nFeature_RNA")) & 
  theme(plot.title = element_text(size=10))

FeaturePlot(NPC_ALL_Scaled_QC,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
VlnPlot(NPC_ALL_Scaled_QC,features = c("S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))

#Since we have performed extensive QC with doublet and empty cell removal, we can now apply SCTransform normalization, 
#that was shown to be beneficial for finding rare cell populations by improving signal/noise ratio
#Single SCTransform command replaces NormalizeData, ScaleData, and FindVariableFeatures
#We will also correct for % MT genes and cell cycle scores using vars.to.regress variables
#our previous exploration has shown that neither cell cycle score nor MT percentage change very dramatically between clusters,
#so we will not remove biological signal, but only some unwanted variation.
NPC_ALL_TRANSFORM <- SCTransform(NPC_ALL_Scaled_QC, method = "glmGamPoi", ncells = 8824, 
                    vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F)

NPC_ALL_TRANSFORM <- RunPCA(NPC_ALL_TRANSFORM, verbose = F)
NPC_ALL_TRANSFORM <- RunUMAP(NPC_ALL_TRANSFORM, dims = 1:30, verbose = F)
NPC_ALL_TRANSFORM <- FindNeighbors(NPC_ALL_TRANSFORM, dims = 1:30, verbose = F)
NPC_ALL_TRANSFORM <- FindClusters(NPC_ALL_TRANSFORM, verbose = F)
table(NPC_ALL_TRANSFORM[[]]$seurat_clusters)
DimPlot(NPC_ALL_TRANSFORM, label = T)
#Schauen, ob wir Cluster noch correct habeb, durch expression von marker genen von kleinen Populationen platelets und DCs hier (einfüllen,w as ich will) außerdem wurde Farbschmea geändert
FeaturePlot(NPC_ALL_TRANSFORM,"PPBP") & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

FeaturePlot(NPC_ALL_TRANSFORM,"LILRA4") & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))