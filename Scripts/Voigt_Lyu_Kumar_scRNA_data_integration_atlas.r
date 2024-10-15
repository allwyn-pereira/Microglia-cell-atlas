.libPaths("/CONDAS/users/pereira-a-1/seurat4/lib/R/library")

library(Seurat)
library(stringr)
library(ggplot2)

#####Seurat object creation for Voigt_et_al_HMG_2022 datasets####

Sample_list = list()
Sample_list_2 = list()
Seurat_object_list = list()
Normalised_Seurat_object_list = list()
final_object_list = list()

input_dir = '/SCRATCH-BIRD/users/pereira-a-1/Voigt_et_al_HMG_2022/dataset/'

plot_dir = '/SCRATCH-BIRD/users/pereira-a-1/retina_atlas/plots/'

samples = list.files('/SCRATCH-BIRD/users/pereira-a-1/Voigt_et_al_HMG_2022/dataset/', pattern = 'counts.csv.gz')

sample_names = str_sub(samples, start = 12, end = 24)

sample_names = str_replace_all(sample_names, pattern = '_$', replacement = '')

file_loc = paste0(input_dir, samples)

for (i in 1:length(file_loc)) {
  dataset = read.csv(gzfile(file_loc[i]))
  dataset_count = dataset[,-1:-8]
  rownames(dataset_count) = dataset[,1]
  dataset_count = data.frame(t(dataset_count))
  Sample_list[[i]] <- dataset_count
}

for (i in 1:length(file_loc)) {
  dataset = read.csv(gzfile(file_loc[i]))
  dataset_meta = dataset[,1:8]
  rownames(dataset_meta) = dataset[,1]
  dataset_meta = dataset_meta[,-1]
  Sample_list_2[[i]] <- dataset_meta
}


for (i in 1:21) {
  seurat.data = CreateSeuratObject(counts = Sample_list[[i]], assay = 'RNA', meta.data = Sample_list_2[[i]], project = sample_names[[i]],  min.cells = 3, min.features = 200)
  Seurat_object_list[[i]] <- seurat.data
}

####Seurat object creation for Lyu_et_al_Scientific_Reports_2021 datasets####

Lyu_dir = '/SCRATCH-BIRD/users/pereira-a-1/Lyu_et_al_Scientific_Reports_2021/temp/'

Lyu_sample_name = c(paste0('M', 1:4, '_78y'), paste0('M', 5:8, '_90y'), paste0('P', 1:4, '_78y'), paste0('P', 5:8, '_90y'))

Lyu_ident = c(paste0('Lyu_M', 1:8), paste0('Lyu_P', 1:8))

Lyu_input_dir = paste0(Lyu_dir, Lyu_sample_name)

Lyu_object_list = list()

for (i in 1:length(Lyu_input_dir)) {
  data <- Read10X(data.dir = Lyu_input_dir[i])
  seurat.data <- CreateSeuratObject(counts = data, project = Lyu_ident[i], min.cells = 3, min.features = 200)
  seurat.data[['percentage.mt']] <- PercentageFeatureSet(seurat.data, pattern = '^MT-')
  seurat.data = subset(seurat.data, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percentage.mt < 15)
  Lyu_object_list[[i]] <- seurat.data
}

####Seurat object creation for Kumar_et_al_Nature_Neuroscience_2022 datasets####

Kumar_dir = '/SCRATCH-BIRD/users/pereira-a-1/Kumar_et_al_Nat_Neuro_2022/unprocessed_dataset/'

Kumar_sample_name = paste0('Sample', 1:11)

Kumar_ident = paste0('Kumar_Patient', 1:11)

Kumar_input_dir = paste0(Kumar_dir, Kumar_sample_name)

Kumar_object_list = list()

for (i in 1:length(Kumar_input_dir)) {
  data <- Read10X(data.dir = Kumar_input_dir[i])
  seurat.data <- CreateSeuratObject(counts = data$'Gene Expression', project = Kumar_ident[i], min.cells = 3, min.features = 200)
  seurat.data[['percentage.mt']] <- PercentageFeatureSet(seurat.data, pattern = '^MT-')
  seurat.data = subset(seurat.data, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percentage.mt < 15)
  Kumar_object_list[[i]] <- seurat.data
}

####Integration of Voigt_et_al_2022, Lyu_et_al_2019 and Kumar_et_al_2022 datasets####

comp_object_list = c(Seurat_object_list, Lyu_object_list, Kumar_object_list)

Normalised_Seurat_object_list = list()

final_object_list = list()

for (i in 1:48) {
  object = NormalizeData(comp_object_list[[i]])
  object = FindVariableFeatures(object, selection.method = 'vst', nfeatures = 2000)
  Normalised_Seurat_object_list[[i]] <- object
}

rm(Sample_list)
rm(Sample_list_2)
rm(Seurat_object_list)
rm(Lyu_object_list)
rm(Kumar_object_list)

features = SelectIntegrationFeatures(object.list = Normalised_Seurat_object_list)

for (i in 1:48) {
  object = ScaleData(Normalised_Seurat_object_list[[i]], features = features)
  object = RunPCA(object, features = features)
  final_object_list[[i]] <- object
}

dataset_anchors = FindIntegrationAnchors(object.list = final_object_list, anchor.features = features, reduction = 'rpca')

dataset_integrated = IntegrateData(anchorset = dataset_anchors)

DefaultAssay(dataset_integrated) <- 'integrated'

dataset_integrated = ScaleData(dataset_integrated)
dataset_integrated = RunPCA(dataset_integrated, npcs = 50)

pdf(paste0(plot_dir, 'CNS_integrated_dataset_PCA.pdf'), width = 30, height = 20)
DimHeatmap(dataset_integrated, dims = 1:25, cells = 500, balanced = TRUE)
DimHeatmap(dataset_integrated, dims = 25:50, cells = 500, balanced = TRUE)
ElbowPlot(dataset_integrated, ndims = 50)
dev.off()

ndims = c('30', '35', '40', '45', '50')

for (dim in ndims) {
  dataset_integrated = RunUMAP(dataset_integrated, reduction = 'pca', dims = 1:dim)
  dataset_integrated = FindNeighbors(dataset_integrated, reduction = 'pca', dims = 1:dim)
  dataset_integrated = FindClusters(dataset_integrated, resolution = 0.5)
  pdf(paste0(plot_dir, 'UMAP_CNS_integrated_dataset_dim', dim, '.pdf'), width = 8, height = 8)
  print(DimPlot(dataset_integrated, reduction = 'umap', repel = T, label = T))
  dev.off()
  pdf(paste0(plot_dir, 'UMAP_orig_ident_CNS_integrated_dataset_dim', dim, '.pdf'), width = 16, height = 16)
  print(DimPlot(dataset_integrated, reduction = 'umap', group.by = 'orig.ident', repel = T))
  dev.off()
}


