.libPaths("/CONDAS/users/pereira-a-1/seurat4/lib/R/library")

library(Seurat)
library(stringr)
library(ggplot2)

data_dir = '/SCRATCH-BIRD/users/pereira-a-1/retina_atlas/data/'

plot_dir = '/SCRATCH-BIRD/users/pereira-a-1/retina_atlas/plots/'

dataset = readRDS(paste0(data_dir, 'CNS_integrated_dataset_UMAP.RDS'))

macroglia_markers = c('VIM', 'HES1', 'GLUL', 'GFAP', 'APOE', 'RLBP1')

microglia_markers = c('AIF1', 'HLA-DRA', 'C1QB', 'TYROBP', 'CCL2', 'CSF1R', 'P2RY12', 'LILRB2')

endothelial_markers = c('CD34', 'PECAM1', 'LY6E', 'VWF')

pericyte_markers = c('COL8A1', 'CFH', 'TIMP3', 'CETP', 'SPRK2')

bipolar_markers = c('VSX2', 'CAMK2B', 'TRPM1', 'OTX2')

horizontal_markers = c('ONECUT1', 'ONECUT2', 'CALB1', 'CNTNAP2')

amarcrine_markers = c('GAD1', 'CALB1', 'CHAT')

cone_markers = c('ARR3', 'GNGT2', 'GUCA1C')

rod_markers = c('PDE6A', 'RHO', 'CNGA1')

retinal_ganglion_markers = c('NEFL', 'GAP43', 'SNCG')

Bcell_markers = c('MS4A1', 'CD79A', 'CD22')

Tcell_markers = c('GZMB', 'IFNG', 'FOXP3')

dendritic_cell_markers = c('RGS1', 'PDGFB', 'MRC1')

MN_phagocyte_markers = c('CD63', 'ABCA1', 'PILRA', 'VEGFA', 'MMP9', 'CCR2')

cell_type_markers = c('macroglia_markers', 'microglia_markers', 'endothelial_markers', 'pericyte_markers', 'bipolar_markers', 'horizontal_markers', 'amarcrine_markers', 'cone_markers', 'rod_markers', 'retinal_ganglion_markers', 'Bcell_markers', 'Tcell_markers', 'dendritic_cell_markers', 'MN_phagocyte_markers')

DefaultAssay(dataset) <- 'RNA'

for (markers in cell_type_markers) {
    pdf(paste0(plot_dir, 'CNS_integrated_dataset_UMAP_', markers, '.pdf'), width=18, height=18)
    print(FeaturePlot(dataset, features=get(markers), col = c('grey', 'black', 'red'), pt.size = 0.1))
    dev.off()
}