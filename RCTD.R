rm(list = ls())
library(spacexr)
library(Seurat)
library(tidyverse)
# dir.create("output/onlyEpi_and_fibRCTD/")
ref <- qs::qread("tmp/OSCC_annotated_updata_STref.qs")###ref <- qs::qread("tmp/OSCC_annotated_STref.qs")
for (gsm_id in 6339632:6339642) {
data_dir <- paste0("data/OSCC_ST/OSCC_GSM", gsm_id, "/")
filename <- "filtered_feature_bc_matrix.h5"
# 读入Spatial数据
seu_vis <- Load10X_Spatial(data.dir = data_dir, filename = filename)
ref$celltype <- sub("/", "_", ref$celltype)
unique(ref$celltype)
Idents(ref) <- ref$celltype
table(ref$celltype) %>% sort()
ref.ds <- subset(ref, downsample = 350)####350
table(ref.ds$celltype) %>% sort()
# ref.ds <- subset(ref.ds, Subset != "Mast")
# extract information to pass to the RCTD Reference function
counts <- ref.ds[["RNA"]]@counts
cluster <- ref.ds$celltype # should be named factors
cluster <- as.factor(cluster)
cluster <- droplevels(cluster)
nUMI <- colSums(counts)
reference <- Reference(counts, cluster, nUMI) # only support <= 10,000 cells
class(reference)
# set up query with the RCTD function SpatialRNA   ###读入st数据
# seu_vis <- Load10X_Spatial(data.dir = "data/OSCC_ST/OSCC_GSM6339631/",
#                            filename = "filtered_feature_bc_matrix.h5")
counts <- seu_vis[["Spatial"]]@counts
coords <- GetTissueCoordinates(seu_vis)
colnames(coords) <- c("x", "y")
query <- SpatialRNA(coords, counts, colSums(counts))
class(query)
# deconvolution
RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# dir.create("tmp")
# qs::qsave(RCTD, "tmp/OSCCall_RCTD.qs")
head(RCTD@results$results_df)
# RCTD@results$weights[1:10,1:10]
dim(RCTD@results$weights)
seu_vis <- AddMetaData(seu_vis, metadata = RCTD@results$results_df)
Idents(seu_vis) <- seu_vis$first_type
table(is.na(Idents(seu_vis)))
# qs::qsave(seu_vis, "tmp/seu_vis.RCTD_OSCC1.qs")
# SpatialDimPlot(seu_vis, cells.highlight = CellsByIdentities(object = seu_vis, idents = levels(Idents(seu_vis))),
#                facet.highlight = T, ncol = 5)
####试试weight
weights.df <- as.data.frame(as.matrix(RCTD@results$weights))

seu_vis <- AddMetaData(seu_vis, metadata = weights.df)
colnames(seu_vis@meta.data)[13:21]

library(ggplot2)  # 确保已加载 ggplot2 包
output_pdf_file_p1 <- paste0("output/onlyEpi_and_fibRCTD/OSCC_GSM", gsm_id, "_P1.pdf")
p1 <- SpatialFeaturePlot(seu_vis, features = colnames(seu_vis@meta.data)[13:21], ncol = 5)
  # 保存为 PDF
ggsave(filename = output_pdf_file_p1, plot = p1,width = 23, height = 7)
###ggsave(filename = output_pdf_file_p1, plot = p1,width = 23, height = 14)


output_pdf_file_p3 <- paste0("output/onlyEpi_and_fibRCTD/OSCC_GSM", gsm_id, "_P3.pdf")
p3 <-  pheatmap::pheatmap(cor(weights.df))
  # 保存为 PDF
ggsave(filename = output_pdf_file_p3, plot = p3,width = 8, height = 5)
}

