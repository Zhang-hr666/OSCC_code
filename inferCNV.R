#inferCNV
# setwd("GBM_Recur/")
library(tidyverse)
library(infercnv)
library(Seurat)
scRNA <- qs::qread("tmp/OSCC_annotated_STref.qs")
table(scRNA$Epithelial_0.2resolution)###table(scRNA$celltype)
scRNA <- subset(scRNA,idents = c("Epithelial_01","Epithelial_02",
                                 "Epithelial_03","Epithelial_04",
                                 "Epithelial_05","Epithelial_06","Epithelial_07",
                                 "T/NK","Macrophages","Fibroblast"))
scRNA <- subset(scRNA, downsample = 200)

counts_matrix = Seurat::GetAssayData(scRNA, slot="counts")

cellanno = FetchData(scRNA, vars = "Epithelial_0.2resolution" ) %>% tibble::rownames_to_column(var = "cellbarcode")
write.table(cellanno, "tmp/cnv_cellanno.txt", sep = "\t", col.names = F,row.names =FALSE, quote =F )
head(cellanno)

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,  # 可以直接提供矩阵对象
                                    annotations_file="./tmp/cnv_cellanno.txt",
                                    delim="\t",
                                    gene_order_file="./tmp/hg38_gencode_v27.txt",
                                    ref_group_names=c("T/NK", "Macrophages"))  # 用于设置参考组，正常的细胞类型组别
# qs::qsave(infercnv_obj, file = "tmp/infercnv_obj.qs")

library(infercnv)
options(scipen = 100)
options(error = function() traceback(2))
# 
# #options(bitmapType="Xlib")
# 
# infercnv_obj <- qs::qread("tmp/infercnv_obj.qs")

infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir= "infercnv_output_EPI_ALL",  # dir is auto-created for storing outputs
                              cluster_by_groups=T ,   # cluster
                              hclust_method="ward.D2",
                              plot_steps=F,
                              output_format="pdf")
qs::qsave(infercnv_obj2, file = "tmp/infercnv_obj_10000.qs")
############################成纤维亚群
seu<- qs::qread("tmp/OSCC_annotated_STref.qs")
DefaultAssay(seu) <- "RNA"
Idents(seu) <- seu$combined
table(seu$combined)
# 确保 'combined' 列是字符变量
seu$combined <- as.character(seu$combined)
# 替换 'combined' 列中的值
seu$combined <- recode(seu$combined, "F04 " = "F04", "Resting" = "OSCC_Normal")
# 将 'combined' 列转换回因子变量，并指定期望的水平
seu$combined <- factor(seu$combined, levels = c("B/Plasma","DC","Endothelial","Epithelial01","Epithelial02","Epithelial03","Epithelial04",
                                                "Epithelial05","Epithelial06","Epithelial07","F00","F01","F02","F03",
                                                "F04","F06","F07","F08","F09","F10","Macrophages",
                                                "Mast","Myocytes","OSCC_Normal","T/NK"))
# 检查更新后的水平
table(seu$combined)
scRNA <- subset(seu,idents = c("F00","F01","F02","F03",
                                 "F04","F06","F07","F08","F09","F10",
                                 "OSCC_Normal","T/NK"))
scRNA <- subset(scRNA, downsample = 200)

counts_matrix = Seurat::GetAssayData(scRNA, slot="counts")

cellanno = FetchData(scRNA, vars = "combined" ) %>% tibble::rownames_to_column(var = "cellbarcode")
write.table(cellanno, "tmp/cnv_cellanno.txt", sep = "\t", col.names = F,row.names =FALSE, quote =F )
head(cellanno)

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,  # 可以直接提供矩阵对象
                                    annotations_file="./tmp/cnv_cellanno.txt",
                                    delim="\t",
                                    gene_order_file="./tmp/hg38_gencode_v27.txt",
                                    ref_group_names=c("T/NK"))  # 用于设置参考组，正常的细胞类型组别
# qs::qsave(infercnv_obj, file = "tmp/infercnv_obj.qs")

library(infercnv)
options(scipen = 100)
options(error = function() traceback(2))
# 
# #options(bitmapType="Xlib")
# 
# infercnv_obj <- qs::qread("tmp/infercnv_obj.qs")

infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir= "infercnv_output_fib_ALL",  # dir is auto-created for storing outputs
                              cluster_by_groups=T ,   # cluster
                              hclust_method="ward.D2",
                              plot_steps=F,
                              denoise = TRUE,
                              output_format="pdf")
qs::qsave(infercnv_obj2, file = "tmp/infercnv_obj_10000.qs")



# infer_CNV_obj <- qs::qread("tmp/infercnv_EPI_ALL.qs")
# infer_CNV_obj<-readRDS('./infercnv_output_EPI_ALL/run.final.infercnv_obj')
# expr<-infer_CNV_obj@expr.data
# expr[1:4,1:4]
# data_cnv<-as.data.frame(expr)
# dim(expr)
# colnames(data_cnv)
# rownames(data_cnv)
# sce.all.int <- qs::qread("tmp/OSCC_annotated_STref.qs")
# sce.all.int
# sce.all.int=sce.all.int[,colnames(sce.all.int) %in% rownames(phe)]
# sce.all.int
# identical(colnames(sce.all.int), rownames(phe))
# sce.all.int$celltype = phe$celltype
# table(sce.all.int$celltype )
# meta <- sce.all.int@meta.data