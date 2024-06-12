library(Seurat)
rm(list = ls())
scdata <- Read10X(data.dir = "data/")
scobj <- CreateSeuratObject(counts = scdata,
                            project = "N1",
                            min.cells = 3,
                            min.features = 200)
### metadata 增加分组信息
metadata = scobj@meta.data
scobj@meta.data$group = "N1"
saveRDS(scobj,file = "output/N1_seurat.rds")
library(harmony)
###remotes::install_github('satijalab/seurat-wrappers@community-vignette')
library(SeuratWrappers)
################################################
data1 <- readRDS(file = "")
data2 <- readRDS(file = "")
data3 <- readRDS(file = "")
data4 <- readRDS(file = "")
data5 <- readRDS(file = "")
data6 <- readRDS(file = "")
data7 <- readRDS(file = "")
data8 <- readRDS(file = "")
data9 <- readRDS(file = "")
data10 <- readRDS(file = "")
data11 <- readRDS(file = "")
data12 <- readRDS(file = "")
data13 <- readRDS(file = "")
data14 <- readRDS(file = "")
data15 <- readRDS(file = "")
data16 <- readRDS(file = "")
data17 <- readRDS(file = "")
data18 <- readRDS(file = "")
### 把数据变成list
data <- list(data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14,data15,data16,data17,data18)

################################################
### 2.Merge 合并数据
scobj <- merge(x=data[[1]], y = data[-1])
## 删除所有data前缀的数据,释放空间
rm(list =  ls(pattern="data.*"))
################################################
### 3.数据质控
### 主要PercentageFeatureSet函数计算线粒体含量
### 人类使用pattern = "^MT-"，小鼠使用pattern = "^mt-"
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
### 该操作会在metadata数据里面增加一列叫做percent.mt
metadata <- scobj@meta.data
### 质控数据可视化，使用VlnPlot函数
### nFeature_RNA, number of Feature, 每个细胞中有多少个基因
### nCount_RNA, number of counts, 每个细胞中有多少个counts
### percent.mt, 我们自己增加的列,  每个细胞中线粒体基因的比例
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
### 正式筛选，筛选的是细胞，最终细胞减少
### nFeature_RNA > 200
### nFeature_RNA < 2500
### percent.mt < 5
scobj <- subset(scobj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
#scobj <- subset(scobj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
#################################################################
scobj <- NormalizeData(scobj)
##scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)
scobj <- ScaleData(scobj, features = rownames(scobj))
scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj),reduction.name = "pca")
####remotes::install_version("Matrix", version = "1.6-1.1")
### 不进行批次矫正
scobj <- RunUMAP(scobj,reduction = "pca", dims = 1:20, reduction.name = "umap_naive")
##scobj <- RunUMAP(scobj,reduction = "pca", dims = 1:15, reduction.name = "umap_naive")
### 使用harmony进行批次矫正
####设置随机种子
set.seed(970408)
### 默认reduction = "pca",group.by.vars参数输入的是批次信息,reduction.save 是结果保存的名称
scobj <- RunHarmony(scobj,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
scobj <- RunUMAP(scobj, reduction = "harmony", dims = 1:30,reduction.name = "umap")
####用group试一下
scobj <- RunHarmony(scobj,reduction = "pca",group.by.vars = "group",reduction.save = "harmony1")
scobj <- RunUMAP(scobj, reduction = "harmony1", dims = 1:30,reduction.name = "umap1")

p1 <- DimPlot(scobj, reduction = "umap_naive",group.by = "group")
p2 <- DimPlot(scobj, reduction = "umap1",group.by = "group")
p3 <- DimPlot(scobj, reduction = "umap",group.by = "group")
p1+p2+p3
scobj <- FindNeighbors(scobj, reduction = "harmony", dims = 1:30)
### 设置多个resolution选择合适的resolution
scobj <- FindClusters(scobj, resolution = seq(0.2,1.2,0.1))
### 多个分辨率的分群信息会保存在metadata中
metadata <- scobj@meta.data
library(clustree)
clustree(scobj)
### 选定分辨率0.8作为分群
scobj@meta.data$seurat_clusters <- scobj@meta.data$RNA_snn_res.0.4
Idents(scobj) <- "seurat_clusters"
### 作图展示
DimPlot(scobj, reduction = "umap", label = T)
### 学习Dimplot的两个重要参数 group.by 和split.by
### group.by, 在一张图中展示多个信息
DimPlot(scobj, reduction = "umap", group.by = "group")
### split.by, 分成多个信息来展示
DimPlot(scobj, reduction = "umap", split.by = "group")
###分样本展示
DimPlot(scobj, reduction = "umap", split.by = "orig.ident")
####
library(ggplot2)
P5<- DimPlot(scobj, reduction = "umap", split.by = "orig.ident")
pdf("P5_plot.pdf", width = 150, height = 10)
# 绘制 DimPlot
print(P5)
# 关闭 PDF 文件
dev.off()
### 保存数据
### 保存前瘦身,删除掉scale.data
scobj@assays$RNA@scale.data <- matrix()
### 当前的reduction信息
names(scobj@reductions)
### FeaturePlot 如果不限定reduction
### 先找umap, 再找tsne, 再找pca
FeaturePlot(scobj, features = "MS4A1", order = TRUE)
FeaturePlot(scobj, features = "MS4A1", order = TRUE,reduction = "umap_naive")
FeaturePlot(scobj, features = "MS4A1", order = TRUE,reduction = "umap")
### 删除掉多余的reduction 信息
scobj@reductions$umap_naive <- NULL
scobj@reductions$umap1 <- NULL
FeaturePlot(scobj, features = "MS4A1", order = TRUE)
### 保存注释前的数据
saveRDS(scobj,file = "output/.rds")

library(Seurat)
scobj <- readRDS(file = "output/.rds")
DimPlot(scobj, reduction = "umap", label = T)
marker_genes <- c("MS4A1", "CD79A","CD19")
marker_genes <- c("FCER1A","LILRA4","PTCRA","PLD4","SMPD3","IRF4")###DC
marker_genes <- c("GNLY", "NKG7")
marker_genes <- c("CD3E","CD8A","CD4","IL7R","PTPRC")
marker_genes <- c("CD14","FCGR3A","LYZ","C15orf48","HLA-DQB1","HLA-DRA")
marker_genes <- c("FCER1A",)###DC
marker_genes <- c("HBB","HBA2")
marker_genes <- c("SFN","FXYD3","TACSTD2","S100A14","KRT6A")
marker_genes <- c("MZB1", "DERL3","IGKC" )
marker_genes <- c("TPSB2", "TPSAB1","CPA3","HPGDS","SLC18A2","FCER1A")##肥大
marker_genes <- c("KRT14", "KRT5","KRT17")
marker_genes <- c("LYZ", "CD14","C1QB")
marker_genes <- c("DCN", "COL1A1","COL3A1")
marker_genes <- c("TM4SF1", "PECAM1","VWF","CLDN5")
marker_genes <- c("CXCL8", "G0S2","CSF3R")
marker_genes <- c("ACTA1", "MYL1","MYH2")
marker_genes <- c("AIF1", "FCER1G","LYZ","TYROBP","LST1")###巨噬
marker_genes <- c("DCN", "COL1A2","C1S","COL1A2","COL6A2","C1R")###成纤维
marker_genes <- c("SFN","TACSTD2","S100A14","KRT6A","KRT14", "KRT5","KRT17")###上皮
marker_genes <- c("COX6A2","CKM","MB","DES","ACTA1")###肌细胞
marker_genes <- c("CD3D","CD3E","TRAC","CD2","TRBC2")###T/NK
marker_genes <- c("CD79A","MS4A1","MZB1","IGHG1","DERL3","IGKC")###B细胞
marker_genes <- c("ECSCR.1","PCAT19","VWF","PECAM1","EGFL7")###内皮
marker_genes <- c("EPCAM","SOX2","PROM1","ALCAM")###癌症干细胞
marker_genes <- c("LILRA4","PTCRA","PLD4","SMPD3","IRF4")###DC
marker_genes <- c("IGHG1","IGLG2","IGHG3")
### 使用VlnPlot画marker小提琴图
VlnPlot(scobj, features = marker_genes)
### 使用FeaturePlot画出特征分布图
FeaturePlot(scobj, features = marker_genes, order = TRUE,ncol=2)
head(Idents(scobj))
Idents(scobj) <- "seurat_clusters"
### B.给每个群添加注释
scobj <- RenameIdents(scobj,
                      "0"="T/NK",
                      "1"="Epithelial",
                      "2"="Epithelial",
                      "3"= "Macrophages",
                      "4"= "Fibroblast",
                      "5"= "Epithelial",
                      "6"= "T/NK",
                      "7"= "T/NK",
                      "8"= "Epithelial",
                      "9"= "Epithelial",
                      "10"= "Cancer stem cell",
                      "11"= "Endothelial",
                      "12"= "Epithelial",
                      "13"= "B/Plasma",
                      "14"= "T/NK",
                      "15"= "Fibroblast",
                      "16"= "B/Plasma",
                      "17"= "Epithelial",
                      "18"= "Macrophages",
                      "19"= "Fibroblast",
                      "20"= "Mast",
                      "21"= "Epithelial",
                      "22"= "DC",
                      "23"= "Epithelial",
                      "24"= "Myocytes",
                      "25"= "Epithelial",
                      "26"= "Macrophages",
                      "27"= "Fibroblast",
                      "28"= "Epithelial",
                      "29"= "T/NK",
                      "30"= "Epithelial",
                      "31"= "T/NK",
                      "32"= "Epithelial"
)
head(Idents(scobj))

z### 汇总画图
marker_genes <- c("MS4A1",
                  "GNLY", "NKG7",
                  "CD3E","CD8A","CD4","IL7R",
                  "CD14", "FCGR3A", "LYZ",
                  "FCER1A",
                  "PPBP"
)
FeaturePlot(scobj, features = marker_genes, order = TRUE,ncol=4)
################################################
### 找出所有的marker
all_markers <- FindAllMarkers(object = scobj)
saveRDS(all_markers,file = "output/Seurat_stim_all_markers.rds")
library(dplyr)
top_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))%>%
  slice(1:15) %>%
  ungroup()
################################################
### 确定注释结果三部曲
### A.确认群的个数
head(Idents(scobj))
Idents(scobj) <- "seurat_clusters"
### B.给每个群添加注释
scobj <- RenameIdents(scobj,
                      "0"="Epithelial",
                      "1"="T/NK",
                      "2"="Macrophages",
                      "3"= "Epithelial",
                      "4"= "T/NK",
                      "5"= "Fibroblast",
                      "6"= "Epithelial",
                      "7"= "Epithelial",
                      "8"= "Epithelial",
                      "9"= "Epithelial",
                      "10"= "Epithelial",
                      "11"= "Fibroblast",
                      "12"= "Endothelial",
                      "13"= "B cell",
                      "14"= "B/Plasma",
                      "15"= "Epithelial",
                      "16"= "Epithelial",
                      "17"= "T/NK",
                      "18"= "Epithelial",
                      "19"= "Epithelial",
                      "20"= "Fibroblast",
                      "21"= "T/NK",
                      "22"= "Mast",
                      "23"= "Macrophages",
                      "24"= "T/NK",
                      "25"= "Myocytes",
                      "26"= "DC",
                      "27"= "Epithelial",
                      "28"= "Plasma",
                      "29"= "Fibroblast"
)
head(Idents(scobj))
### C.保存注释的结果
DimPlot(scobj, reduction = "umap", label = T)
metadata <- scobj@meta.data
scobj@meta.data$celltype = Idents(scobj)
saveRDS(scobj,file = "output/annotaion.rds")
##########################################################
### 注释结果可视化
rm(list = ls())
library(Seurat)
scobj <- readRDS(file = "output/annotaion.rds")
################################################
### Dimplot
DimPlot(scobj, reduction = "umap")
DimPlot(scobj, reduction = "umap", label = T)+ ggsci::scale_color_d3("category20")
DimPlot(scobj, reduction = "umap", label = T)+NoLegend()
scCustomize::DimPlot_scCustom(scobj, figure_plot = TRUE,color_seed = 3)

DimPlot(scobj, reduction = "umap",split.by = "group")
DimPlot(scobj, reduction = "umap",split.by = "group",label = T) + NoLegend()
################################################
### 比例改变
## 计算比例
data <- as.data.frame(table(scobj$group,scobj$celltype))
colnames(data) <- c("group","CellType","Freq")
library(dplyr)
df <- data %>%
  group_by(group) %>%
  mutate(Total = sum(Freq)) %>%
  ungroup() %>%
  mutate(Percent = Freq/Total) %>%
  as.data.frame()
df$CellType  <- factor(df$CellType,levels = unique(df$CellType))
write.csv(df,file = "output/cell_percent.csv",row.names = F,quote = F)
library(ggplot2)
p <- ggplot(df, aes(x = group, y = Percent, fill = CellType)) +
  geom_bar(position = "fill", stat="identity", color = 'white', alpha = 1, width = 0.95) +
  #scale_fill_manual(values = mycol) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()
p
### 换种形式
### facet_wrap 可以设置行列参数
s<- ggplot(df,aes(x = group, y = Percent,fill=group))+
  geom_bar(stat="identity")+
  facet_wrap(~ CellType,nrow = 2)+
  theme_bw()
s + scale_fill_manual(values = c("#33CCCC", "#FF9999"))
ggplot(df,aes(x = CellType, y = Percent,fill=group))+
  geom_bar(stat="identity")+
  facet_wrap(~ group,nrow = 2)+
  scale_fill_manual(values = c("#33CCCC", "#FF9999"))+
  theme_bw()#+

################################################
### marker 热图,需要scale data
scobj <- ScaleData(scobj, features = rownames(scobj))
top_markers <- markers.to.plot
DoHeatmap(scobj, features = top_markers)

### Identity的大小修改，通过size参数
DoHeatmap(scobj, features = top_markers,size = 3)

### 基因的大小修改theme(axis.text.y = element_text(size = 8))
DoHeatmap(scobj, features = top_markers,size = 3)+
  theme(axis.text.y = element_text(size = 8))

### subset和downsample 可以随机取每个群的细胞数
DoHeatmap(subset(scobj, downsample = 50), features = top_markers,size = 3)+
  theme(axis.text.y = element_text(size = 8))

library(scRNAtoolVis)
AverageHeatmap(object = scobj,markerGene = top_markers)
## 修改基因的字号
AverageHeatmap(object = scobj,markerGene = top_markers,fontsize = 8)
## 组间比较
AverageHeatmap(object = subset(scobj,group=="CTRL"),markerGene = top_markers,fontsize = 8)+
  AverageHeatmap(object = subset(scobj,group=="STIM"),markerGene = top_markers,fontsize = 8)

expected_name <- c("ITGB4","ITGA6","IGF1R","IGF1")
############################
rm(list = ls())
library(Seurat)
### 读入数据
scobj <- readRDS(file = "output/annotaion.rds")
DimPlot(scobj, label = TRUE)
DimPlot(scobj, label = TRUE,split.by = "group")
########################################################
### 神奇操作, 制作新的分组
scobj$celltype.Normal <- paste(scobj$celltype, scobj$group, sep = "_")
metadata <- scobj@meta.data
Idents(scobj) <- "celltype.Normal"
table(scobj$celltype.Normal)
### 找差异，每个亚群的marker
sce.markers <- FindAllMarkers(object = scobj, 
                              only.pos = TRUE, 
                              min.pct = 0.25,
                              thresh.use = 0.25)
saveRDS(sce.markers,file = "output/OSCC_Seurat_sce_split_markers.rds")
sce.markers <- readRDS(file = "output/Seurat_sce_split_markers.rds")
library(dplyr)
markers <- sce.markers %>%
  filter(p_val_adj < 0.001) 
#devtools::install_github("YuLab-SMU/clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
### 名称转换
gid <- bitr(unique(markers$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
### 交叉合并
colnames(gid)[1] <- "gene"
markers <- merge(markers, gid, by='gene')
### 数据预处理
### https://mp.weixin.qq.com/s/4cV8R4NxNklW4jIzsqLGbg
library(tidyr)
markers <- markers  %>%
  separate(cluster,into = c("celltype","group"),sep = "_",remove = F)%>%
  filter(!celltype %in% c("Eryth","Mk","Mono/Mk Doublets"))
### 多组富集分析
x = compareCluster(ENTREZID ~ celltype+group, data = markers, fun='enrichKEGG')
### 绘图
library(ggplot2)
dotplot(x, label_format=60,x="group") + 
  facet_grid(~celltype)+
  theme(axis.text.x = element_text(angle=45, hjust=1)) 
##################################################
### 基于差异分析的GSEA分析策略(对接GZ00技能)
### 同一群细胞中, 处理和对照的差异
interferon.response <- FindMarkers(scobj, 
                                   ident.1 = "Macrophages_OSCC",
                                   ident.2 = "Macrophages_Normal", 
                                   logfc.threshold = 0)
head(interferon.response, n = 15)
### GSEA 分析
## geneList 三部曲
gene_df <- interferon.response
## 1.获取基因logFC
geneList <- gene_df$avg_log2FC
## 2.命名
names(geneList) =  rownames(gene_df)
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)
head(geneList)
library(clusterProfiler)
#################################################################
### 1.kegg 通路
## 读入kegg gene set
genesets  <- read.gmt("resource/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
y <- GSEA(geneList,TERM2GENE = genesets)
yd <- as.data.frame(y)
### 看整体分布
dotplot(y,showCategory=12,split=".sign")+facet_grid(~.sign)
library(enrichplot)
gseaplot2(y,"KEGG_RIBOSOME",color = "red",pvalue_table = T)
gseaplot2(y,"KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",color = "red",pvalue_table = T)

#################################################################
### 2.hallmarks gene set
genesets <- read.gmt("resource/h.all.v2022.1.Hs.symbols.gmt")
### 主程序GSEA
y <- GSEA(geneList,TERM2GENE = genesets)
yd <- as.data.frame(y)
library(ggplot2)
dotplot(y,showCategory=30,split=".sign")+facet_grid(~.sign)
library(enrichplot)
gseaplot2(y, "HALLMARK_INTERFERON_ALPHA_RESPONSE",color = "red",pvalue_table = T)

#################################################################
