setwd('D:/work/other/single_cell/')
library(Seurat)
library(dplyr)
library(patchwork)
library(SingleR)
library(celldex)
library(pheatmap)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)

add = read.csv('add/mstn related genes.csv', T,
               stringsAsFactors=F)
add.genes = toupper(add$Gene.Symbol)

#### mono ####
dir.path1 = 'mat/MAS0_1.outs/filtered_feature_bc_matrix/'

dat1 = Read10X(dir.path1)
row.names(dat1) = toupper(row.names(dat1))

seu1 = CreateSeuratObject(dat1, project='MAS0',
                          min.cells=3,
                          min.features=200)
seu1

# filter
seu1[["percent.mt"]] = PercentageFeatureSet(seu1, pattern="^MT-")
VlnPlot(seu1, pt.size=0,
        features=c("nCount_RNA", "nFeature_RNA", "percent.mt"))
seu1 = subset(seu1, subset=nCount_RNA > 1000 &
                 nFeature_RNA < 6000 &
                 nFeature_RNA > 200 &
                 percent.mt < 25)

# normalize
seu1 = NormalizeData(seu1, 
                     normalization.method="LogNormalize", 
                     scale.factor=10000)

# variable
seu1 = FindVariableFeatures(seu1, 
                            selection.method="vst", 
                            nfeatures=2000)
top10 = head(VariableFeatures(seu1), 10)
plot1 = VariableFeaturePlot(seu1)
plot2 = LabelPoints(plot=plot1, points=top10, repel=T)
plot1+plot2

# scale
all.genes = rownames(seu1)
seu1 = ScaleData(seu1, features=all.genes)

# PCA
seu1 = RunPCA(seu1, features=VariableFeatures(seu1))
print(seu1[["pca"]], dims=1:5, nfeatures=5)
VizDimLoadings(seu1, dims=1:2, reduction="pca")
DimPlot(seu1, reduction="pca")
DimHeatmap(seu1, dims=1:15, cells=500, balanced=T)

ElbowPlot(seu1)

# cluster
seu1 = FindNeighbors(seu1, dims=1:10)
seu1 = FindClusters(seu1, resolution=0.5)
head(Idents(seu1), 5)

# umap & tsne
seu1 = RunUMAP(seu1, dims=1:10)
DimPlot(seu1, reduction='umap', label=T)

seu1 = RunTSNE(seu1, dims=1:10)
DimPlot(seu1, reduction='tsne', label=T)

# vis
VlnPlot(seu1, features=add.genes)
FeaturePlot(seu1, features=add.genes)

#### ploy ####
dir.path1 = 'mat/MAS0_1.outs/filtered_feature_bc_matrix/'
dir.path2 = 'mat/MAS7dpi_1.outs/filtered_feature_bc_matrix/'
dir.path3 = 'mat/TA0_1.outs//filtered_feature_bc_matrix/'
dir.path4 = 'mat/TA7dpi_1.outs//filtered_feature_bc_matrix/'

dat1 = Read10X(dir.path1)
row.names(dat1) = toupper(row.names(dat1))
seu1 = CreateSeuratObject(dat1, project='MAS0',
                          min.cells=3,
                          min.features=200)
dat2 = Read10X(dir.path2)
row.names(dat2) = toupper(row.names(dat2))
seu2 = CreateSeuratObject(dat2, project='MAS1',
                          min.cells=3,
                          min.features=200)
dat3 = Read10X(dir.path3)
row.names(dat3) = toupper(row.names(dat3))
seu3 = CreateSeuratObject(dat3, project='TA0',
                          min.cells=3,
                          min.features=200)
dat4 = Read10X(dir.path4)
row.names(dat4) = toupper(row.names(dat4))
seu4 = CreateSeuratObject(dat4, project='TA1',
                          min.cells=3,
                          min.features=200)
rm(dat1, dat2, dat3, dat4)
seu1
seu2
seu3
seu4

# QC
seu1[["percent.mt"]] = PercentageFeatureSet(seu1, pattern="^MT-")
VlnPlot(seu1, pt.size=0,
        features=c("nCount_RNA", "nFeature_RNA", "percent.mt"))
seu1 = subset(seu1, subset=nCount_RNA > 1000 &
                nFeature_RNA < 6000 &
                nFeature_RNA > 200 &
                percent.mt < 25)

seu2[["percent.mt"]] = PercentageFeatureSet(seu2, pattern="^MT-")
VlnPlot(seu2, pt.size=0,
        features=c("nCount_RNA", "nFeature_RNA", "percent.mt"))
seu2 = subset(seu2, subset=nCount_RNA > 1000 &
                nFeature_RNA < 8000 &
                nFeature_RNA > 200 &
                percent.mt < 20)

seu3[["percent.mt"]] = PercentageFeatureSet(seu3, pattern="^MT-")
VlnPlot(seu3, pt.size=0,
        features=c("nCount_RNA", "nFeature_RNA", "percent.mt"))
seu3 = subset(seu3, subset=nCount_RNA > 1000 &
                nFeature_RNA < 7000 &
                nFeature_RNA > 200 &
                percent.mt < 25)

seu4[["percent.mt"]] = PercentageFeatureSet(seu4, pattern="^MT-")
VlnPlot(seu4, pt.size=0,
        features=c("nCount_RNA", "nFeature_RNA", "percent.mt"))
seu4 = subset(seu4, subset=nCount_RNA > 1000 &
                nFeature_RNA < 8000 &
                nFeature_RNA > 200 &
                percent.mt < 20)

# list
seu.list = list(seu1, seu2, seu3, seu4)


seu.list = lapply(X=seu.list,
                  FUN=function(x) {
                    # normalize
                    x=NormalizeData(x,
                                    normalization.method="LogNormalize",
                                    scale.factor=10000)
                    # variable
                    x=FindVariableFeatures(x,
                                           selection.method="vst",
                                           nfeatures=2000)
                  })
features = SelectIntegrationFeatures(object.list=seu.list)
top10 = features[1:10]
features = union(features, add.genes)

# batch effect processing
anchors = FindIntegrationAnchors(object.list=seu.list,
                                 anchor.features=features)
combined = IntegrateData(anchorset=anchors)
saveRDS(seu.list, file='seulist.rds')
saveRDS(anchors, file='anchors.rds')
saveRDS(combined, file='combined.rds')

## seu.list = readRDS('seulist.rds')
## combined = readRDS('combined.rds')

rm(seu1, seu2, seu3, seu4, seu.list, anchors)

# continue
DefaultAssay(combined) = "integrated"
## DefaultAssay(combined) = "RNA"
combined = ScaleData(combined, features=features)

# PCA
combined = RunPCA(combined, features=VariableFeatures(combined))
print(combined[["pca"]], dims=1:5, nfeatures=5)
VizDimLoadings(combined, dims=1:2, reduction="pca")
DimPlot(combined, reduction="pca")
DimHeatmap(combined, dims=1:15, cells=500, balanced=T)

p0 = ElbowPlot(combined)
p0 = p0 + geom_line(data=p0$data, mapping=aes(x=dims, y=stdev))
tiff('figures/ElbowPlot.tiff', 
     res=300, width=1200, height=1200)
print(p0)
dev.off()

# cluster
combined = FindNeighbors(combined, dims=1:12)
combined = FindClusters(combined, resolution=0.1)
head(Idents(combined), 5)
combined$seurat_clusters = as.factor(as.numeric(as.character(combined$seurat_clusters)) + 1)

# umap & tsne
combined = RunUMAP(combined, dims=1:12)
DimPlot(combined, reduction='umap', label=T)
DimPlot(combined, reduction='umap', label=F, group.by='orig.ident')

combined = RunTSNE(combined, dims=1:12)
DimPlot(combined, reduction='tsne', label=T)
DimPlot(combined, reduction='tsne', label=F, group.by='orig.ident')

# vis
VlnPlot(combined, features=add.genes)
FeaturePlot(combined, features=add.genes)

#### markers ####


#### identifying cells ####

## mouse = MouseRNAseqData()
## saveRDS(mouse, 'mouse.rds')
## immune = ImmGenData()
## saveRDS(immune, 'immune.rds')
mouse = readRDS('mouse.rds')
immune = readRDS('immune.rds')

mat = GetAssayData(combined, slot='data')
pred = SingleR(test=mat, ref=mouse,
               assay.type.test=1, labels=mouse$label.main)
## pred = SingleR(test=mat, ref=immune,
##                assay.type.test=1, labels=immune$label.main)
cell.table = as.data.frame.matrix(table(pred$labels, 
                                        combined$seurat_clusters))
## cell.table[1,] = cell.table[1,] + cell.table[2,]
## cell.table[15,] = cell.table[15,] + cell.table[16,]
## cell.table[19,] = cell.table[19,] + cell.table[20,]
cell.table = cell.table[!(rownames(cell.table)%in%
                            c('B cells, pro',
                              'NKT',
                              'Tgd')),]
write.csv(cell.table, 'tables/CellTable.csv', 
          quote=F, row.names=T)

rownames(cell.table)[4] = 'Myocytes'
rownames(cell.table)[11] = 'Unidentified'
cell.table = cell.table[order(rownames(cell.table)),]

plotScoreHeatmap(pred)
plotDeltaDistribution(pred, ncol = 3)
pheatmap(log10(cell.table + 10), angle_col=0)

pheatmap(log10(cell.table + 10), angle_col=0,
         color=colorRampPalette(c('white', 'red'))(50),
         cellwidth=18, cellheight=18, 
         border_color='black',
         scale='none', 
         filename='figures/HeatMapCellAnnotation.tiff')
dev.off()

# update plot
cell.labels = pred$labels
cell.labels[cell.labels == 'B cells, pro'] = 'B cells'
cell.labels[cell.labels == 'NKT'] = 'NK cells'
cell.labels[cell.labels == 'Tgd'] = 'T cells'

cell.labels[cell.labels == 'Cardiomyocytes'] = 'Myocytes'
cell.labels[cell.labels == 'Hepatocytes'] = 'Unidentified'
combined$labels = cell.labels

cell_type_cols = c(brewer.pal(9, "Set1"), 
                   "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
                   "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
                   "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00")
tiff('figures/DimPlot.tiff',
     res=300, width=4000, height=2400)
print(DimPlot(combined, group.by=c("seurat_clusters", "labels"),
              reduction="tsne", label=T, repel=T,
              cols=cell_type_cols))
dev.off()
tiff('figures/DimPlot_UMAP.tiff',
     res=300, width=4000, height=2400)
print(DimPlot(combined, group.by=c("seurat_clusters", "labels"),
              reduction="umap", label=T, repel=T,
              cols=cell_type_cols))
dev.off()
## DimPlot(combined, group.by=c("seurat_clusters", "labels"),
##         reduction="tsne")
## DotPlot(combined, features=add.genes)


# diff exp 1
cols = c('#303030','blue','#303030','blue')
shapes = c(16, 16, 17, 17)
comps = list(c('MAS1', 'MAS0'),
             c('TA1', 'TA0'),
             c('TA0', 'MAS0'),
             c('TA1', 'MAS1'))
p = VlnPlot(combined, features=add.genes[1], group.by='orig.ident',
            same.y.lims=F, pt.size=0, slot='data',
            fill.by='feature', y.max=1,
            cols=cols) + 
  stat_compare_means()
p + geom_jitter(mapping=aes(color=ident, shape=ident), 
                data=p$data,
                width=0.2) +
  scale_color_manual(values=cols) +
  scale_shape_manual(values=shapes)
  

# diff exp 2
mat9 = mat[add.genes, ]
ggdata = data.frame(Group=combined$orig.ident)
for(i in 1:9) {
  temp = ggdata
  temp$Gene = row.names(mat9)[i]
  temp$Exp = as.numeric(mat9[i,])
  if(i == 1) {
    dat = temp
  } else {
    dat = rbind(dat, temp)
  }
}
ggdata = dat

cols = c('#d121f7','#1c3bd5','#d121f7','#1c3bd5')
shapes = c(16, 16, 17, 17)
comps = list(c('MAS1', 'MAS0'),
             c('TA1', 'TA0'),
             c('TA0', 'MAS0'),
             c('TA1', 'MAS1'))

p1 = ggplot(data=subset(ggdata, ggdata$Gene==add.genes[1]), 
            aes(x=Group, y=Exp, fill=Group, shape=Group)) + 
  geom_jitter(aes(color=Group, fill=Group), width=0.2) +
  geom_violin(alpha=0.5) + 
  xlab('Group') + ylab('Expression Level') +
  ggtitle(add.genes[1]) +
  stat_compare_means(comparisons=comps) + 
  stat_compare_means(label.y=6.5) +
  theme_pubr() +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=15),
        legend.position='none',
        axis.title=element_text(size=15),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values=shapes)
tiff(paste0('figures/Violin_All_', add.genes[1], '.tiff'),
     res=300, width=1800, height=1800)
print(p1)
dev.off()

#### patches ####
p1 = ggplot(data=subset(ggdata, ggdata$Gene==add.genes[1]), 
            aes(x=Group, y=Exp, fill=Group, shape=Group)) + 
  geom_jitter(aes(color=Group, fill=Group), width=0.2) +
  geom_violin(alpha=0.5) + 
  xlab('') + ylab('Expression Level') +
  ggtitle(add.genes[1]) +
  stat_compare_means(comparisons=comps, size=6) + 
  stat_compare_means(label.y=6.5, size=6) +
  theme_pubr() +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=25),
        legend.position='none',
        text=element_text(face='bold'),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x=element_text(angle=45, hjust=1, size=20)) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values=shapes)
p1

p2 = ggplot(data=subset(ggdata, ggdata$Gene==add.genes[2]), 
            aes(x=Group, y=Exp, fill=Group, shape=Group)) + 
  geom_jitter(aes(color=Group, fill=Group), width=0.2) +
  geom_violin(alpha=0.5) + 
  xlab('') + ylab('Expression Level') +
  ggtitle(add.genes[2]) +
  stat_compare_means(comparisons=comps, size=6) + 
  stat_compare_means(label.y=6.5, size=6) +
  theme_pubr() +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=25),
        legend.position='none',
        text=element_text(face='bold'),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x=element_text(angle=45, hjust=1, size=20)) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values=shapes)
p2

p3 = ggplot(data=subset(ggdata, ggdata$Gene==add.genes[3]), 
            aes(x=Group, y=Exp, fill=Group, shape=Group)) + 
  geom_jitter(aes(color=Group, fill=Group), width=0.2) +
  geom_violin(alpha=0.5) + 
  xlab('') + ylab('Expression Level') +
  ggtitle(add.genes[3]) +
  stat_compare_means(comparisons=comps, size=6) + 
  stat_compare_means(label.y=5.5, size=6) +
  theme_pubr() +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=25),
        legend.position='none',
        text=element_text(face='bold'),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x=element_text(angle=45, hjust=1, size=20)) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values=shapes)
p3

p4 = ggplot(data=subset(ggdata, ggdata$Gene==add.genes[4]), 
            aes(x=Group, y=Exp, fill=Group, shape=Group)) + 
  geom_jitter(aes(color=Group, fill=Group), width=0.2) +
  geom_violin(alpha=0.5) + 
  xlab('') + ylab('Expression Level') +
  ggtitle(add.genes[4]) +
  stat_compare_means(comparisons=comps, size=6) + 
  stat_compare_means(label.y=6.5, size=6) +
  theme_pubr() +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=25),
        legend.position='none',
        text=element_text(face='bold'),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x=element_text(angle=45, hjust=1, size=20)) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values=shapes)
p4

p5 = ggplot(data=subset(ggdata, ggdata$Gene==add.genes[5]), 
            aes(x=Group, y=Exp, fill=Group, shape=Group)) + 
  geom_jitter(aes(color=Group, fill=Group), width=0.2) +
  geom_violin(alpha=0.5) + 
  xlab('') + ylab('Expression Level') +
  ggtitle(add.genes[5]) +
  stat_compare_means(comparisons=comps, size=6) + 
  stat_compare_means(label.y=6.5, size=6) +
  theme_pubr() +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=25),
        legend.position='none',
        text=element_text(face='bold'),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x=element_text(angle=45, hjust=1, size=20)) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values=shapes)
p5


p6 = ggplot(data=subset(ggdata, ggdata$Gene==add.genes[6]), 
            aes(x=Group, y=Exp, fill=Group, shape=Group)) + 
  geom_jitter(aes(color=Group, fill=Group), width=0.2) +
  geom_violin(alpha=0.5) + 
  xlab('') + ylab('Expression Level') +
  ggtitle(add.genes[6]) +
  stat_compare_means(comparisons=comps, size=6) + 
  stat_compare_means(label.y=5, size=6) +
  theme_pubr() +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=25),
        legend.position='none',
        text=element_text(face='bold'),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x=element_text(angle=45, hjust=1, size=20)) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values=shapes)
p6

p7 = ggplot(data=subset(ggdata, ggdata$Gene==add.genes[7]), 
            aes(x=Group, y=Exp, fill=Group, shape=Group)) + 
  geom_jitter(aes(color=Group, fill=Group), width=0.2) +
  geom_violin(alpha=0.5) + 
  xlab('') + ylab('Expression Level') +
  ggtitle(add.genes[7]) +
  stat_compare_means(comparisons=comps, size=6) + 
  stat_compare_means(label.y=4.5, size=6) +
  theme_pubr() +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=25),
        legend.position='none',
        text=element_text(face='bold'),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x=element_text(angle=45, hjust=1, size=20)) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values=shapes)
p7

p8 = ggplot(data=subset(ggdata, ggdata$Gene==add.genes[8]), 
            aes(x=Group, y=Exp, fill=Group, shape=Group)) + 
  geom_jitter(aes(color=Group, fill=Group), width=0.2) +
  geom_violin(alpha=0.5) + 
  xlab('') + ylab('Expression Level') +
  ggtitle(add.genes[8]) +
  stat_compare_means(comparisons=comps, size=6) + 
  stat_compare_means(label.y=6, size=6) +
  theme_pubr() +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=25),
        legend.position='none',
        text=element_text(face='bold'),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x=element_text(angle=45, hjust=1, size=20)) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values=shapes)
p8

p9 = ggplot(data=subset(ggdata, ggdata$Gene==add.genes[9]), 
            aes(x=Group, y=Exp, fill=Group, shape=Group)) + 
  geom_jitter(aes(color=Group, fill=Group), width=0.2) +
  geom_violin(alpha=0.5) + 
  xlab('') + ylab('Expression Level') +
  ggtitle(add.genes[9]) +
  stat_compare_means(comparisons=comps, size=6) + 
  stat_compare_means(label.y=5, size=6) +
  theme_pubr() +
  theme(plot.title=element_text(hjust=0.5, face='bold', size=25),
        legend.position='none',
        text=element_text(face='bold'),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x=element_text(angle=45, hjust=1, size=20)) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values=shapes)
p9


p_violin = p1+p2+p3+p4+p5+p6+p7+p8+p9+
  plot_layout(ncol=3, nrow=3)
tiff(paste0('figures/Violin_All_', 'Genes', '.tiff'),
     res=300, width=5400, height=5400)
print(p_violin)
dev.off()

# diff exp 3
for(i in 1:9) {
  tiff(paste0('figures/FeaturePlot_tsne_', 
              as.character(i),
              '_', add.genes[i],
              '.tiff'),
       res=300, width=1800, height=1800)
  print(FeaturePlot(combined, features=add.genes[i],
                    split.by='orig.ident',reduction='tsne') +
          plot_layout(ncol=2, nrow=2))
  dev.off()
}

# diff exp 4
tiff('figures/DotPlot.tiff', 
     res=300, width=1800, height=2400)
print(DotPlot(combined, features=add.genes, group.by='orig.ident',
              dot.scale=14) +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        coord_flip())
dev.off()

# tab
exp.table = data.frame(gene=add.genes)
# mean
agg.mean = aggregate(x=dat$Exp,
                     by=list(dat$Group, dat$Gene),
                     FUN=mean)
agg.mean1 = agg.mean[agg.mean$Group.1 %in% 'MAS0',]
agg.mean2 = agg.mean[agg.mean$Group.1 %in% 'MAS1',]
agg.mean3 = agg.mean[agg.mean$Group.1 %in% 'TA0',]
agg.mean4 = agg.mean[agg.mean$Group.1 %in% 'TA1',]
exp.table$mean.MAS0 = agg.mean1$x[match(exp.table$gene, agg.mean1$Group.2)]
exp.table$mean.MAS1 = agg.mean2$x[match(exp.table$gene, agg.mean2$Group.2)]
exp.table$mean.TA0 = agg.mean3$x[match(exp.table$gene, agg.mean3$Group.2)]
exp.table$mean.TA1 = agg.mean4$x[match(exp.table$gene, agg.mean4$Group.2)]
rm(agg.mean, agg.mean1, agg.mean2, agg.mean3, agg.mean4)

# median
agg.median = aggregate(x=dat$Exp,
                     by=list(dat$Group, dat$Gene),
                     FUN=median)
agg.median1 = agg.median[agg.median$Group.1 %in% 'MAS0',]
agg.median2 = agg.median[agg.median$Group.1 %in% 'MAS1',]
agg.median3 = agg.median[agg.median$Group.1 %in% 'TA0',]
agg.median4 = agg.median[agg.median$Group.1 %in% 'TA1',]
exp.table$median.MAS0 = agg.median1$x[match(exp.table$gene, agg.median1$Group.2)]
exp.table$median.MAS1 = agg.median2$x[match(exp.table$gene, agg.median2$Group.2)]
exp.table$median.TA0 = agg.median3$x[match(exp.table$gene, agg.median3$Group.2)]
exp.table$median.TA1 = agg.median4$x[match(exp.table$gene, agg.median4$Group.2)]
rm(agg.median, agg.median1, agg.median2, agg.median3, agg.median4)

exp.table$MAS1vsMAS0 = ifelse(exp.table$mean.MAS1 >
                                exp.table$mean.MAS0,
                              'high', 'low')
exp.table$TA1vsTA0 = ifelse(exp.table$mean.TA1 >
                              exp.table$mean.TA0,
                            'high', 'low')
exp.table$TA0vsMAS0 = ifelse(exp.table$mean.TA0 >
                               exp.table$mean.MAS0,
                             'high', 'low')
exp.table$TA1vsMAS1 = ifelse(exp.table$mean.TA1 >
                               exp.table$mean.MAS1,
                             'high', 'low')
exp.table$MAS1vsMAS0.sig = T
exp.table$TA1vsTA0.sig = T
exp.table$TA0vsMAS0.sig = c(T,F,T,F,T,T,F,T,T)
exp.table$TA1vsMAS1.sig = c(F,T,F,F,T,T,T,F,F)
write.csv(exp.table, 'tables/ExpTable.csv',
          quote=F, row.names=F)
write.csv(dat, 'tables/dat.csv', 
          quote=F, row.names=F)

table(combined$labels)
temp = combined
Idents(temp) = 'orig.ident'
avg.cells = as.data.frame(log1p(AverageExpression(temp, 
                                                  verbose = F)$RNA))
avg.cells$gene = rownames(avg.cells)
avg.cells$label = ifelse(avg.cells$gene %in% add.genes,
                         T, F)
rm(temp)

#### Dot ####

p1 = ggplot(avg.cells, aes(x=MAS0, y=MAS1, color=label)) +
  geom_point() + theme_pubr() +
  ggtitle('Average Expression') +
  scale_color_manual(values=c('#303030','#1c3bd5')) +
  theme(legend.position='none',
        plot.title=element_text(face='bold', 
                                hjust=0.5, 
                                size=20))
p1 = LabelPoints(plot=p1, points=add.genes, repel=T)
p1

p2 = ggplot(avg.cells, aes(x=TA0, y=TA1, color=label)) +
  geom_point() + theme_pubr() +
  ggtitle('Average Expression') +
  scale_color_manual(values=c('#303030','#1c3bd5')) +
  theme(legend.position='none',
        plot.title=element_text(face='bold', 
                                hjust=0.5, 
                                size=20))
p2 = LabelPoints(plot=p2, points=add.genes, repel=T)
p2

p3 = ggplot(avg.cells, aes(x=MAS0, y=TA0, color=label)) +
  geom_point() + theme_pubr() +
  ggtitle('Average Expression') +
  scale_color_manual(values=c('#303030','#1c3bd5')) +
  theme(legend.position='none',
        plot.title=element_text(face='bold', 
                                hjust=0.5, 
                                size=20))
p3 = LabelPoints(plot=p3, points=add.genes, repel=T)
p3

p4 = ggplot(avg.cells, aes(x=MAS1, y=TA1, color=label)) +
  geom_point() + theme_pubr() +
  ggtitle('Average Expression') +
  scale_color_manual(values=c('#303030','#1c3bd5')) +
  theme(legend.position='none',
        plot.title=element_text(face='bold', 
                                hjust=0.5, 
                                size=20))
p4 = LabelPoints(plot=p4, points=add.genes, repel=T)
p4

tiff('figures/Scatter_All_1_MAS0_MAS1.tiff',
     res=300, width=1800, height=1800)
print(p1)
dev.off()

tiff('figures/Scatter_All_2_TA0_TA1.tiff',
     res=300, width=1800, height=1800)
print(p2)
dev.off()

tiff('figures/Scatter_All_3_MAS0_TA0.tiff',
     res=300, width=1800, height=1800)
print(p3)
dev.off()

tiff('figures/Scatter_All_4_MAS1_TA1.tiff',
     res=300, width=1800, height=1800)
print(p4)
dev.off()
