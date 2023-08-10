setwd('D:/work/Pathway_Area/')
library(org.Hs.eg.db)
library(clusterProfiler)
library(Hmisc)

GOKEGG <- function(geneset=NULL, output, filename, trend) {
  if(output == 'GO') {
    ego <- enrichGO(gene          = bitr(geneset,
                                         fromType = 'SYMBOL', toType = 'ENTREZID',
                                         OrgDb = 'org.Hs.eg.db')$ENTREZID,
                    # universe      = names(geneList),
                    OrgDb         = org.Hs.eg.db,
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    GO_ <- as.data.frame(ego)
    write.table(GO_, paste0(filename, '_GO_', trend, '.txt'), 
                sep = '\t', row.names = F, quote = F)
    return(GO_)}
  if(output == 'KEGG') {
    kk <- enrichKEGG(gene          = bitr(geneset,
                                          fromType = 'SYMBOL', toType = 'ENTREZID',
                                          OrgDb = 'org.Hs.eg.db')$ENTREZID,
                     organism      = 'hsa',
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,)
    KEGG_ <- as.data.frame(setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID"))
    write.table(KEGG_, paste0(filename, '_KEGG_', trend, '.txt'), 
                sep = '\t', row.names = F, quote = F)
    return(KEGG_)}
}

DrawEnrich <- function(dat, type, top.number = 5, col="blue", trend){
  if(type == 'BP') {tit <- 'Biological Process of '}
  if(type == 'CC') {tit <- 'Cellular Component of '}
  if(type == 'MF') {tit <- 'Molecular Function of '}
  if(type == 'KEGG') {tit <- 'KEGG Pathway of '}
  if(type == 'GO') {tit <- 'Gene Ontology Enrichment of '}
  dat1 = dat[c(1:top.number),]
  dat1$Description = capitalize(dat1$Description)
  dat1 = dat1[order(dat1$p.adjust),,drop=F]
  dat1$Description = factor(dat1$Description,levels=dat1$Description[length(dat1$Description):1])
  dat1$PValue = -log10(dat1$p.adjust)
  dat1$GeneRatio <- dat1$Count / as.numeric(gsub('^.*/', '', dat1$GeneRatio))
  p = ggplot(dat1,aes(GeneRatio, Description)) +
    geom_point(aes(size=Count,colour=PValue)) +
    scale_colour_gradient(low=col,high="red") + 
    labs(colour=expression(-log[10]("P Value")),size="Gene counts",  
         x="Gene Ratio",y="",title=paste0(tit, trend, '-regulated Genes')) +
    theme_bw() + theme(axis.text.x = element_text(size = 14), axis.text.y=element_text(size=12), 
                       plot.title=element_text(size=20, hjust=1), 
                       legend.text=element_text(size=12), legend.title=element_text(size=14),
                       axis.title.x=element_text(size=18)) +
    scale_x_continuous(limits = c(0,max(dat1$GeneRatio) * 1.2)) 
  return(p)
}

Deg <- read.csv('./clustere3.markers.csv', T, stringsAsFactors=F)
Deg <- read.csv('./clustere63.markers.csv', T, stringsAsFactors=F)
Deg <- read.csv('./clustere83.markers.csv', T, stringsAsFactors=F)
Deg <- read.csv('./clustert3.markers.csv', T, stringsAsFactors=F)
Deg <- read.csv('./clustert63.markers.csv', T, stringsAsFactors=F)
Deg <- read.csv('./sox2.markers.csv', T, stringsAsFactors=F)
Deg <- read.csv('./sox2_2.csv', T, stringsAsFactors=F)
Deg <- read.csv('./clustere.markers.sox2.csv', T, stringsAsFactors=F)
Deg <- read.csv('./olk.markers.csv', T, stringsAsFactors=F)
Deg <- read.csv('./olk.markers2.csv', T, stringsAsFactors=F)
Deg <- read.csv('./tE_ACE2.csv', T, stringsAsFactors=F)
Deg <- read.csv('./0-1.csv', T, stringsAsFactors=F)
Deg <- read.csv('./3-4.csv', T, stringsAsFactors=F)
Deg <- read.csv('./sqh2021.csv', T, stringsAsFactors=F)

setup <- Deg$gene[Deg$Pvalue < 0.05 & Deg$log2_FC > 0]

GO <- GOKEGG(geneset=setup, output='GO', filename='sqh2021', trend='UP')
KEGG <- GOKEGG(geneset=setup, output='KEGG', filename='sqh2021', trend='UP')
DrawEnrich(GO, 'GO', 8, 'blue', 'Up')
DrawEnrich(KEGG, 'KEGG', 5, 'blue', 'Up')

GO <- GOKEGG(geneset=setdown, output='GO', filename='geo_down', trend='DOWN')
KEGG <- GOKEGG(geneset=setdown, output='KEGG', filename='geo_down', trend='DOWN')
DrawEnrich(GO, 'GO', 5, 'blue', 'Down')
DrawEnrich(KEGG, 'KEGG', 5, 'blue', 'Down')
