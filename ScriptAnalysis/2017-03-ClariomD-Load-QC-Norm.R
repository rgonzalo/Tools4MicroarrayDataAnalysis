## ----setDirs,echo=FALSE--------------------------------------------------
mainDir <-getwd()
workingDir <- mainDir
dataDir <-file.path(mainDir, "dades")
celDir <-  file.path(workingDir, "celfiles")
resultsDir <- file.path(workingDir, "results")
imagesDir<-file.path(mainDir,"images")

## ----loadpackages, echo=FALSE, results='hide',message=FALSE--------------
library(xtable)
library(Biobase)
library(oligo)
library(arrayQualityMetrics)
library(ggplot2)
library(pd.clariom.d.human)
library(ggrepel)

## ----phenoData1, echo=FALSE, results='asis'------------------------------
targets <-read.table(file.path(dataDir,"targets.PA0481.txt"), 
  header = TRUE, row.names = 1) 

stopifnot(require(xtable))
x.big<-xtable(targets,caption="Targets file showing samples and covariates",label="targettable")
print(x.big,tabular.environment='longtable',floating=FALSE,size="tiny")

## ----readcels, echo=FALSE,results='hide'---------------------------------
celFiles<-list.celfiles(celDir,full.names=TRUE)
rawData<-read.celfiles(celFiles)

## ----preajustes, echo=FALSE,results='hide'-------------------------------
colores <- as.character(targets$Colores)
grupos <- targets$Grupo
batch<-targets$Batch
numSamples <- nrow(targets)
sampleNames <-targets$ShortName
#forma2pca<-c(15,16,17,18,rep(c(15,16),2),rep(c(17,18),2),rep(c(15,16),2),17,18, rep(c(15,16),2))

## ----boxplot2pdf,results='hide',message=FALSE----------------------------
pdf(file.path(resultsDir,"BoxplotRaw.pdf"))
boxplot(rawData, which="all",cex.axis=0.6, col=colores,  las=2, names=sampleNames, main="Boxplot for arrays intensity: Raw Data")
dev.off()

## ----plotPCAfunction, echo=FALSE-------------------------------------------------

plotPCA2 <- function ( datos, labels, factor,title,scale) {
  data <- prcomp(t(datos),scale=scale)
  #ajustos del gràfic
  dataDf <- data.frame(data$x)
  Grupo <- factor
  loads <- round(data$sdev^2/sum(data$sdev^2)*100,1)
  # the graphic
  p1 <- ggplot(dataDf,aes(x=PC1, y=PC2)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(aes(color = Grupo), alpha = 0.55, size = 3) +
    coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5))
  # the graphic with ggrepel
  p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels)) + 
    labs(x = c(paste("PC1",loads[1],"%")),y=c(paste("PC2",loads[2],"%"))) +  
    ggtitle(paste("Principal Component Analysis for:",title,sep=" "))+ 
    theme(plot.title = element_text(hjust = 0.5))
}

## ----plotPCA2Dpdf,results='hide',message=FALSE---------------------------
pdf(file.path(resultsDir,"PCAraw.pdf"))
plotPCA2(exprs(rawData),labels =sampleNames, factor=grupos,title="Grupo",scale = FALSE )
dev.off()

## ----distAnalisis, echo=FALSE,results='hide'-----------------------------
manDist <-  dist(t(exprs(rawData))) 
clust.euclid.average <- hclust(dist(t(exprs(rawData))),method="average")

pdf(file.path(resultsDir,"QCRawData.pdf"))
heatmap (as.matrix(manDist),  col=heat.colors(16))  
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of samples",  hang=-1)
dev.off()

## ----arrayQuality, warning=FALSE,echo=FALSE------------------------------
arrayQualityMetrics(rawData, outdir = file.path(resultsDir, "QCDir.Raw"), 
                   force=TRUE)

## ----normalization.rma,echo=FALSE,results='hide'-------------------------
eset_rma <- rma(rawData)

## ----normBoxPlot,echo=FALSE,results='hide',message=FALSE-----------------
pdf(file.path(resultsDir,"BoxplotNorm.pdf"))
boxplot(eset_rma,main="Boxplot of Normalized data", names=sampleNames, cex.axis=0.6, col=colores,las=2)
dev.off()

## ----plotPCA2DNorm,echo=FALSE,results='hide',message=FALSE---------------
pdf(file.path(resultsDir,"PCAnorm.pdf"))
plotPCA2(exprs(eset_rma),labels =sampleNames, factor=grupos,title="Grupo",scale = FALSE )
dev.off() 

## ----distAnalisis2, echo=FALSE,results='hide'----------------------------
manDist <-  dist(t(exprs(eset_rma))) 
clust.euclid.average <- hclust(dist(t(exprs(eset_rma))),method="average")

pdf(file.path(resultsDir,"QCNormData.pdf"))
heatmap (as.matrix(manDist),  col=heat.colors(16))  
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of samples",  hang=-1)
dev.off()

## ----arrayQuality2,  warning=FALSE,echo=FALSE----------------------------
arrayQualityMetrics(eset_rma, outdir = file.path(resultsDir, "QCDir.Norm"), 
                   force=TRUE)

## ----savedata,echo=FALSE,results='hide'----------------------------------
save(rawData,eset_rma,targets,file="normData.Rda")

