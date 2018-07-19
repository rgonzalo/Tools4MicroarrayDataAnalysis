## ----setDirs,echo=FALSE--------------------------------------------------
mainDir <-getwd()
workingDir <- mainDir
dataDir <-file.path(mainDir, "dades")
celDir <-  file.path(workingDir, "celfiles")
celDir2 <-  file.path(workingDir, "celfiles2")
#resultsDir <- file.path(workingDir, "results")
resultsDir <- file.path(workingDir, "results/ABS")
imagesDir<-file.path(mainDir,"images")

## ----loadData2,echo=FALSE------------------------------------------------
load("afterTopTabs.Rda")

<<GOstats1,echo=FALSE,eval=FALSE>>=
library(genefilter)
library(pd.clariom.s.mouse.ht)
library(clariomsmousehttranscriptcluster.db)
library(org.Mm.eg.db)
library(GOstats)
library(annotate)
library(xtable)
library(ReportingTools)
library(KEGG.db)
library(DBI)#necesaría para correr en el cluster (sino error x dbGetQuery)

#creamos una lista con unos cutoffs
selected1<-csv2topTab_FT.AvsCTL[which(csv2topTab_FT.AvsCTL$P.Val<0.01 & 
                                                   abs(csv2topTab_FT.AvsCTL$logFC)>0.5),]
dim(selected1) #464 8
selected2<-csv2topTab_FT.A.C9vsCTL[which(csv2topTab_FT.A.C9vsCTL$P.Val<0.05 & 
                                            abs(csv2topTab_FT.A.C9vsCTL$logFC)>0),]
dim(selected2) #495 8
selected3<-csv2topTab_FT.A.C9vsFT.A[which(csv2topTab_FT.A.C9vsFT.A$P.Val<0.05 & 
                                                   abs(csv2topTab_FT.A.C9vsFT.A$logFC)>0.5),]
dim(selected3) #544 8
selected4<-csv2topTab_DementiavsControls[which(csv2topTab_DementiavsControls$P.Val<0.025 & 
                                            abs(csv2topTab_DementiavsControls$logFC)>0),]
dim(selected4) #495 8

#creem el UNIVERSO a partir de tots els gens que han passat el filtratge
#filtramos para quitar los controles y los que no tengan ENTREZ
filtered2universe <- featureFilter(eset_filtered, require.entrez=TRUE,
    require.GOBP=TRUE, require.GOCC=TRUE,
    require.GOMF=TRUE, require.CytoBand=FALSE,
    remove.dupEntrez=TRUE, feature.exclude="^AFFX")
#dim(filtered2universe) #8042
uni<-exprs(filtered2universe)
EntrezUni <- getEG (rownames(uni), annotation(eset_rma))

#creem les llistes necessaries
listofcomparisons<-list(topTab_FT.AvsCTL,topTab_FT.A.C9vsCTL,topTab_FT.A.C9vsFT.A,topTab_DementiavsControls)
listofcomparisonsnames<-c("topTab_FT.AvsCTL","topTab_FT.A.C9vsCTL","FT.A.C9vsFT.A","topTab_DementiavsControls") 
listsOfResults <- list(GOResults=list(), KEGGResults=list())
listOfCategories<-list("CC","BP","MF")
listOfData<-list(selected1,selected2,selected3,selected4)

#------------------------------------------------ 
#ANALISIS DE LA GO PRIMERO Y DESPUES DE LA KEGG POR SEPARADO.
#------------------------------------------------ 


for (i in 1:length(listofcomparisons)){
  data<-listOfData[[i]]
  geneList<-data$Entrez
  symbolsList <- data$Gene.Symbol
  # Seleccionamos las listas para el an?lisis
  comparison = listofcomparisonsnames[i]
  myLists <- listofcomparisons[[i]]
  entrezUniverse <-  as.character(EntrezUni)
  geneIds <-as.character(geneList)
 
  #Loop para las diferentes categorías de la GO
  for (u in 1:length(listOfCategories)){
  # Creamos los "hiperparametros" en que se basa el an?lisis
   GOparams = new("GOHyperGParams",
     geneIds=geneIds, 
     universeGeneIds=entrezUniverse,
     annotation="org.Hs.eg.db", 
     ontology=listOfCategories[[u]],
     pvalueCutoff=0.05, 
     conditional=FALSE,
     testDirection="over")
  KEGGparams = new("KEGGHyperGParams",
     geneIds=geneIds, 
     universeGeneIds=entrezUniverse,
     annotation="org.Hs.eg.db",  
     pvalueCutoff=0.05, 
     testDirection="over")
   # Ejecutamos el análisis de la GO 
   GOhyper = hyperGTest(GOparams)   
    cat("\nComparison: ", comparison,"\n")
    cat("GO: ",listOfCategories[[u]],"\n")
    print(head(summary(GOhyper)))

      #Creamos un informe html con los resultados
        ####Lo comentamos si se crea el informe con el paquete Reporting tools
        # GOfilename =file.path(resultsDir,
        # paste("GOResults",listOfCategories[[u]],listofcomparisonsnames[[i]],"html", sep="."))
        # htmlReport(GOhyper, file = GOfilename, summary.args=list("htmlLinks"=TRUE))
        #   } 
      ####Lo comentamos si se crea el informe con el paquete Reporting tools 
      #   KEGGfilename =file.path(resultsDir, 
      #     paste("KEGGResults",listofcomparisonsnames[[i]],"html", sep="."))
      #   htmlReport(KEGGhyper, file = KEGGfilename, summary.args=list("htmlLinks"=TRUE))

  #------------------------------------------------  
  #crea un report html para GO amb links amb Reporting Tools
  #---------------------------------------------------
  
  GOfilename =paste("GOResults",listOfCategories[[u]],listofcomparisonsnames[[i]],"html", sep=".")
    
  GOReport <- HTMLReport(shortName = GOfilename,
                        title =paste("Analisys of Biological Significance with GO for ",
                        listofcomparisonsnames[[i]],sep=""),reportDirectory = "./results")
  
  publish(GOhyper, GOReport, selectedIDs = geneIds,
          annotation.db = "org.Hs.eg.db", categorySize = 3)
  
  finish(GOReport)
   }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ara es faria la KEGG
for (i in 1:length(listofcomparisons)){
data<-listOfData[[i]]
  geneList<-data$Entrez
  symbolsList <- data$Gene.Symbol
  # Seleccionamos las listas para el an?lisis
  comparison = listofcomparisonsnames[i]
  myLists <- listofcomparisons[[i]]
  entrezUniverse <-  as.character(EntrezUni)
  geneIds <-as.character(geneList)
KEGGparams = new("KEGGHyperGParams",
     geneIds=geneIds, 
     universeGeneIds=entrezUniverse,
     annotation="org.Hs.eg.db",  
     pvalueCutoff=0.05, 
     testDirection="over")
  KEGGhyper = hyperGTest(KEGGparams)
   cat("KEGG\n")
   print(head(summary(KEGGhyper)))

runMulticore<-"0"
GOTerms2Genes.sql <- function(hgResult, anotPackage)
{
  selectedGOTerms <- intersect(names(geneIdUniverse(hgResult)), summary(hgResult)[, 1])
  selectedGO<- geneIdUniverse(hgResult)[selectedGOTerms]
       if (runMulticore ==1 || runMulticore ==3) { 
           selectedGenes <- mclapply(selectedGO, function(x) {intersect(geneIds(hgResult),x)})
         } else {
           selectedGenes <- lapply(selectedGO, function(x) {intersect(geneIds(hgResult),x)})
       }
  
  sql.ENTREZSYMBOL <- "SELECT gene_id, symbol FROM genes, gene_info WHERE genes._id=gene_info._id"

  if (regexpr(".db", anotPackage) < 0)
  {
     genesAnnot <- dbGetQuery(eval(parse(text = paste(anotPackage, "_dbconn()", sep = ""))), sql.ENTREZSYMBOL)
  }else{
    genesAnnot <- dbGetQuery(eval(parse(text = paste(substr(anotPackage, 1, nchar(anotPackage) - 3), "_dbconn()", sep = ""))), sql.ENTREZSYMBOL)
  }
  
  if (runMulticore ==1 || runMulticore ==3) { 
      selectedSymb <- mclapply(selectedGenes, function(x) {genesAnnot[which(unlist(genesAnnot) %in% x), ]})
      selectedSymb <- mclapply(selectedSymb, function(x) {x <- x[, -1]})
    } else {
      selectedSymb <- lapply(selectedGenes, function(x) {genesAnnot[which(unlist(genesAnnot) %in% x), ]})
      selectedSymb <- lapply(selectedSymb, function(x) {x <- x[, -1]})         
   }

  return(selectedSymb)
}

####CANVIAR EL NOM DE L'ORGANISME
organisme<-"hsa"
#####!!!!!!!!!!!!!!!!!!!!!!!!

addGeneNames<-TRUE

   getSymbol <- function (x)
  {
    if (length(x)>0)
    {
      simbols <- getSYMBOL(x, anotPack)
    }else{
      simbols <- NULL
    }

    return(simbols)
  }
  
  hyperRes <- KEGGhyper
  sumari <- summary(hyperRes)
  fName <- "KEGG Enrichment Analysis" # Informe en HTML
     
  if (addGeneNames)
  {    
    EnrichedKEGGTerms <- as.character(sumari[, 1])
        
    if (length(EnrichedKEGGTerms) > 0)
    {
      ###CANVIAR EL NOM DE LES ANOTACIONS
      anotPackage<-"clariomdhumantranscriptcluster.db"
      #####!!!!!!!!!!!!!!!!!!!!!!!!
      selectedSymbols <- GOTerms2Genes.sql(hyperRes, anotPackage)
      genesInHg <- sapply(selectedSymbols, function(x) paste(x, collapse = ", "))
      Report <- cbind(sumari, GeneNames = genesInHg)
    }else{
      Report <- sumari
    }
  }else{
    Report <- sumari
  }

  Report[, 1] <- paste("<a href=\"http://www.genome.jp/dbget-bin/www_bget?pathway+", organisme, Report[, 1], "\">", Report[, 1], "</a>", sep="") 

  ReportSig <- Report[1:nrow(sumari),]

HTML(ReportSig,file=file.path(resultsDir,
                              paste("KEGGResults",listofcomparisonsnames[[i]],"html",sep = ".")))
}

