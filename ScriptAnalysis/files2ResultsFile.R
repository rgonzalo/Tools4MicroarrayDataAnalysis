#library(devtools)
source("https://raw.githubusercontent.com/alexsanchezpla/scripts/master/Rpackages/Create_Results_Browser/writeResults2Browser.R")

mainDir <-getwd()
workingDir <- mainDir
dataDir <-file.path(mainDir, "dades")
resultsDir <- file.path(workingDir, "results")

#(1) listem els arxius i directoris de la carpeta results
#RECORDAR ficar a la carpeta results la proposta,main report i el targets
FileName <- list.files(resultsDir,pattern=NULL,all.files=FALSE,recursive = FALSE,
                  full.names=FALSE,include.dirs = TRUE,no..=FALSE)

#(2)eliminem aquells arxius que no ens interessen (en principi totes les carpetes i arxius tipus .log o .css, etc)
FileName
File<-FileName[-c(9, 23:113, 157)]
File

#(3) creem la resta de variables que ha de tenir l'arxiu
Category<-rep("a",length(File))
Subcateg<-rep("NA",length(File))
Description<-rep("c",length(File))

#(4) creem un data frame amb les quatre columnes que ens interessen i ho gravem a csv
lFile<-data.frame(File,Category,Subcateg,Description) 
head(lFile) 
write.csv2(lFile,file.path(workingDir,"lFile.csv"), row.names = FALSE) 

#(5) es modifica en excel el csv d'acord amb els arxius que vols incloure 
#i les categories que es poden veure més avall.
#s'agafa un arxiu d'un estudi anterior com a model
##simplifica les coses si l'arxiu antic s'ordena per la primera columna "File" així tots dos estaran igual
#recordar afegir \index.html en els QCDir
#Ordenar per Category
#guardar el fitxer separat per tabulacions

##definim els parametres de la funció
resultsDir<-"results/"
workingDir<-getwd()
htmlInfo <- list(To = "Ma Hernandez (VHIR). D1957",
                 Description = "Differential expression study in mice with amyloid beta protein blood vessels deposits in four time points",
                 Analysts = "Ricardo Gonzalo, Ferran Briansó and Alex Sanchez",           
                 Contact = "ueb@vhir.org")
linksFileName<-"lFile.csv"

#es defineixen les categories que volem que apareguin al results files
Category <-c("Report and result summaries","Data", "Quality Control","Analysis","Annotations",
             "Multiple Comparison","GO Analysis","Reactome")
#es defineixen les categories que hem posat al csv
Names<-c("INFO", "DATA","QC", "ANALYSIS", 
         "ANNOT", "MULTCOMP",  "GO", "REACTOME")

LinksFile2Html(linksFileName, workingDir, htmlInfo,categs.descs = Category,
                IndexDir = resultsDir, UEB = TRUE, categs.names = Names,
               resultsFileName="ResultFiles")

