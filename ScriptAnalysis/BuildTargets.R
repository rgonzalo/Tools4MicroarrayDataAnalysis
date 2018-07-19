#De moment el que fa aquest script es llegir el nom dels fitxers CEL i crear un dataframe amb ells
#la resta de columnes s'han de omplir de moment amb la fulla de càlcul

mainDir <-getwd()
workingDir <- mainDir
dataDir <- file.path(workingDir, "dades")
celDir <-  file.path(workingDir, "celfiles")

#listem els arxius i directoris de la carpeta celDir
FileName<-list.files(celDir,pattern=NULL,all.files=FALSE,recursive = FALSE,
                  full.names=FALSE,include.dirs = TRUE,no..=FALSE)

FileName

#creem la resta de variables que ha de tenir l'arxiu
Group <- rep("g",length(FileName))
ShortName <- rep("sn",length(FileName))
colores <- rep("c",length(FileName))
Batch <- rep("b", length(FileName))

#creem un data frame amb les columnes que ens interessen i ho gravem a csv
lFile <- data.frame(FileName, Group, ShortName, colores) 
head(lFile) 
client <- "MHG"
ID <- "D1957"
write.csv2(lFile, file.path(dataDir, paste("targets", client, ID, "csv", sep = ".")), 
                            row.names = FALSE)
