#Plot Venn diagramm from 2 to 5 comparisons
#Function based in previous work of Ferran Briansó and Miriam Mota
#Parameters available to tune:
##topTabs: list of toptables
##compNames: vector with the names of the comparisons to appear in the venn diagram
##label: Name to appear in the pdf file
##colFeat: column to seek the common objects. Usually "Gene.Symbol"
##colPVal: used adjustd pvalue or raw pvalue to select the objects to count
##pval: cutoff to use of pvalue
##pltR: show plot in R window
##pltPdf: save in a pdf file
##venn: plot venn diagram
##eul: plot euler diagram
##csv: create a csv file with the genes in each diagram
##colors: vector with the colors of each diagram
##trans: transparency of the diagrams
##cex1: size of the names of comparisons
##cex2: size of the numbers
##rotation: degrees to rotate the diagrams
##position: vector with the degrees of to write the names of comparisons. If two 0,-20 will be ok. If four c( 0, -20, -220, 150, 10)




createVennEuler2 <- function(topTabs, compNames, 
                            label = "selected", 
                            colFeat = "X", colPVal = "P.Value", pval = 0.05, 
                            pltR = TRUE, pltPdf = TRUE, venn = TRUE, eul = TRUE, csv = TRUE, 
                            colors = rainbow(length(compNames)), trans = 0.5, cex1 = 0.5, 
                            rotation = 0, position, cex2 =0.5){
  
## Initializing lists
list_genes_sel <- list()
  
## Reading input data
for (i in 1:length(topTabs)) {
    colpval <- which(names(topTabs[[i]]) == colPVal)
    colFeature <- which(names(topTabs[[i]]) == colFeat)
    list_genes_sel[[i]] <- as.character(topTabs[[i]][, colFeat][topTabs[[i]][, colpval] < pval])
  }
  
## Creating Venn Diagram
if (venn) {
  venn.plot <- venn.diagram(list_genes_sel,
                            category.names = compNames,
                            fill = colors, 
                            alpha = trans,
                            resolution = 800,
                            cat.cex = cex1,
                            cex = cex2,
                            main = paste0("Venn diagram (" , colPVal, " < ", pval,")"),
                            filename = NULL,
                            rotation.degree = rotation,
                            cat.pos = position)
  if (pltPdf) {
    pdf(paste0("VennDiagram.", label, ".", colPVal, pval, ".pdf"))
    grid.draw(venn.plot)
    dev.off()
    }
  if (pltR) {grid.draw(venn.plot)}
  }
  
  ## Creating Euler Diagram
  if (eul) {
    set <- NULL
    for (i in 1:length(compNames)) {
      set <- c(set, rep(compNames[i],length(list_genes_sel[[i]])))
    }
    v <- venneuler(data.frame(elements = c(unlist(list_genes_sel)),
                              sets = set))
    if (pltPdf) {
      pdf(paste0("EulerDiagram.", label, ".", colPVal, pval, ".pdf"))
      plot(v, main = paste0("Euler diagram (" , colPVal, " < ", pval,")"))
      dev.off()
    }
    if (pltR) {plot(v, main = paste0("Euler diagram (" , colPVal, " < ", pval,")"))}
  }
  
  ## Obtaining lists of combined elements and computing shared elements
  names(list_genes_sel) <- paste0(compNames, "_")
  combs <-  unlist(lapply(1:length(list_genes_sel),
                          function(j) combn(names(list_genes_sel), j, simplify = FALSE)),
                   recursive = FALSE)
  names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
  #str(combs)
  
  elements <- lapply(combs, function(i) Setdiff(list_genes_sel[i], list_genes_sel[setdiff(names(list_genes_sel), i)]))
  n.elements <- sapply(elements, length)
  list_res <- list(elements = elements, n.elements = n.elements)
  
  seq.max <- seq_len(max(n.elements))
  mat <- sapply(elements, "[", i = seq.max)
  mat[is.na(mat)] <- ""
  sharedElements <- rbind(t(data.frame(n.elements)),data.frame(mat))
  
  ## Writing table of shared elements as csv
  if (csv) {
    write.csv(sharedElements,
              file = paste0("sharedElements.", label, ".", colPVal, pval, ".csv"), 
              row.names = FALSE)
  }
  
  ## Returning table as data.frame
  return(sharedElements)
}

###################################################################
###### Internal functions used to extract the lists of shared elements 
###### both x and y are, in all cases, lists of character strings
###################################################################
Intersect <- function(x) {
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function(x) {
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function(x, y) {
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}
