#!/bin/Rscript
# Differential expression script

suppressMessages(library(edgeR))
suppressMessages(library(KEGGREST))
suppressMessages(library(pheatmap))

method <- NULL
output <- NULL
species <- NULL
gff <- NULL

# Read the file with count files and its condition/group
loadData <- function(countMetaFile) {
  countMeta <- read.csv(countMetaFile, header=FALSE)
  return(countMeta)
}


# Loads the countdata, filters the counts and calculates the
# normalization factor for every sample.
getCountData <- function(countFiles, conditions) {
  counts <- readDGE(countFiles, group=conditions)
  keep <- rowSums(cpm(counts)>1) >= 2
  filteredCounts <- counts[keep, , keep.lib.sizes=FALSE]
  normFactor <- calcNormFactors(filteredCounts)
  return(normFactor)
}


# classic exact test
eTest <- function(countsData, design.matrix) {
  dispersion <- estimateDisp(countsData)
  
  et <- exactTest(dispersion)
  
  makeResults(dispersion, et)
}

# generalized linear model likelihood ratiotest
glmLRTest <- function(countsData, design.matrix) {
  dispersion <- estimateDisp(countsData, design.matrix)
  
  fit <- glmFit(dispersion, design.matrix)
  lrt <- glmLRT(fit)
  
  makeResults(dispersion, lrt)
}


# generalized linear model quasi likelihood f-test
glmQLFT <- function(countsData, design.matrix) {
  dispersion <- estimateGLMCommonDisp(countsData, design.matrix)
  dispersion <- estimateGLMTrendedDisp(countsData, design.matrix)
  
  fit <- glmQLFit(dispersion, design.matrix)
  qlf <- glmQLFTest(fit)
  
  makeResults(dispersion, qlf)
}


# Processes the results, formats it and writes it away.
makeResults <- function(dispersion, DGE) {
  # If a KEGG species identifier is supplied, pathway analysis will be performed,
  # if not, it is skipped.
  if (!is.null(species)) {
    # Object containing the NCBI protein ID's and the corresponding KEGG protein id.
    keggIDs <- keggConv(species, "ncbi-proteinid")
    pathwayAnnotation(DGE, keggIDs)
  }
  
  gffInfo <- loadCSV()
  
  countsPM <- cpm(dispersion)
  fdr <- p.adjust(DGE$table$PValue, method="BH")
  results <- cbind(countsPM, DGE$table, fdr)
  results$Accession <- rownames(results)
  results <- merge(results, gffInfo, by.x="Accession", by.y="GffAccession")
  results <- geneAnnotation(results)
  makeGraphs(cpm(dispersion, log=TRUE), DGE, results)
  
  exportData(results, "Genes.csv")
}


loadCSV <- function() {
  gffInfo <- read.csv2(gff, sep="\t")
  return(gffInfo)
}


# Makes grapgs
makeGraphs <- function(cpm, DGE, results) {
  topDEG <- rownames(topTags(DGE, n=30))
  GeneNames <- unname(sapply(topDEG, function(x)
  {gsub("_", ".", results$NCBIGeneName[results$Accession==x])}))
  png(paste(output, "Heatmap.png"))
  # Log cpm is manually calculated so it becomes available for every sample
  # instead of just the conditions (which is the case if we used 
  # DGE$table$logcpm)
  par(mar=c(4,4.5,2,2))
  par(oma=c(0,0,0,0))
  colLabels <- as.character(DGE$samples$group)
  pheatmap(cpm[topDEG,], labels_col=colLabels, labels_row=GeneNames) # log 2 value
  dev.off()
  
  DGE$colour[p.adjust(DGE$table$PValue, "BH")<=0.05] <- "red"
  DGE$colour[p.adjust(DGE$table$PValue, "BH")>0.05] <- "black"
  png(paste(output, "Volcano-plot.png"))
  # Log 10 value of pvalue is used (for now)
  plot(DGE$table$logFC, -log(p.adjust(DGE$table$PValue, "BH")), xlab="log2 fold change", 
       ylab="-log2 FDR", main="Volcano plot", col=DGE$colour, cex=0.2)
  dev.off()
  
  de <- decideTestsDGE(DGE, adjust.method="fdr")
  detags <- rownames(DGE)[as.logical(de)]
  png(paste(output, "MA-plot.png"))
  # Plots log FC (log 10) against average log cpm (log 2)
  plotSmear(DGE, de.tags=detags)
  dev.off()
}


# Annotate the genes.
# Function adds KEGG protein ID, protein function, KEGG pathway ID and the pathway
# name/function to each gene, if available.
geneAnnotation <- function(results) {
  keggIDs <- keggConv(species, "ncbi-proteinid")
  # Adds the KEGG protein ID's to the corresponding NCBI protein ID. The row names
  # in the supplied data frame are NCBI protein ID's, these are supplied to a function
  # which returns the corresponding KEGG protein ID.
  results$keggProteinID <- sapply(results$Accession, function(x) {
    paste("smu:", unlist(ncbiToKegg(strsplit(x, ".", fixed=TRUE)[[1]][1], keggIDs)), sep="")
  })
  
  keggProteins <- keggList(species)
  genePathwayIDs <- keggLink("pathway", species)
  speciesPathwayNames <- keggList(unique(unname(genePathwayIDs)))
  
  results$proteinFunction <- noName(results$keggProteinID, keggProteins)
  results$pathwayID <- noName(results$keggProteinID, genePathwayIDs)
  results$pathwayFunction <- noName(results$pathwayID, speciesPathwayNames)
  return(results)
}


# Function for converting NCBI protein ID's to KEGG protein ID's.
# Accepts a NCBI protein ID and a conversion object, containing NCBI protein ID's and their
# corresponding KEGG ID's
ncbiToKegg <- function(id, keggIDs) {
  proteinPrefix <- c("AP_", "NP_", "YP_", "XP_", "ZP_")
  
  # lapply(df, function(x) {
  if ("TRUE" %in% sapply(proteinPrefix, grepl, id)) {
    return(unlist(strsplit(keggIDs[paste("ncbi-proteinid:", id, sep="")][[1]], ":", fixed=TRUE))[2])
  } else {
    return(NA)
  }
}

# Used for conversion
# A source data frame column and target object is supplied.
# The target contains something like a list with ID's and corresponding data,
# pathway ID's and their function for example. A data frame column containing
# these ID's, for example, is supplied, which will be used to search the corresponding value
# to return and add to the original data frame.
noName <- function(df, target) {
  sapply(df, function(x) {
    if (!is.na(x)) {
      return(unname(target[x]))
    } else {
      return(NA)
    }
  })  
}


# Pathway annotation.
pathwayAnnotation <- function(resultsDf, keggIDs) {
  row.names(resultsDf) <- lapply(rownames(resultsDf), function(x) {
    unlist(ncbiToKegg(strsplit(x, ".", fixed=TRUE)[[1]][1], keggIDs))
  })
  
  kegg <- kegga(resultsDf, species.KEGG=species)
  exportData(kegg, "Pathways.csv")
}


# Write a data frame to a specified file.
exportData <- function(results, file) {
  write.csv2(results, paste(output, file, sep=""), dec=".")
}


# Main function.
main <- function(args) {
  if (length(args) < 6) {
    cat("Too few arguments given, aborting...\n")
    quit()
  }
  for (i in 1:length(args)) {
    if (args[i] %in% c("-f", "--file")) {
      countMetaFile <- args[i+1]
    } else if (args[i] %in% c("-m", "--method")) {
      method <- tolower(args[i+1])
    } else if (args[i] %in% c("-o", "--output")) {
      output <<- args[i+1]
    } else if (args[i] %in% c("-s", "--species")) {
      species <<- tolower(args[i+1])
    } else if (args[i] %in% c("-g", "--gff")) {
      gff <<- args[i+1]
    }
  }
  countMetaData <- loadData(countMetaFile)
  conditions <- as.vector(countMetaData$V1)
  countFiles <- as.vector(countMetaData$V2)
  
  countsData <- getCountData(countFiles, conditions)
  design.matrix <- model.matrix(~factor(conditions))
  
  # Make if statement for different DGE methods
  if (method == "et" || method == "exacttest") {
    eTest(countsData, design.matrix)
  } else if (method == "glmqlf" || method == "glmqlftest") {
    glmQLFT(countsData, design.matrix)
  } else if (method == "glmlrt" || method == "qlmlrttest") {
    glmLRTest(countsData, design.matrix)
  }
}

main(commandArgs(T))

