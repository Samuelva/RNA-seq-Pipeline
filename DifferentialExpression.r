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
eTest <- function(dispersion) {
  dispersion <- estimateCommonDisp(dispersion)
  dispersion <- estimateTrendedDisp(dispersion)
  dispersion <- estimateTagwiseDisp(dispersion)
  
  et <- exactTest(dispersion)
  
  makeResults(dispersion, et)
}

# generalized linear model likelihood ratiotest
glmLRTest <- function(dispersion, design.matrix) {
  dispersion <- calcDispersions(dispersion)
  
  fit <- glmFit(dispersion, design.matrix)
  lrt <- glmLRT(fit)
  
  makeResults(dispersion, lrt)
}


# generalized linear model quasi likelihood f-test
glmQLFT <- function(dispersion, design.matrix) {
  dispersion <- calcDispersions(dispersion)
  
  fit <- glmQLFit(dispersion, design.matrix)
  qlf <- glmQLFTest(fit)
  
  makeResults(dispersion, qlf)
}


# Calculates the different types of dispersions for the glm tests.
calcDispersions <- function(dispersion) {
  dispersion <- estimateGLMCommonDisp(dispersion)
  dispersion <- estimateGLMTrendedDisp(dispersion)
  dispersion <- estimateGLMTagwiseDisp(dispersion)
  
  return(dispersion)
}


# Processes the results, formats it and writes it away.
makeResults <- function(dispersion, results) {
  # If a KEGG species identifier is supplied, pathway analysis will be performed,
  # if not, it is skipped.
  if (!is.null(species)) {
    # Object containing the NCBI protein ID's and the corresponding KEGG protein id.
    keggIDs <- keggConv(species, "ncbi-proteinid")
    pathwayAnnotation(results, keggIDs)
  }
  
  gffInfo <- loadCSV()
  
  countsPM <- cpm(dispersion)
  makeGraphs(cpm(dispersion, log=TRUE), results)
  fdr <- p.adjust(results$table$PValue, method="BH")
  results <- cbind(countsPM, results$table, fdr)
  results$Accession <- rownames(results)
  results <- merge(results, gffInfo, by.x="Accession", by.y="V2")
  results <- geneAnnotation(results)
  
  exportData(results, "Genes.csv")
}


loadCSV <- function() {
  gffInfo <- read.csv(gff, sep="\t", header=FALSE)
  return(gffInfo)
}


# Makes grapgs
makeGraphs <- function(cpm, results) {
  topDEG <- rownames(topTags(results, n=30))
  topDEGG <- unname(sapply(topDEG, function(x)
    {results$NCBIGeneName[results$Accession==x]}))
  png(paste(output, "Heatmap.png"))
  # Log cpm is manually calculated so it becomes available for every sample
  # instead of just the conditions (which is the case if we used 
  # results$table$logcpm)
  par(mar=c(4,4.5,2,1))
  par(oma=c(0,0,0,0))
  pheatmap(cpm[topDEG,], labels_col=results$samples$group) # log 2 value
  dev.off()
  
  results$colour[results$table$PValue<=0.05] <- "red"
  results$colour[results$table$PValue>0.05] <- "black"
  png(paste(output, "Volcano-plot.png"))
  # Log 10 value of pvalue is used (for now)
  plot(results$table$logFC, -log(results$table$PValue), xlab="log fold change", 
       ylab="-log p-value", main="Volcano plot", col=results$colour, cex=0.2)
  dev.off()
  
  de <- decideTestsDGE(results, adjust.method="fdr")
  detags <- rownames(results)[as.logical(de)]
  png(paste(output, "MA-plot.png"))
  # Plots log FC (log 10) against average log cpm (log 2)
  plotSmear(results, de.tags=detags)
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
  results$keggProteinID <- sapply(rownames(results), function(x) {
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
  write.csv2(results, paste(output, file, sep=""))
}


# Main function.
main <- function(args) {
  if (length(args) < 6) {
    cat("Too few arguments given, abortion...\n")
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
  
  dispersion <- estimateDisp(countsData, design.matrix)
  
  # Make if statement for different DGE methods
  if (method == "et" || method == "exacttest") {
    eTest(dispersion)
  } else if (method == "glmqlf" || method == "glmqlftest") {
    glmQLFT(dispersion, design.matrix)
  } else if (method == "glmlrt" || method == "qlmlrttest") {
    glmLRTest(dispersion, design.matrix)
  }
}

main(commandArgs(T))

