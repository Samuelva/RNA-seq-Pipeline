#!/bin/Rscript
# Differential gene expression analysis script.

suppressMessages(library(edgeR))
suppressMessages(library(KEGGREST))
suppressMessages(library(pheatmap))

# Global variables.
output <- NULL
species <- NULL
gff <- NULL

# Read the file with count files and its condition/group.
loadConditionsFile <- function(conditionsFileLoc) {
  conditionsFile <- read.csv(conditionsFileLoc, header=FALSE)
  return(conditionsFile)
}


# Loads the countdata, filters the counts and calculates the
# normalization factor for every sample.
loadCounts <- function(countFiles, readConditions) {
  counts <- readDGE(countFiles, group=readConditions)
  keep <- rowSums(cpm(counts)>1) >= 2
  filteredCounts <- counts[keep, , keep.lib.sizes=FALSE]
  normFactor <- calcNormFactors(filteredCounts)
  return(normFactor)
}


# Classic exact test.
eTest <- function(counts) {
  dispersion <- estimateDisp(counts)
  
  et <- exactTest(dispersion)
  
  generateResults(dispersion, et)
}


# Generalized linear model likelihood ratiotest.
glmLRTest <- function(counts, designMatrix) {
  dispersion <- estimateDisp(counts, designMatrix)
  
  fit <- glmFit(dispersion, designMatrix)
  lrt <- glmLRT(fit)
  
  generateResults(dispersion, lrt)
}


# Generalized linear model quasi likelihood f-test.
glmQLFT <- function(counts, designMatrix) {
  dispersion <- estimateGLMCommonDisp(counts, designMatrix)
  dispersion <- estimateGLMTrendedDisp(counts, designMatrix)
  
  fit <- glmQLFit(dispersion, designMatrix)
  qlf <- glmQLFTest(fit)
  
  generateResults(dispersion, qlf)
}


# Processes the results, formats it and writes it away.
generateResults <- function(dispersion, DGETestResults) {
  # If a KEGG species identifier is supplied, pathway analysis will be performed,
  # if not, it is skipped.
  if (!is.null(species)) {
    # Object containing the NCBI protein ID's and the corresponding KEGG protein id.
    keggIDs <- keggConv(species, "ncbi-proteinid")
    pathwayEnrichmentAnalysis(DGETestResults, keggIDs)
  }
  
  # Variable containing the contents of the gff file.
  gffInfo <- loadGFF()
  countsPM <- cpm(dispersion)
  fdr <- p.adjust(DGETestResults$table$PValue, method="BH")
  
  # Results is a dataframe containing information such as counts, countspm,
  # p-value, fdr, and annotations for every gene.
  results <- cbind(countsPM, DGETestResults$table, fdr)
  results$Accession <- rownames(results)
  results <- merge(results, gffInfo, by.x="Accession", by.y="GffAccession")
  results <- geneAnnotation(results)
  
  generateGraphs(cpm(dispersion, log=TRUE), DGETestResults, results)
  
  exportData(results, "Genes.csv")
}


# Opens a gff file and returns its content.
loadGFF <- function() {
  gffInfo <- read.csv2(gff, sep="\t")
  return(gffInfo)
}


# Generates graphs.
generateGraphs <- function(cpm, DGETestResults, results) {
  # determines the top differentially expressed genes with the built-in function topTags
  topDEGAccession <- rownames(topTags(DGETestResults, n=30))
  # Converts NCBI accession numbers to gene names which will be displayed in the
  # heatmap.
  topDEGNames <- unname(sapply(topDEGAccession, function(x)
  {gsub("_", ".", results$NCBIGeneName[results$Accession==x])}))
  
  # Heatmap generation.
  png(paste(output, "Heatmap.png", sep=""))
  par(mar=c(4,4.5,2,2))
  par(oma=c(0,0,0,0))
  colLabels <- as.character(DGETestResults$samples$group)
  pheatmap(cpm[topDEGAccession,], labels_col=colLabels, labels_row=topDEGNames) # log 2 value
  dev.off()
  
  # Volcano plot generation.
  DGETestResults$colour[p.adjust(DGETestResults$table$PValue, "BH")<=0.05] <- "red"
  DGETestResults$colour[p.adjust(DGETestResults$table$PValue, "BH")>0.05] <- "black"
  
  png(paste(output, "Volcano-plot.png", sep=""))
  plot(DGETestResults$table$logFC, -log(p.adjust(DGETestResults$table$PValue, "BH")), xlab="log2 fold change", 
       ylab="-log2 FDR", main="Volcano plot", col=DGETestResults$colour, cex=0.2)
  dev.off()
  
  # MA-plot generation.
  de <- decideTestsDGE(DGETestResults, adjust.method="fdr")
  detags <- rownames(DGETestResults)[as.logical(de)]
  png(paste(output, "MA-plot.png", sep=""))
  # Plots log FC (log 10) against average log cpm (log 2)
  plotSmear(DGETestResults, de.tags=detags)
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
    paste(paste(species, ":", sep=""), unlist(ncbiToKegg(strsplit(x, ".", fixed=TRUE)[[1]][1], keggIDs)), sep="")
  })
  
  keggProteins <- keggList(species)
  genePathwayIDs <- keggLink("pathway", species)
  speciesPathwayNames <- keggList(unique(unname(genePathwayIDs)))
  
  results$proteinFunction <- convert(results$keggProteinID, keggProteins)
  results$pathwayID <- convert(results$keggProteinID, genePathwayIDs)
  results$pathwayFunction <- convert(results$pathwayID, speciesPathwayNames)
  return(results)
}


# Function for converting NCBI protein ID's to KEGG protein ID's.
# Accepts a NCBI protein ID and a conversion object, containing NCBI protein ID's and their
# corresponding KEGG ID's.
ncbiToKegg <- function(id, keggIDs) {
  proteinPrefix <- c("AP_", "NP_", "YP_", "XP_", "ZP_")
  
  if ("TRUE" %in% sapply(proteinPrefix, grepl, id)) {
    return(unlist(strsplit(keggIDs[paste("ncbi-proteinid:", id, sep="")][[1]], ":", fixed=TRUE))[2])
  } else {
    return(NA)
  }
}


# Used for conversion.
# A source data frame column and target object is supplied.
# The target contains something like a list with ID's and corresponding data,
# e.g. pathway ID's and their function. A data frame column containing
# these pathway ID's, is supplied, which will be used to get the corresponding value
# (pathway function) which is returned and added to the original data frame.
convert <- function(df, target) {
  sapply(df, function(x) {
    if (!is.na(x)) {
      return(unname(target[x]))
    } else {
      return(NA)
    }
  })  
}


# Pathway enrichment analysis. Performed using the built in function kegga.
pathwayEnrichmentAnalysis <- function(resultsDf, keggIDs) {
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
      conditionsFileLoc <- args[i+1]
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
  
  conditionsFile <- loadConditionsFile(conditionsFileLoc)
  readConditions <- as.vector(conditionsFile$V1)
  countFiles <- as.vector(conditionsFile$V2)
  
  counts <- loadCounts(countFiles, readConditions)
  designMatrix <- model.matrix(~factor(readConditions))
  
  # Make if statement for different DGE methods.
  if (method == "et" || method == "exacttest") {
    eTest(counts)
  } else if (method == "glmqlf" || method == "glmqlftest") {
    glmQLFT(counts, designMatrix)
  } else if (method == "glmlrt" || method == "qlmlrttest") {
    glmLRTest(counts, designMatrix)
  }
}


main(commandArgs(T))
