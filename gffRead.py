#!/bin/python
# Last edit: 02-06-2017
# Script to parse the gene ID, accession, name and product for every gene in a GFF file
# and export it in csv format. The output is printed to the terminal so use ">" to redirect the output
# to a file.
# Usage: $ python gffRead.py gffFile.gff

import sys

def main(file):
    gff = readFile(file)
    
    # Stores all gene attributes
    attributeDict = {}


    for line in gff:
        lineSplit = line.split("\t")
        # Ignores lines with no tab delimitation as these are probably comments
        if len(lineSplit) > 1:
            # The 9th column of of the gff file (9th position in lineSplit list) contains attributes
            # such as the ID, name and product for a type (CDS, gene, exon, etc) 
            attributes = lineSplit[8].split(";")
            # Only features with a gene type will be taken from te GFF as these are the most important
            # to RNA-Seq analysis. From the gene types, the ID and name is extracted and put into a
            # dictionary which will later be expanded.
            if lineSplit[2] == "gene":
                geneID = extractAttribute(attributes, "ID").split("=")[1]
                geneName = extractAttribute(attributes, "Name").split("=")[1]
                attributeDict[geneID] = {"Name" : geneName, "Product": "NA", "Accession": "NA"}
            # Extracts parent, gene product and gene accession from the CDS type and adds these to the
            # previously made dictionary of the corresponding gene type.
            elif lineSplit[2] == "CDS":
                parent = extractAttribute(attributes, "Parent").split("=")[1]
                product = extractAttribute(attributes, "product").split("=")[1]
                accession = extractAttribute(attributes, "Name").split("=")[1]
                attributeDict[parent]["Product"] = product
                attributeDict[parent]["Accession"] = accession

    # Header for the csv file
    print("GeneID\tGffAccession\tNCBIGeneName\tNCBIGeneProduct")
    for gene in attributeDict:
        # Prints the attributes for every gene to the terminal 
        print(gene + "\t" + attributeDict[gene]["Accession"] + "\t" + attributeDict[gene]["Name"]
              + "\t" + attributeDict[gene]["Product"])


def extractAttribute(attributes, match):
    """
    Function that loops over the given list (attributes) and returns
    the item in the list that contains "match". Used to extract attributes from
    the attributes variable.
    """
    for attribute in attributes:
        if match in attribute:
            return(attribute)


def readFile(file):
    """
    Opens the given file and returns its content.
    """
    with open(file, "r") as f:
        content = f.read().split("\n")
    f.closed
    return content


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        print("Too few arguments...")