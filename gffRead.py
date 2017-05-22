#!/bin/python

import sys

def main(file):
    gff = readFile(file)
    
    # lijstje = []
    dictTest = {}


    for line in gff:
        lineSplit = line.split("\t")
        if len(lineSplit) > 1:
            # if lineSplit[2] in ["region", "sequence_feature"]:
            #     continue
            featureSplit = lineSplit[8].split(";")
            if lineSplit[2] == "gene":
                geneID = extractFeature(featureSplit, "ID").split("=")[1]
                geneName = extractFeature(featureSplit, "Name").split("=")[1]
                dictTest[geneID] = {"Name" : geneName, "Product": "NA", "Accession": "NA"}
                # print(extractFeature(featuresplit, "Name").split("=")[1])

            elif lineSplit[2] == "CDS":
                
                parent = extractFeature(featureSplit, "Parent").split("=")[1]
                product = extractFeature(featureSplit, "product").split("=")[1]
                name = extractFeature(featureSplit, "Name").split("=")[1]
                dictTest[parent]["Product"] = product
                dictTest[parent]["Accession"] = name

    print("GeneID\tGffAccession\tNCBIGeneName\tNCBIGeneProduct")
    
    for x in dictTest:
        print(x + "\t" + dictTest[x]["Accession"] + "\t" + dictTest[x]["Name"]
              + "\t" + dictTest[x]["Product"])


def extractFeature(features, match):
    for feature in features:
        if match in feature:
            return(feature)


def readFile(file):
    with open(file, "r") as f:
        content = f.read().split("\n")
    f.closed
    return content


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        print("Too few arguments...")