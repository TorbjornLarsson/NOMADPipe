require(NOMAD)

#inPath <- cat(getwd(), "/", sep="")
#BalfPeptidesSubset <- read.table(paste(inPath, "BalfPeptidesSubset.txt", sep=""),sep="\t")
#peptides <- read.table(paste(inPath, "peptides.txt", sep=""),sep="\t")
#proteinGroups <- read.table(paste(inPath, "proteinGroups.txt", sep=""),sep="\t")
#Peptides <- read.table(paste(inPath, "Peptides.txt", sep=""),sep="\t")

#inPath <- cat(getwd(), "/Original data/", sep="")
#peptides <- read.table(paste(inPath, "peptides.txt", sep=""), sep="\t")
# peptides <- read.table(paste(inPath, "FilteredPeptidesTest.txt", sep=""), sep="\t")

inPath <- cat(getwd(), "/", sep="")
peptides <- read.table(paste(inPath, "FilteredPeptides.txt", sep=""), sep="\t")

