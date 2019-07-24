require(NOMAD)

inPath <- cat(getwd(), "/", sep="")
peptides <- read.table(paste(inPath, "FilteredPeptides.txt", sep=""), sep="\t")
