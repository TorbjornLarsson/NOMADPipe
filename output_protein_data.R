## Output Proteins dataset.

require(NOMAD)

## The log2 protein data is ready for output.
## It has already run iTRAQ and NOMAD correction on abundances.
## We have filtered out rows with peptides identified to just one proteinGroup, or ~ 90 % of peptides.

## Protein normalized data is an n*(exps*tags) matrix.
#protDat <- protaov$scores[1:2,]

#rownonzeros <- rowSums(protaov$scores != 0)
#cat("Number of proteins not expressed in any sample: ", length(which(rownonzeros == 50)), "\n")
# cat("Number of proteins not expressed in any sample: ", length(which(rownonzeros == 45)), "\n")

protDat <- protaov$scores

inPath <- cat(getwd(), "/", sep="")
# write.table(protDat, paste(inPath, "NOMADNormalizedProteinScores.txt", sep=""), sep="\t")
write.table(protDat, paste(inPath, "NonzeroAbundancePeptidesNOMADProteinScores.txt", sep=""), sep="\t")

# newscores <- read.table(paste(inPath, "NOMADNormalizedProteinScores.txt", sep=""), sep="\t")
# inPath <- cat(getwd(), "/2a testet/", sep="")
# oldscores <- read.table(paste(inPath, "NOMADNormalizedProteinScores.txt", sep=""), sep="\t")