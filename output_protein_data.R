## Output Proteins dataset.

require(NOMAD)

## The log2 protein data is ready for output.
## It has already run iTRAQ and NOMAD correction on abundances.
## We have filtered out rows with peptides identified to just one proteinGroup, or ~ 90 % of peptides.

## Protein normalized data is an n*(exps*tags) matrix.

protDat <- protaov$scores
inPath <- cat(getwd(), "/", sep="")
write.table(protDat, paste(inPath, "NonzeroAbundancePeptidesNOMADProteinScores.txt", sep=""), sep="\t")
