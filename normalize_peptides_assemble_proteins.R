require(NOMAD)

## read in data
y <- addDat[,5]
x <- as.data.frame(addDat[,c(1:4)])

## default call
#results <- nomadNormalization(y, x)

# ## default call
#proteinScores <- nomadAssembleProteins(results$y, results$x)