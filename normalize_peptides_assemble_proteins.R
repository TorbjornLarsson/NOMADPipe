require(NOMAD)

## Format data
y <- addDat[,5]
x <- as.data.frame(addDat[,c(1:4)])

## Default call
results <- nomadNormalization(y, x)

## Default call
proteinScores <- nomadAssembleProteins(results$y, results$x)