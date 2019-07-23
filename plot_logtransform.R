require(NOMAD)

## read in data
#data(BalfPeptides)

## convert to NOMAD inputs
#y <- BalfPeptides$Abundance
#x <- BalfPeptides[,c(1:4,6)]

#y1 <- addDat[,5]
#x1 <- addDat[,c(1:4,4)]
y <- addDat[,5]
x <- as.data.frame(addDat[,c(1:4)])

#y <- normDat[,5]
#x <- as.data.frame(normDat[,c(1:4)])

## check log transform function
nomadCheckLogTransform(y, x, rawTrim=0.9)
