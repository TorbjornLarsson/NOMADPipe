require(NOMAD)

## Read data
y <- addDat[,5]
x <- as.data.frame(addDat[,c(1:4)])

## check log transform function
nomadCheckLogTransform(y, x, rawTrim=0.9)
