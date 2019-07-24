require(NOMAD)

## Filtered, all non-zero abundance, peptides
y <- filtDat$Abundance
x <- filtDat[,c(1:4)]

## Apply suggested NOMAD normalization; log2 transform for biological samples
res <- nomadNormalization(y, x, doRobust=TRUE, doLog=TRUE)

## Apply protein assembly
protaov <- nomadAssembleProteins(y=res$y, x=res$x, method="TukeyBW",
                                 format="BySample", combineDupPeptides=FALSE,
                                 standardizeAbundance=FALSE)

## Apply protein assembly to non-normalized data
protaovRaw <- nomadAssembleProteins(y=log2(y), x=x, method="TukeyBW",
                                    format="BySample", combineDupPeptides=FALSE,
                                    standardizeAbundance=FALSE)

## Our samples have 5 runs and 10 tags
runInd <- rep(1:5, each=10) ## look at bias due to run effect
itraqInd <- rep(1:10, 5) ## look at bias due to iTRAQ channel

## plot run and itraq biases for normalized and non-normalized data.
par(mfrow=c(2,2))
nomadCheckBias(dat=protaovRaw$scores, fact=runInd, numNAs=12,
               label="MS run: Raw")
nomadCheckBias(dat=protaov$scores, fact=runInd, numNAs=12,
               label="MS run: NOMAD")
nomadCheckBias(dat=protaovRaw$scores, fact=itraqInd, numNAs=12,
               label="iTRAQ: Raw")
nomadCheckBias(dat=protaov$scores, fact=itraqInd, numNAs=12,
               label="iTRAQ: NOMAD")