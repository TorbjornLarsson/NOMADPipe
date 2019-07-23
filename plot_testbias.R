require(NOMAD)

# ## read in data
# data(BalfPeptides)
# 
# ## convert to NOMAD inputs
# y <- BalfPeptides$Abundance
# x <- BalfPeptides[,c(1:4,6)]

# 13e testet med filtrerade, all non-zero abundance, peptides
y <- filtDat$Abundance
x <- filtDat[,c(1:4)]

# Default
# apply NOMAD normalization
res <- nomadNormalization(y, x, doRobust=TRUE, doLog=TRUE)

# # 6e & 9e & 11e & 12e testet
#res <- nomadNormalization(y, x, doRobust=FALSE, doLog=TRUE)

# Default 2a & 3e & 6e & 11e & 12e testet
## apply protein assembly
protaov <- nomadAssembleProteins(y=res$y, x=res$x, method="TukeyBW",
                                 format="BySample", combineDupPeptides=FALSE,
                                 standardizeAbundance=FALSE)

# # 4e testet
# ## apply protein assembly
# protaov <- nomadAssembleProteins(y=res$y, x=res$x, method="TukeyBW",
#                                  format="BySample", combineDupPeptides=TRUE,
#                                  standardizeAbundance=FALSE)

# # 5e testet
# ## apply protein assembly
# protaov <- nomadAssembleProteins(y=res$y, x=res$x, method="TukeyBW",
#                                 format="BySample", combineDupPeptides=FALSE,
#                                 standardizeAbundance=TRUE)

# 7e & 9e testet
## apply protein assembly
# protaov <- nomadAssembleProteins(y=res$y, x=res$x, method="Mean",
#                                  format="BySample", combineDupPeptides=FALSE,
#                                  standardizeAbundance=FALSE)

# # 8e testet
# ## apply protein assembly
# protaov <- nomadAssembleProteins(y=res$y, x=res$x, method="Median",
#                                  format="BySample", combineDupPeptides=FALSE,
#                                  standardizeAbundance=FALSE)

# Default
## apply protein assembly to non-normalized data
protaovRaw <- nomadAssembleProteins(y=log2(y), x=x, method="TukeyBW",
                                    format="BySample", combineDupPeptides=FALSE,
                                    standardizeAbundance=FALSE)

# ## protaov and protaovRaw files are organized such that each column is
# ## one iTRAQ channel (8 iTRAQ channels for each of the 3
# ## runs, or days).
# runInd <- rep(1:3, each=8) ## look at bias due to run effect
# itraqInd <- rep(1:8, 3) ## look at bias due to iTRAQ channel

## Our samples have 5 runs and 10 tags
runInd <- rep(1:5, each=10) ## look at bias due to run effect
itraqInd <- rep(1:10, 5) ## look at bias due to iTRAQ channel

# ## Our normalized samples have 5 runs and 9 tags
# runInd <- rep(1:5, each=9) ## look at bias due to run effect
# itraqInd <- rep(1:9, 5) ## look at bias due to iTRAQ channel


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