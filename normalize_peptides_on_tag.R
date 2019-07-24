## Read in x,y and normalize peptide abundances per run in y on tag 10.
## Since tag 10 abundance can be 0, do nothing with that run.
## Reformat x, y and check data with a log plot.

mkRow <- function(row) {
  # Normalizing on tag
  base <- 10*(row-1)+1
  x9 <- x[base:(base+8),]
  if (identical(y[(base+9)], 0)) (y9 <- y[base:(base+8)]) else (y9 <- y[base:(base+8)]/y[(base+9)])
  r <- cbind(x9,y9)
  
  return(r)
}

mkFrameList <- function(nRow) {
  d <- sapply(seq_len(nRow),function(i) {
    ri <- mkRow(i)
    data.frame(t(ri))
  })
  #print(d)
  do.call(rbind,d)
}

tryCatch({
  n <- 73913*5
  normDat <- mkFrameList(n)
  colnames(normDat) <- c("Protein", "Peptide", "Run", "iTRAQ", "Abundance")
  inPath <- cat(getwd(), "/", sep="")

  write.table(normDat, paste(inPath, "NormedPeptidesOnNOMADformat.txt", sep=""), sep="\t")
  
}, warning = function(w) {
  print(w)
}, error = function(e) {
  print(e)
})


## check log transform function
nomadCheckLogTransform(y, x, rawTrim=0.9)