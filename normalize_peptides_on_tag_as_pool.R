## Read in x,y and normalize peptide abundances per run in y on tag 10.
## Since tag 10 abundance can be 0, do nothing with that run.
## Reformat x, y and check data with a log plot.

mkRow <- function(row) {

  # Normalizing on pool of tags, such that run_n -> run_n*mean_tag/tag_n
  base <- 50*(row-1)+1
  
  # x matrix
  x01_09 <- x[base:(base+8),]
  x11_19 <- x[(base+10):(base+18),]
  x21_29 <- x[(base+20):(base+28),]
  x31_39 <- x[(base+30):(base+38),]
  x41_49 <- x[(base+40):(base+48),]
  
  # y normalization tags
  y10 <- y[(base+9)]
  y20 <- y[(base+19)]
  y30 <- y[(base+29)]
  y40 <- y[(base+39)]
  y50 <- y[(base+49)]
  
  # Tag pool
  y_pool <- (y10 + y20 + y30 + y40 + y50)/5 
  
  # If tag_n is zero, do nothing; else normalize
  if ( identical(y10, 0)) (y01_09 <- y[base:(base+8)]) else (y01_09 <- y[base:(base+8)]*(y_pool/y10))
  if ( identical(y20, 0)) (y11_19 <- y[(base+10):(base+18)]) else (y11_19 <- y[(base+10):(base+18)]*(y_pool/y20))
  if ( identical(y30, 0)) (y21_29 <- y[(base+20):(base+28)]) else (y21_29 <- y[(base+20):(base+28)]*(y_pool/y30))
  if ( identical(y40, 0)) (y31_39 <- y[(base+30):(base+38)]) else (y31_39 <- y[(base+30):(base+38)]*(y_pool/y40))
  if ( identical(y50, 0)) (y41_49 <- y[(base+40):(base+48)]) else (y41_49 <- y[(base+40):(base+48)]*(y_pool/y50))
  
  # Create result matrix
  r1 <- cbind(x01_09,y01_09)
  r2 <- cbind(x11_19,y11_19)
  r3 <- cbind(x21_29,y21_29)
  r4 <- cbind(x31_39,y31_39)
  r5 <- cbind(x41_49,y41_49)

  r <- rbind(as.matrix(r1),as.matrix(r2),as.matrix(r3),as.matrix(r4),as.matrix(r5) )

  return(r)
}

mkFrameList <- function(nRow) {
  d <- sapply(seq_len(nRow),function(i) {
    ri <- mkRow(i)
    data.frame(t(ri))
  })
  do.call(rbind,d)
}

tryCatch({
  n <- 73913
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