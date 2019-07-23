## Filter Peptides dataset for all abundances non-zero (active peptides) and write to file.

require(NOMAD)

mkRow <- function(row) {
  
  ## Filter on peptides not expressed in any sample.
  
  base <- 50*(row-1)+1

  # x matrix
  xpep <- x[base:(base+49),]

  # y column
  Abundance <- y[base:(base+49)]

  # Create x,y result matrix
  r <- cbind(xpep,Abundance)

  # If y column is all zeros, return null object
  if(any(r$Abundance == 0)) (r <- NULL)

  return(r)
}

mkFrameList <- function(nRow) {
  d <- sapply(seq_len(nRow),function(i) {
    ri <- mkRow(i)
    data.frame(ri)
  })
  do.call(rbind,d)
}

tryCatch({
  n <- 73913
  #n <- 10
  filtDat <- mkFrameList(n)
  
  #colnames(filtDat) <- c("Protein", "Peptide", "Run", "iTRAQ", "Abundance")
  inPath <- cat(getwd(), "/", sep="")
  
  write.table(filtDat, paste(inPath, "NonzeroAbundancePeptidesOnNOMADformat.txt", sep=""), sep="\t")
  
}, warning = function(w) {
  print(w)
}, error = function(e) {
  print(e)
})
