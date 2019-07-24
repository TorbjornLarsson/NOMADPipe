## Format Peptides dataset.
## First without filtering on MaxQuant score 
## (which should be > 100, e.g. https://www.researchgate.net/topic/MaxQuant).

require(NOMAD)

## The peptide file is ready for parsing.
## It has already run iTRAQ correction on abundances; we don't need an index for that.
## We have filtered out rows with peptides identified to just one proteinGroup, or ~ 90 % of peptides.

## Parse file with unique peptide vs protein ID.
## Loop on rows in peptides with lapply for speed. For each protein copy the run m, tag channel n abundances 
## from "Reporter intensity corrected n m" cells. No score filtering yet. Hardcoded to save run time.

mkRow <- function(row) {
  r <- matrix(
    c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 1, as.numeric(as.character(peptides$V94[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 2, as.numeric(as.character(peptides$V95[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 3, as.numeric(as.character(peptides$V96[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 4, as.numeric(as.character(peptides$V97[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 5, as.numeric(as.character(peptides$V98[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 6, as.numeric(as.character(peptides$V99[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 7, as.numeric(as.character(peptides$V100[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 8, as.numeric(as.character(peptides$V101[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 9, as.numeric(as.character(peptides$V102[row])),as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 10, as.numeric(as.character(peptides$V103[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 1, as.numeric(as.character(peptides$V124[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 2, as.numeric(as.character(peptides$V125[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 3, as.numeric(as.character(peptides$V126[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 4, as.numeric(as.character(peptides$V127[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 5, as.numeric(as.character(peptides$V128[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 6, as.numeric(as.character(peptides$V129[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 7, as.numeric(as.character(peptides$V130[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 8, as.numeric(as.character(peptides$V131[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 9, as.numeric(as.character(peptides$V132[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 10, as.numeric(as.character(peptides$V133[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 1, as.numeric(as.character(peptides$V154[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 2, as.numeric(as.character(peptides$V155[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 3, as.numeric(as.character(peptides$V156[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 4, as.numeric(as.character(peptides$V157[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 5, as.numeric(as.character(peptides$V158[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 6, as.numeric(as.character(peptides$V159[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 7, as.numeric(as.character(peptides$V160[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 8, as.numeric(as.character(peptides$V161[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 9, as.numeric(as.character(peptides$V162[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 10, as.numeric(as.character(peptides$V163[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 1, as.numeric(as.character(peptides$V184[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 2, as.numeric(as.character(peptides$V185[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 3, as.numeric(as.character(peptides$V186[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 4, as.numeric(as.character(peptides$V187[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 5, as.numeric(as.character(peptides$V188[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 6, as.numeric(as.character(peptides$V189[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 7, as.numeric(as.character(peptides$V190[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 8, as.numeric(as.character(peptides$V191[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 9, as.numeric(as.character(peptides$V192[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 10, as.numeric(as.character(peptides$V193[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 1, as.numeric(as.character(peptides$V214[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 2, as.numeric(as.character(peptides$V215[row])),as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 3, as.numeric(as.character(peptides$V216[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 4, as.numeric(as.character(peptides$V217[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 5, as.numeric(as.character(peptides$V218[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 6, as.numeric(as.character(peptides$V219[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 7, as.numeric(as.character(peptides$V220[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 8, as.numeric(as.character(peptides$V221[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 9, as.numeric(as.character(peptides$V222[row])), as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 10, as.numeric(as.character(peptides$V223[row]))), nrow=5)
      return(r)
}



mkFrameList <- function(nRow) {
  d <- sapply(seq_len(nRow),function(i) {
    # Start in 2nd row in peptides
    ri <- mkRow(i+1)
    data.frame(ri)
  })
  do.call(rbind,d)
}

tryCatch({
  n <- 73913
  addDat <- mkFrameList(n)
  
  colnames(addDat) <- c("Protein", "Peptide", "Run", "iTRAQ", "Abundance")
  inPath <- cat(getwd(), "/", sep="")
  
  write.table(addDat, paste(inPath, "PeptidesOnNOMADformat.txt", sep=""), sep="\t")
  
}, warning = function(w) {
  print(w)
}, error = function(e) {
  print(e)
})
