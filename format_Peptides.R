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


# mkRow <- function(j, nCol) {
#   print(j)
#   r <- c(j:(j+2))
#   print(r)
#   return(r)
# }



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


##########################

# row <- 2
# 
# protein <- as.numeric(as.character(peptides$V253[row])))
# print(protein)


# ## Preallocate the large and small data frames
# addDat <- data.frame(matrix(nrow=50,ncol=5, dimnames=list(c(), c("Protein", "Peptide", "Run", "iTRAQ", "Abundance"))))
# sumDat <- data.frame(matrix(nrow=73914,ncol=5, dimnames=list(c(1:73914), c("Protein", "Peptide", "Run", "iTRAQ", "Abundance"))))
# 
# tryCatch({
# 
#      #for(row in 2:dim(peptides)[1]) {
# 
#    for(row in 2:3) {
# 
#      #if (identical((row/100), floor(row/100))) print(row)
# 
  # # Run 1
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 1, as.numeric(as.character(peptides$V94[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 2, as.numeric(as.character(peptides$V95[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 3, as.numeric(as.character(peptides$V96[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 4, as.numeric(as.character(peptides$V97[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 5, as.numeric(as.character(peptides$V98[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 6, as.numeric(as.character(peptides$V99[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 7, as.numeric(as.character(peptides$V100[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 8, as.numeric(as.character(peptides$V101[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 9, as.numeric(as.character(peptides$V102[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 10, as.numeric(as.character(peptides$V103[row]))))
  # 
  # # Run 2
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 1, as.numeric(as.character(peptides$V124[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 2, as.numeric(as.character(peptides$V125[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 3, as.numeric(as.character(peptides$V126[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 4, as.numeric(as.character(peptides$V127[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 5, as.numeric(as.character(peptides$V128[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 6, as.numeric(as.character(peptides$V129[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 7, as.numeric(as.character(peptides$V130[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 8, as.numeric(as.character(peptides$V131[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 9, as.numeric(as.character(peptides$V132[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 10, as.numeric(as.character(peptides$V133[row]))))
  # 
  # # Run 3
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 1, as.numeric(as.character(peptides$V154[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 2, as.numeric(as.character(peptides$V155[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 3, as.numeric(as.character(peptides$V156[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 4, as.numeric(as.character(peptides$V157[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 5, as.numeric(as.character(peptides$V158[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 6, as.numeric(as.character(peptides$V159[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 7, as.numeric(as.character(peptides$V160[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 8, as.numeric(as.character(peptides$V161[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 9, as.numeric(as.character(peptides$V162[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 10, as.numeric(as.character(peptides$V163[row]))))
  # 
  # # Run 4
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 1, as.numeric(as.character(peptides$V184[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 2, as.numeric(as.character(peptides$V185[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 3, as.numeric(as.character(peptides$V186[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 4, as.numeric(as.character(peptides$V187[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 5, as.numeric(as.character(peptides$V188[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 6, as.numeric(as.character(peptides$V189[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 7, as.numeric(as.character(peptides$V190[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 8, as.numeric(as.character(peptides$V191[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 9, as.numeric(as.character(peptides$V192[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 10, as.numeric(as.character(peptides$V193[row]))))
  # 
  # # Run 5
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 1, as.numeric(as.character(peptides$V214[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 2, as.numeric(as.character(peptides$V215[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 3, as.numeric(as.character(peptides$V216[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 4, as.numeric(as.character(peptides$V217[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 5, as.numeric(as.character(peptides$V218[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 6, as.numeric(as.character(peptides$V219[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 7, as.numeric(as.character(peptides$V220[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 8, as.numeric(as.character(peptides$V221[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 9, as.numeric(as.character(peptides$V222[row]))))
  # addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 10, as.numeric(as.character(peptides$V223[row]))))
  # 
  
  
# }
# 
# colnames(addDat) <- c("Protein", "Peptide", "Run", "iTRAQ", "Abundance")
# inPath <- cat(getwd(), "/", sep="")
# 
# write.table(addDat, paste(inPath, "PeptidesOnNOMADformat.txt", sep=""),
#             sep="\t")
# 
# }, warning = function(w) {
#   print(w)
# }, error = function(e) {
#   print(e)
# })

#############################################################################
## Old file contained ~ 10 % peptides identified to several proteinGroup IDs.
## Note that find/replace has worked through here.

# # protgr_string <- as.character(peptides$V253[8])
# # protgroups <- as.list(as.numeric(unlist(strsplit(protgr_string, ";"))))
# # print(protgroups)
# # 
# # 
# # print(data.frame("Peptide" = c(1, 1)))
# #  sapply(protgroups, function(x) {
# #    print(x)
# #    a <- rbind(a, data.frame("Protein" = unlist(x)))
# #   a <- cbind(a, data.frame("Peptide" = c(1, 1)))
# #   b <<- a
# #   #a <<- rbind(a, x, as.numeric(as.character(peptides$V252[8])), 1, 1, as.numeric(as.character(peptides$V94[8])))
# #   #a <<- cbind(q, x, as.numeric(as.character(peptides$V252[8])), 1, 2, as.numeric(as.character(peptides$V95[8])))
# #
# # })
# 
# for(row in 2:dim(peptides)[1]) {
# #for(row in 2:10) {
# 
#   #row <- 8
# 
#   protgr_string <- as.character(peptides$V253[row])
#   protgroups <- as.list(as.numeric(unlist(strsplit(protgr_string, ";"))))
#   # print(protgroups)
# 
#   for (i in 1:length(protgroups)) {
#     # Run 1
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 1, as.numeric(as.character(peptides$V94[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 2, as.numeric(as.character(peptides$V95[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 3, as.numeric(as.character(peptides$V96[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 4, as.numeric(as.character(peptides$V97[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 5, as.numeric(as.character(peptides$V98[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 6, as.numeric(as.character(peptides$V99[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 7, as.numeric(as.character(peptides$V100[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 8, as.numeric(as.character(peptides$V101[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 9, as.numeric(as.character(peptides$V102[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 1, 10, as.numeric(as.character(peptides$V103[row]))))
#     
#     # Run 2
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 1, as.numeric(as.character(peptides$V124[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 2, as.numeric(as.character(peptides$V125[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 3, as.numeric(as.character(peptides$V126[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 4, as.numeric(as.character(peptides$V127[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 5, as.numeric(as.character(peptides$V128[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 6, as.numeric(as.character(peptides$V129[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 7, as.numeric(as.character(peptides$V130[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 8, as.numeric(as.character(peptides$V131[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 9, as.numeric(as.character(peptides$V132[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 2, 10, as.numeric(as.character(peptides$V133[row]))))
#     
#     # Run 3
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 1, as.numeric(as.character(peptides$V154[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 2, as.numeric(as.character(peptides$V155[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 3, as.numeric(as.character(peptides$V156[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 4, as.numeric(as.character(peptides$V157[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 5, as.numeric(as.character(peptides$V158[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 6, as.numeric(as.character(peptides$V159[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 7, as.numeric(as.character(peptides$V160[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 8, as.numeric(as.character(peptides$V161[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 9, as.numeric(as.character(peptides$V162[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 3, 10, as.numeric(as.character(peptides$V163[row]))))
#   
#     # Run 4
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 1, as.numeric(as.character(peptides$V184[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 2, as.numeric(as.character(peptides$V185[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 3, as.numeric(as.character(peptides$V186[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 4, as.numeric(as.character(peptides$V187[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 5, as.numeric(as.character(peptides$V188[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 6, as.numeric(as.character(peptides$V189[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 7, as.numeric(as.character(peptides$V190[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 8, as.numeric(as.character(peptides$V191[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 9, as.numeric(as.character(peptides$V192[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 4, 10, as.numeric(as.character(peptides$V193[row]))))
#   
#     # Run 5
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 1, as.numeric(as.character(peptides$V214[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 2, as.numeric(as.character(peptides$V215[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 3, as.numeric(as.character(peptides$V216[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 4, as.numeric(as.character(peptides$V217[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 5, as.numeric(as.character(peptides$V218[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 6, as.numeric(as.character(peptides$V219[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 7, as.numeric(as.character(peptides$V220[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 8, as.numeric(as.character(peptides$V221[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 9, as.numeric(as.character(peptides$V222[row]))))
#     addDat <- rbind(addDat, c(as.numeric(as.character(peptides$V253[row])), as.numeric(as.character(peptides$V252[row])), 5, 10, as.numeric(as.character(peptides$V223[row]))))
#   }
# }
# 
# colnames(addDat) <- c("Protein", "Peptide", "Run", "iTRAQ", "Abundance")
# inPath <- cat(getwd(), "/", sep="")
# 
# write.table(addDat, paste(inPath, "Peptides.txt", sep=""),
#             sep="\t")
# 

# 
#   ## Loop through peptides for the ones clustered to each protein group.
#   ## Loop through experiment runs.
#   ## Loop through tag channel.
# 
#   # Number of runs is 5
#   n <- 5
# 
#   # Start column for run is 94 + 30*(run-1).
#   # Run 1: V94
#   # Run 2: V124 ...
# 
#   # Peptides
#   for(i in 1:length(pep)) {
# #    print(i)
#     peptide <- pep[i]
#     #peptide <- 3
# #    print(peptide)
# 
#     # Experiment runs
#     for (j in 1:n) {
#       exp <- j
# #      print(exp)
# 
#       # Tag channel
#       for (k in 1:10) {
# #        print(k)
#         itraq <- k
#         #print(itraq)
#         abundances <- as.numeric(as.character(tmp2[unlist(peptide), (94-1) + 30*(j-1) + k]))
#         thisRow <- cbind(protein, peptide, exp, itraq, abundances)
#         sumDat <- rbind(sumDat, thisRow)
#       }
#     }
#   }
# 
# }
# colnames(sumDat) <- c("Protein", "Peptide", "Run", "iTRAQ", "Abundance")
# inPath <- cat(getwd(), "/", sep="")
# 
# write.table(sumDat, paste(inPath, "Peptides.txt", sep=""),
#             sep="\t")



#######################################
## Loop code: works, but inefficiently.

# ## Parse files with peptides and the proteinGroups they are identified to be in.
# ## Loop on experiment run m, then tag channel n from "Reporter intensity corrected n m" labels.
# 
# tmp1 <- proteinGroups
# tmp2 <- peptides
# 
# # ## filter peptides according to whatever rules one has. In
# # ## this case we require a confidence greater than 0.95 and
# # ## a "used" value of 1
# # rowInd <- as.numeric(tmp[,8]) >= 0.95 & as.numeric(tmp[,5]) == 1
# # 
# # ## keep only the rows and columns we want
# # filtDat <- tmp[rowInd, colInd]
# 
# ## Reformat so each row is a single peptide with columns defining MQ numerical
# ## protein group ID, peptide ID, run ID, tag ID, and peptide abundance.
# ## Inefficient re-formatting to show steps. Loop through
# ## rows and reformat. 10 channels for a single peptide.
# 
# # protein <- rep(tmp1[2,1], 10) ## protein label
# # print(protein)
# 
# #print(tmp1[2,244])
# #pepgr <- proteinGroups[2,244]
# #pepgr <- as.list(tmp1[2, 244])
# #pepgr <- as.array(tmp1[2, 244])
# 
# #protein <- tmp1$V243[2]
# #print(protein)
# #pepgr <- as.character(tmp1$V244[2])
# # print(pepgr)
# 
# sumDat <- NULL
# 
# ## First row is column names
# 
# for(row in 2:dim(tmp)[1]) {
# #for(row in 2:4) {
#   
# #  row <- 2
# #  print(row)
#   protein <- tmp1$V243[row]
#   pepgr <- as.character(tmp1$V244[row])
#   pep <- as.list(as.numeric(unlist(strsplit(pepgr, ";"))))
#   #print(pep[1])
#   
#   ## Loop through peptides for the ones clustered to each protein group.
#   ## Loop through experiment runs.
#   ## Loop through tag channel.
#   
#   # Number of runs is 5
#   n <- 5
#   
#   # Start column for run is 94 + 30*(run-1).
#   # Run 1: V94
#   # Run 2: V124 ...
#   
#   # Peptides
#   for(i in 1:length(pep)) {
# #    print(i)
#     peptide <- pep[i]
#     #peptide <- 3
# #    print(peptide)
#     
#     # Experiment runs
#     for (j in 1:n) {
#       exp <- j
# #      print(exp)
#       
#       # Tag channel
#       for (k in 1:10) {
# #        print(k)
#         itraq <- k
#         #print(itraq)
#         abundances <- as.numeric(as.character(tmp2[unlist(peptide), (94-1) + 30*(j-1) + k]))
#         thisRow <- cbind(protein, peptide, exp, itraq, abundances) 
#         sumDat <- rbind(sumDat, thisRow)
#       }
#     }
#   }
#  
# }
# colnames(sumDat) <- c("Protein", "Peptide", "Run", "iTRAQ", "Abundance")
# inPath <- cat(getwd(), "/", sep="")
# 
# write.table(sumDat, paste(inPath, "Peptides.txt", sep=""),
#             sep="\t")
# 
# #for(row in 1:dim(tmp)[1]) {
# #   iTRAQCount <- iTRAQCount + 1 # unique id for this row
# #   protein <- rep(filtDat[row,2], 8) ## protein label
# #   peptide <- rep(filtDat[row,3], 8) ## peptide label
# #   exp <- rep(i, 8) ## run label
# #   itraq <- 1:8 ## itraqs
# #   itraqCount <- rep(iTRAQCount, 8) ## peptide itraq ID
# #   ## abundance for each itraq (1 to 8)
# #   abundances <- as.numeric(filtDat[row, 5:12])
# #   thisRow <- cbind(protein, peptide, exp, itraq, itraqCount,
# #                    abundances)
# #   formDat <- rbind(formDat, thisRow)
# #} ## end for row
# 
# # ## append each file to final data variable
# # allDat <- rbind(allDat, formDat)
# # 
# # ## create column names. The column names must be labeled with the
# # ## following column names (except Abundance)
# # colnames(allDat) <- c("Protein", "Peptide", "Run", "iTRAQ",
# #                       "iTRAQCorrectionIndex", "Abundance")
# # 
# # inPath <- cat(getwd(), "/", sep="")
# # 
# # write.table(allDat, paste(inPath, "BalfPeptidesSubset.txt", sep=""),
# #             sep="\t")

