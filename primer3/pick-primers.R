# A script to pick potential primers from a multi-fasta file
# 
pick.primers = function(fastaFile) {
  library("Biostrings")
  dna.sequences = readDNAStringSet(filepath = fastaFile)
  #   
  # Load function calling Primer3
  source("/Users/nmcglincy/Documents/computing/github/ingolia_lab/ht-primer-design/callPrimer3-NJM.R")
  #
  # Specify values of:
  # @param size_range = '151-500'
  # @param Tm = c(58, 60, 62)
  l.primers = lapply(as.list(dna.sequences), .callP3NreadOrg)
  #
  # Converting into a data.frame I can write to a csv
  library(plyr)
  df.primers = ldply(l.primers)
  write.csv(df.primers, 
            file = "primer-table.csv", 
            row.names = FALSE,
            eol = "\n")
  #
  # Reformatting into a multi-fasta file that is easy to paste into BLAST etc
  system("awk -f table-to-fasta.awk primer-table.csv | tr -d '\"' > primers.fasta")
}
