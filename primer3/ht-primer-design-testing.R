# Testing ht-primer design using callPrimer3 function
# 
# First I need a sequence
library("Biostrings")
?readDNAStringSet
ydr114c.seq = readDNAStringSet(filepath = "ydr117c.fasta")
ydr114c.seq
?DNAString
ydr114c.seqA = DNAString(ydr114c.seq[[1]])
ydr114c.seqA

# # Load the function
source("/Users/nmcglincy/Documents/computing/github/ingolia_lab/ht-primer-design/callPrimer3-NJM.R")
ls()

foo = .callP3NreadOrg(seq = ydr114c.seqA, 
                      Tm = c(58, 60, 62),
                      name = "ydr117c")
foo

# I wonder if biostrings has a lapply like function for DNAstrings sets?

?DNAStringSet

tmaSeqs = readDNAStringSet(filepath = "tma-orf-sequencesA copy.fasta")
tmaSeqs
class(as.list(tmaSeqs))
str(as.list(tmaSeqs))
# try this with taking out my gene name addition
tma.primers = lapply(as.list(tmaSeqs), .callP3NreadOrg, Tm = c(58, 60, 62), name = "Tom")
str(tma.primers)
library(plyr)
tma.primers.df = ldply(tma.primers)
head(tma.primers.df)
write.csv(tma.primers.df, file = "tma-orf-primes.csv")

# Tweeking script to enable passing table to awk

?write.csv

source("pick-primers.R")
pick.primers
ls()
list.files()
pick.primers("tma-orf-sequencesA copy.fasta")
?system
??writelines
rm(list = ls())
