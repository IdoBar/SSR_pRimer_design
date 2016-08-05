# Install missing CRAN packages if needed

# git.packages <- c("docopt/docopt.R")
# github.packages <- c("argparser")
install.deps <- function(p, repo="cran"){
  call_install <- switch(repo,cran=c("install.packages(package, repos=\"http://cloud.r-project.org/\""),
                         bioc=c("biocLite(package, suppressUpdates=TRUE"),
                         git=c("install.github(package"))
  if (repo=="bioc") eval(parse(text = getURL("http://bioconductor.org/biocLite.R", ssl.verifypeer=FALSE)))
  for (package in p) {
    if (!package %in% installed.packages()) {
      cat(sprintf("Please wait, installing and loading missing package %s and its dependencies from %s.\n",
                  package, repo), file=stderr())
      suppressWarnings(eval(parse(text = sprintf("%s, quiet = TRUE)",call_install))))
      if (!package %in% installed.packages()) {
        ifelse(!interactive(),q("no", 2, TRUE),
               stop(sprintf("Unable to install missing package %s from %s.\n",
                            package, repo), call. = FALSE))
      }
    }
    require(package, character.only=TRUE, quietly=TRUE, warn.conflicts=TRUE)
  }
}


# Install and load cran packages
CRAN_packages <- c("RCurl", "dplyr", "ggplot2")
install.deps(CRAN_packages)

# Install and load bioconductor packages
bioc_packages <- c("GenomicRanges")
install.deps(bioc_packages, repo="bioc")

# Reload dplyr (must be last package loaded)
detach("package:dplyr", unload=TRUE)
library(dplyr)

# Load Primer3 functions (make sure primer3_core.exe and settings path are configured)
source("callPrimer3-IB.R")

# Fetch transcript data from Trinotate.sqlite
Trinotate <- src_sqlite("I:/Academic/Bioinformatics/Projects/Paspaley/Mantle/Assembly/P_maxima_mantle_Trinotate.sqlite")
transcripts <- tbl(Trinotate, "Transcript") %>% collect()

# Load PSR results file and match the length of each sequence from the matching transcript table at the db
PSR_table <- read.delim("tests/PSR_polySSR_2016-07-09_15_33_45.txt") %>% mutate(len=sapply(transcripts$sequence[match(.$Seq_ID, transcripts$transcript_id)], nchar, USE.NAMES=FALSE), dist_from_end=len-stop)

# define minimum flanking region
distFromEnd <- 50



# Apply filters (minimum flanking region length, remove single nucleotide repeats, remove non-ploymorphic, remove missing calls)
filtered_PSR <- PSR_table %>% filter(dist_from_end>=distFromEnd, start>=distFromEnd, !apply(.[grepl("poly_Family",colnames(.))], 1, function(r) "0" %in% r), apply(.[grepl("poly_Family",colnames(.))], 1, function(r) length(unique(r))>1), !grepl("\\([ACGT]\\)", SSR))

# Try to get primers for the first SSR
SSR1 <- filtered_PSR[1, c(1:4, 13:14)]
target <- paste(SSR1$start, (SSR1$stop-SSR1$start), sep=",")

# A test to make sure the function is working for the first microsatellite
SSR1_primers <- .callP3NreadOrg(seq = transcripts$sequence[match(SSR1$Seq_ID, transcripts$transcript_id)],
                name = SSR1$Seq_ID, size_range = '100-350', sequence_target = target, report = file.path("tests/output/", sprintf("%s_primer3_report_%s.txt", SSR1$Seq_ID, format(Sys.Date(), "%d_%m_%y"))))

# Then to all the rest

primers=mapply(.callP3NreadOrg,seq=transcripts$sequence[match(filtered_PSR$Seq_ID, transcripts$transcript_id)],
               name = filtered_PSR$Seq_ID, report = paste0(filtered_PSR$Seq_ID, sprintf("_primer3_report_%s.txt", format(Sys.Date(), "%d_%m_%y"))),
               sequence_target = paste(filtered_PSR$start, (filtered_PSR$stop-filtered_PSR$start), sep=","),
                     MoreArgs=list(size_range='100-350',drive = "I")
                     ,SIMPLIFY = FALSE,USE.NAMES = FALSE)
# USE.NAMES = FALSE

# Combine all results to one table
primers_table <- bind_rows(primers[!is.na(primers)]) %>% mutate(SSR=filtered_PSR$SSR[match(Seq_ID, filtered_PSR$Seq_ID)])
write.table(primers_table, sprintf("primer3_report_%s.txt", format(Sys.Date(), "%d_%m_%y")), quote = FALSE, sep = "\t",row.names = FALSE)


