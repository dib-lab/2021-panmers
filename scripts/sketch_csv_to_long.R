## Get file names
files <- snakemake@input
files <- unlist(files, use.names=FALSE)

## read files into list with each row labelled by accession. 
tab_long <- list()
for(i in 1:length(files)){
  print(i)
  sig <- read.csv(files[i])           # read in signature csv as df
  acc <- colnames(sig)[1]             # set lib name using sample
  acc <- gsub("*.sig", "", acc)       # remove file suffix
  sig$acc <- acc                      # set libname as col
  colnames(sig) <- c("minhash", "acc")
  tab_long[[i]] <- sig
}

## bind into one dataframe
tab_long <- do.call(rbind, tab_long)
tab_long$present <- 1 # add a col with "1" for presence/absence
write.csv(tab_long, snakemake@output[['csv']],
          quote = F, row.names = F)
