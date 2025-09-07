library("dplyr")
library("readr")

# Parameters ----
args = commandArgs(trailingOnly=TRUE)
INPUT_DATA <- args[1]
OUT_PATH <- args[2]
if (is.na(OUT_PATH)){
  OUT_PATH <- ''
}
MIN_EVALUE_LOW <- 0.0001
MAX_HITS_LOW <- 5
MIN_EVALUE_HIGH <- 0.5
MAX_HITS_HIGH <- 20
LENGTH_CUTOFF <- 200


# file name template
filename <- basename(INPUT_DATA) # get the file name with path
filename <- tools::file_path_sans_ext(filename) # remove the file extension
out_file <- paste(OUT_PATH, filename, sep ='')
blast_filename <- paste(out_file,'blast.tsv',sep = '.')

# read datafiles ----
blast_result <- read_tsv(blast_filename, col_names = FALSE, show_col_types = FALSE)
df <- read.csv(INPUT_DATA, header = FALSE)

# set headers ----
colnames(blast_result) <- c('qseqid', 'sseqid', 'pident', 'align_length(aa)', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')

# group results by query - number of hits and minimal evalue ----
queries <- blast_result %>% group_by(qseqid) %>% summarise(hits_at_low = sum(evalue < MIN_EVALUE_LOW) ,hits_at_high = sum(evalue < MIN_EVALUE_HIGH), .groups = 'keep')
queries <- merge(queries, df[,c(2,8)], by.x = 'qseqid', by.y = 1) 
queries$V8 <- abs(queries$V8)

# quantile cutoff----
length_dropped <- df[df$V6-df$V5 > LENGTH_CUTOFF,]$V2

#filter results ----
dropped <- queries[queries$hits_at_low >= 5 | queries$hits_at_high >= 20,]$qseqid
dropped <- base::union(dropped,length_dropped)
true_cns <- queries[!queries$qseqid %in% dropped,]
dropped_df <- df[df$V2 %in% dropped,]
true_cns_df <- df[!(df$V2 %in% dropped),]

# write ouput files ----
write.table(dropped_df, paste(out_file,'.dropped.csv', sep= ""),sep = ',', col.names = F, row.names = F, quote = F)
write.table(true_cns_df, paste(out_file,'.non_ORFs.csv', sep= ""), sep = ',', col.names = F, row.names = F, quote = F)
