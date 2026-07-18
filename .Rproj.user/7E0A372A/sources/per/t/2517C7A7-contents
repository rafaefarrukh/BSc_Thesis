###################################################################################################
# Step 3: VALIDATION OF PHB SYNTHESIS
###################################################################################################

# libraries
library(tidyverse)
library(ape)
library(rBLAST)

# work dir
getwd()

# IMPORT DATA ##################################

proteins <- read.FASTA("2_bakta/annotated_genome.faa", type = "AA")

annotation <- read_tsv(
  "2_bakta/annotated_genome.tsv", 
  col_names = c("query_id", "subject_id", "identity", "alignment_length", "mismatches", "gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score")
  )

embl <- read.gff("2_bakta/annotated_genome.gff3")

raw <- embl$attributes
split <- str_split(raw, "Dbxref=")
raw2 <- unlist(map(split, 2))
split2 <- str_split(raw2, ";gene=")
db <- unlist(map(split2, 1))
database <- str_replace(db, ",", "\n")

# SEARCH ########################################

# extract annotation scores
interest <- data.frame(rbind(
  cbind(proposed_gene = "phaA", annotation[which(annotation$query_id == "MPDKNC_03513"),]),
  cbind(proposed_gene = "phaB", annotation[which(annotation$query_id == "MPDKNC_04271"),][1,]),
  cbind(proposed_gene = "phaC", annotation[which(annotation$query_id == "MPDKNC_03777"),])
  ))

# extract protein name
interest$protein_name <- c(
  str_trim(str_sub(labels(genes)[which(str_detect(labels(genes), "MPDKNC_03513"))], 13)),
  str_trim(str_sub(labels(genes)[which(str_detect(labels(genes), "MPDKNC_04271"))], 13)),
  str_trim(str_sub(labels(genes)[which(str_detect(labels(genes), "MPDKNC_03777"))], 13))
)

# extract database IDs
interest$databases <- c(
  database[which(str_detect(embl$attributes, "MPDKNC_03513"))],
  database[which(str_detect(embl$attributes, "MPDKNC_04271"))],
  database[which(str_detect(embl$attributes, "MPDKNC_03777"))]
)

# extract gene IDs
interest$Gene <- c(
  unlist(str_split(embl$attributes[which(str_detect(embl$attributes, "MPDKNC_03513"))], "gene="))[2],
  unlist(str_split(embl$attributes[which(str_detect(embl$attributes, "MPDKNC_04271"))], "gene="))[2],
  unlist(str_split(embl$attributes[which(str_detect(embl$attributes, "MPDKNC_03777"))], "gene="))[2]
)

# EXPORT ########################################

# restructure
interest <- interest[,c(1,16,14,2:13,15)]
View(interest)

# export 
write.csv(interest, "3_validation/info.csv", row.names = FALSE)
