###################################################################################################
# DATA PREPARATION: GENOME ANNOTATION TABLE
###################################################################################################

# LIBRARIES #####################################

library(tidyverse)
library(ape)

# IMPORT DATA ###################################

# tsv file
tsv <- read_tsv(
  "2_bakta/annotated_genome.tsv", 
  col_names = c("query_id", "subject_id", "identity", "alignment_length", "mismatches", "gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score")
)

# gff3 file
gff <- read.gff("2_bakta/annotated_genome.gff3")

# EXTRACT ATTRIBUTES ############################

# attributes 
raw <- gff$attributes

# convert to dataframe
df <- tibble(raw) %>%
  rowid_to_column("row_id") %>%
  separate_longer_delim(raw, delim = ";") %>%
  separate_wider_delim(raw, delim = "=", names = c("key", "value"), too_many = "merge") %>%
  pivot_wider(names_from = key, values_from = value) %>%
  select(-row_id)

# split databases
df <- df %>%
  mutate(row_id = row_number()) %>%
  separate_longer_delim(Dbxref, delim = ",") %>%
  separate_wider_delim(Dbxref, delim = ":", names = c("db_name", "db_val"), too_many = "merge") %>%
  group_by(row_id, db_name) %>%
  summarise(db_val = paste(db_val, collapse = ","), .groups = "drop") %>%
  pivot_wider(names_from = db_name, values_from = db_val) %>%
  right_join(df %>% mutate(row_id = row_number()), by = "row_id") %>%
  select(-Dbxref, -row_id)

# merge to get complete annotation table
df <- left_join(df, tsv, by = c("ID" = "query_id"))

# DATA FOR DASHBOARD ############################

overview <- df %>%
  select(ID, gene, Name, product)

databases <- df %>%
  select(ID, BlastRules, COG, GO, KEGG, SO, UniRef, EC, RefSeq, UniParc, RFAM, NCBIFam, PFAM)

blast <- df %>%
  select(ID, subject_id, query_start, query_end, subject_start, subject_end, identity, alignment_length, mismatches, gap_openings, evalue, bit_score)

# EXPORT DATA ###################################

saveRDS(
  list(
    overview = overview,
    databases = databases,
    blast = blast
  ),
  "dashboard/data/annotated_genome.rds"
)

# CLEAN ENVR ####################################

rm(list=ls())

###################################################################################################