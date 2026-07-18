###################################################################################################
# STEP 4.1: CARVEME CONSTRAINTS
###################################################################################################

# library
library(tidyverse)

# work dir
getwd()

# read kbase reactions
kb_data <- read.delim("4_kbase/kbase_model_tsv/annotated_genome.xml-reactions.tsv", stringsAsFactors = FALSE)

# data manipulation
soft <- kb_data %>%
  filter(!is.na(bigg.id)) %>%
  mutate(
    clean_id = str_replace_all(bigg.id, "-", "_"), # change "-" --> "_"
    clean_id = str_remove_all(clean_id, "\\s*\\(.*?\\)"), # remove brackets, e.g. "PFK (atp") --> "PFK"
    clean_id = str_remove(clean_id, "_[0-9]+$"), # remove kbase numbers, e.g. "PGI_1" --> "PGI"
  ) %>%
  # mapping
  mutate(value = ifelse(direction == "<", -1, 1)) %>%
  select(clean_id, value) %>%
  distinct()

# check data for any missing info
soft[is.na(soft),] # NA values
soft[soft$clean_id == "", ] # missing ids
soft[!unlist(lapply(soft$value, is.numeric)),] # non numeric values

# clean accordingly
soft <- soft[-c(25,85),]

# write as tsv file
write.table(
  soft, 
  "5_carveme/soft_constraints.tsv",
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
  )

###################################################################################################