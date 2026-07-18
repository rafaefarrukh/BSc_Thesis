#!/usr/bin/env bash

###################################################################################################
# STEP 1: ACQUIRING GENOME
###################################################################################################

# define acession ID
id="GCF_022170785.1"

# change work dir
cd 1_genome

# install genome assembly
datasets download genome accession "$id"

# unzip it
unzip "ncbi_dataset.zip"

# move genome to parent folder and rename it to "sequence"
mv ncbi_dataset/data/${id}*/*.fna ./sequence.fna

# change work dir to parent
cd ..

###################################################################################################
