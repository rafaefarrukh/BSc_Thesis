#!/usr/bin/env bash

###################################################################################################
# STEP 5: CARVEME
###################################################################################################

# check work dir
cd "5_carveme"

# use carveme to make GEM
carve "../2_bakta/annotated_genome.faa" \
--soft soft_constraints.tsv \
--fbc2 \
-u grampos \
-o ../models/model_1.xml \
--debug > stdout.log 2> stderr.log

# work dir
cd ..

# get memote file (validation)
memote report snapshot "models/model_1.xml" --filename "memotes/memote_1.html"

###################################################################################################
