#!/usr/bin/env bash

###################################################################################################
# COMPLETE PIPELINE
###################################################################################################

# working directory
wd="/GEM/"
cd "$wd"

# create mamba envr
mamba create -n thesis_gem -c conda-forge -c bioconda \
  python=3.10 \
  pillow \
  tqdm \
  regex \
  platformdirs \
  sniffio \
  websocket-client \
  seaborn \
  networkx \
  SALib \
  highspy \
  ncbi-datasets-cli \
  carveme \
  rdkit \
  memote

# activate mamba
mamba activate thesis_gem

# install other python libraries
pip install cameo straindesign

# STEP 1-5 ######################################

# Step 1: Download Genome
1_genome/step_1.sh

# Step 2: Annotation with BAKTA (online)

# Step 3: Validation with BLAST and Domain Analysis with ColabFold
Rscript 3_validation/step_3.1.R
Jupyter step_3.2 (ColabFold).ipynb
python3 step_3.3.py 2>&1 | tee step_3.3.log

# Step 4: KBase Model (online)

# Step 5: CarveMe Model
## convert KBase model into soft constraints for carveme
Rscript 5_carveme/step_5.1.R
## use carveme to make model_1.xml
5_carveme/step_5.2.sh

#################################################

# STEP 6: CURATION ##############################

cd 6_curation

# model modifications
Rscript step_6.1.R
memote report snapshot "../models/model_2.xml" --filename "../memotes/memote_2.html"

# fix annotations
python3 step_6.2.1.py 2>&1 | tee step_6.2.1.log               # BIGG metabolite and reaction
python3 step_6.2.3.helper.py 2>&1 | tee step_6.2.3.helper.log # gene annotation helper
python3 step_6.2.3.py 2>&1 | tee step_6.2.3.log               # gene annotation
memote report snapshot "../models/model_3.xml" --filename "../memotes/memote_3.html"

# fix unbound flux
jupyter step_6.3.ipynb
memote report snapshot "../models/model_4.xml" --filename "../memotes/memote_4.html"

memote report snapshot "../models/model_5.xml" --filename "../memotes/memote_5.html"



cd ..

#################################################

# STEP 7: ANALYSES ##############################

cd 7_analyses

python3 step_7_config.py 2>&1 | tee step_7_config.log # setup

python3 step_7.01.py 2>&1 | tee step_7.01.log     # FBA
echo "==================================================================================================="
python3 step_7.02.py 2>&1 | tee step_7.02.log     # pFBA
echo "==================================================================================================="
python3 step_7.03.py 2>&1 | tee step_7.03.log     # Production Envelope
echo "==================================================================================================="
python3 step_7.04.py 2>&1 | tee step_7.04.log     # PhPP
echo "==================================================================================================="
python3 step_7.05.py 2>&1 | tee step_7.05.log     # dFBA (combine with maltose dfba)
python3 step_7.05_maltose.py 2>&1 | tee step_7.05_maltose.log     # dFBA (combine with maltose dfba)
echo "==================================================================================================="
python3 step_7.06.py 2>&1 | tee step_7.06.log     # Medium Optimization
echo "==================================================================================================="
python3 step_7.07.py 2>&1 | tee step_7.07.log     # Robustness Analysis - 1D
echo "==================================================================================================="
python3 step_7.08.py 2>&1 | tee step_7.08.log     # Robustness Analysis - 2D
echo "==================================================================================================="
python3 step_7.09.py 2>&1 | tee step_7.09.log     # Aeration and Carbon Source Screen
echo "==================================================================================================="
python3 step_7.10.py 2>&1 | tee step_7.10.log     # Shadow Price and Reduced Cost Analyses
echo "==================================================================================================="
python3 step_7.11.py 2>&1 | tee step_7.11.log     # PHB Pathway
echo "==================================================================================================="
python3 step_7.12.py 2>&1 | tee step_7.12.log     # Pathway Bottleneck Analysis
echo "==================================================================================================="
python3 step_7.13.py 2>&1 | tee step_7.13.log     # PHB Pathway Flexibility
echo "==================================================================================================="
python3 step_7.14.py 2>&1 | tee step_7.14.log     # Flux Shift Analysis
echo "==================================================================================================="
python3 step_7.15.py 2>&1 | tee step_7.15.log     # Multi-fraction and Loopless FVA
echo "==================================================================================================="
python3 step_7.16.py 2>&1 | tee step_7.16.log     # Essential Reaction Knockout Screening
echo "==================================================================================================="
python3 step_7.17.py 2>&1 | tee step_7.17.log     # Strain Design with Optknock
echo "==================================================================================================="

#################################################

# deactivate mamba
mamba deactivate

###################################################################################################
