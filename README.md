# BSc Thesis Code

##  Thesis Description

Title: _In-Silico_ Optimization of Polyhydroxyalkanoate Production From Gram Positive Bacteria Using Genome-Scale Metabolic Modelling.

Student Name: M Rafae Farrukh.

Supervisor: Dr. Samina Mehnaz.

Co-Supervisor: Dr. Rameez Hassan.

Institute: Forman Christian College University, Lahore, Pakistan.

Research Duration: 2025-26

<hr>

# Pipeline

## Step 1: Acquiring the Genome

We know that the bacteria of interest is *Rosellomora marisflavi*, and by searching it on NCBI Genome we found 45 genomes. We selected the reference genome for the study as it should best represent the species.

- **NCBI Genome:** https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=189381
- **Rosellomora marisflavi UTDF19-31A (ASM2217078v1):** https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_022170785.1/

---

## Step 2: Genome Annotation

RefSeq annotation is good, but BAKTA is better because it crosslinks to other databases like BRENDA, which is important for GEM construction.

There is a CLI, but it requires downloading a 32GB database, so the online version was used instead.

The built-in annotation table was used to check the presence of PHB-producing enzymes:

- **phaA:** MPDKNC_03513 and MPDKNC_04194 (paaJ) — acetyl-CoA C-acetyltransferase
- **phaB:** MPDKNC_04271 (fabG) — 3-hydroxybutyrate dehydrogenase
- **phaC:** MPDKNC_03777 (phaC) — alpha/beta hydrolase

All three (main) enzymes are present in the annotated genome.

- **Website:** https://bakta.computational.bio/job/eyJqb2JJRCI6IjZkMmFkYmY2LTlkZmUtNGFhNi05ZDYwLWU0YzBiZDFhMzIwZSIsInNlY3JldCI6ImdXVUg5c1pUMkV6S2lDdG1iYkFaanR4S05aSWpaQWhTIn0=
- **Input**: `sequence.fna` and the parameter settings are saved in `bakta_settings.png`
- **Output**: `2_bakta` folder contains the fna, ffn, faa, gff3, and other files

---

## Step 3: Validation of PHB Synthesis Genes

- **Step 3.1:** Information Extraction
- **Step 3.2:** ColabFold
- **Step 3.3:** Geometric Analysis

---

## Step 4: KBase

KBase provides a suite of tools for various bioinformatics purposes, such as constructing a GEM. It structures an experiment as a "narrative notebook" with "code cells," each cell representing a function.

- **Website:** https://narrative.kbase.us/narratives

### Cell 1: Import from Staging Area
- **FASTA file:** `sequence.fna`
- **GFF3 file:** `sequence.gff3` (the GFF3 file generated has the complete genome at the end of the file, which needs to be removed, otherwise import will fail)
- **Genome object name:** Genome
- **Genotype type:** Finished Isolate
- **Scientific name:** *Bacillus marisflavi* (old name for the bacteria)
- **Source of GFF3:** Other

### Cell 2: MS2 - Build Prokaryotic Metabolic Models with OMEGGA
- **Genome:** Genome
- **Suffix for Output Model:** `.xml`
- **Gapfill:** True
- **Change default media to complete:** True
- **Apply ATP check?:** True
- **Media to ensure ATP production:** Aerobic glucose
- **Template for reconstruction:** Gram positive
- **Gapfilling Type:** Sequential

The FBAModel was then exported as both SBML and TSV.

---

## Step 5: CarveMe

KBase is a bottom-up approach and builds the model by adding reactions. CarveMe is a top-down approach and builds the model by removing reactions from a template. By using KBase as soft constraints, we can make a better model using CarveMe.

Since KBase and CarveMe use different notations, we need to change the KBase TSV file so CarveMe recognizes them.

- **File:** `step_5.1.R`

To construct the model, we use CarveMe.

- **File:** `step_5.2.sh`
- **Output:** `model_1.xml` stored in the `models` folder

---

## Step 6: Curation

The draft GEM always has errors/mistakes which need to be fixed. After each step, MEMOTE is used to validate the model.

### Step 6.1: Model Modification
- PHB synthesis metabolites and reactions were checked and added.
- Use insights from experimental data to remove reactions, such as sugar utilization.

### Step 6.2: Fix Annotation
- Use BiGG database to update metabolite and reaction annotations.
- Use GFF3 file to update gene annotations.

### Step 6.3: Fix Unbound Flux
- Identify and fix unbound reactions.

---

## Step 7: Analyses

### Model Characterization
- 7.01 FBA
- 7.02 pFBA
- 7.06 Medium Optimization
- 7.09 Aeration and Carbon Source Screen

### PHB Pathway Characterization
- 7.03 Production Envelope
- 7.11 PHB Pathway
- 7.13 PHB Pathway Flexibility

### Parametric Sensitivity and Environmental Perturbations
- 7.07 Robustness Analysis - 1D
- 7.08 Robustness Analysis - 2D
- 7.04 PhPP
- 7.05 dFBA

### System Diagnostics and FVA
- 7.10 Shadow Price and Reduced Cost Analyses
- 7.12 Pathway Bottleneck Analysis
- 7.14 Flux Shift Analysis
- 7.15 Multi-fraction and Loopless FVA

### In Silico Strain Design
- 7.16 Essential Reaction Knockout Screening
- 7.17 Strain Design
