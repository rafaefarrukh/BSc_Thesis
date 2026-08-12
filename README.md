# In-Silico Work Pipeline for BSc Thesis

## Thesis Description

**Title:** _In-Silico_ Optimization of Polyhydroxyalkanoate Production From Gram Positive Bacteria Using Genome-Scale Metabolic Modelling.

**Description:** Polyhydroxyalkanoates (PHAs) are a class of biopolymers produced by bacteria which are being considered as the leading alternative to conventional plastic. Poly-3-hydroxybutyrate (PHB) is the most commonly produced bioplastic. The thesis screens several bacterial strains for PHB production and the selected strain was extensively characterized, including the optimization of growth conditions using response surface methodology (RSM). Additionally, the downstream recovery of PHB was optimized as well. Then using bioinformatics tools, a genome-scale metabolic (GEM) model was constructed and curated which describes all biochemical reactions taking place inside the bacterial strain. The model was studied and analyzed to better understand PHB biosynthesis and to provide suggestions for improving the yield. Finally, a cost analysis was performed to evaluate the commercial viability of PHB

**Duration:** 2025-26.

## Pipeline Description

The entire pipeline is orchestrated by a single master script, **`0_Pipeline.sh`**. Running this one script executes every step below **in order**, from environment setup through final analyses: it builds the mamba environment, then calls each step's script sequentially (`1_genome/step_1.sh` → `5_carveme/step_5.1.R` → `5_carveme/step_5.2.sh` → the `6_curation/` scripts → the `7_analyses/` scripts), `cd`-ing into each stage's folder as needed. Every Python step is piped through `tee` to save a `.log` file alongside its output, and `echo "==="` separators are printed between analysis steps so progress is easy to follow in the terminal. Steps 2, 3, and 4 (BAKTA, BLAST/domain validation, and KBase) are performed online/manually and are not called by the script — they sit between the automated steps.

The results are organized in a shiny dashboard which can be accessed [here](https://rafaefarrukh.shinyapps.io/bsc_thesis/).

---

## Step 0: Environment Setup

Before anything else, a dedicated `mamba` environment is created and activated, and the working directory is set to the project folder. This ensures every downstream script has the exact tool versions it needs (CarveMe, MEMOTE, COBRA-adjacent libraries, etc.).

| | |
|---|---|
| **Script** | `0_Pipeline.sh` (top section) |
| **Input** | None (environment specification is hard-coded in the script) |
| **Output** | Mamba environment `thesis_gem` containing: python=3.10, pillow, tqdm, regex, platformdirs, sniffio, websocket-client, seaborn, networkx, SALib, highspy, ncbi-datasets-cli, carveme, rdkit, memote — plus `cameo` and `straindesign` installed via pip |

```bash
wd=".../GEM/"
cd "$wd"

mamba create -n thesis_gem -c conda-forge -c bioconda \
  python=3.10 pillow tqdm regex platformdirs sniffio websocket-client \
  seaborn networkx SALib highspy ncbi-datasets-cli carveme rdkit memote

mamba activate thesis_gem
pip install cameo straindesign
```

---

## Step 1: Acquiring the Genome

We know that the bacteria of interest is *Rosellomora marisflavi*, and by searching it on NCBI Genome we found 45 genomes. We selected the reference genome for the study as it should best represent the species.

- **NCBI Genome:** https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=189381
- **Rosellomora marisflavi UTDF19-31A (ASM2217078v1):** https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_022170785.1/

| | |
|---|---|
| **Script** | `1_genome/step_1.sh` |
| **Input** | Sequence Accession ID |
| **Output** | `sequence.fna` (genome) |

---

## Step 2: Genome Annotation

RefSeq annotation is good, but BAKTA is better because it crosslinks to other databases like BRENDA, which is important for GEM construction.

There is a CLI, but it requires downloading a 32GB database, so the online version was used instead. This step is performed manually/online and is not called by `0_Pipeline.sh`.

The built-in annotation table was used to check the presence of PHB-producing enzymes:

- **phaA:** MPDKNC_03513 and MPDKNC_04194 (paaJ) — acetyl-CoA C-acetyltransferase
- **phaB:** MPDKNC_04271 (fabG) — 3-hydroxybutyrate dehydrogenase
- **phaC:** MPDKNC_03777 (phaC) — alpha/beta hydrolase

All three (main) enzymes are present in the annotated genome.

- **Website:** https://bakta.computational.bio/job/eyJqb2JJRCI6IjZkMmFkYmY2LTlkZmUtNGFhNi05ZDYwLWU0YzBiZDFhMzIwZSIsInNlY3JldCI6ImdXVUg5c1pUMkV6S2lDdG1iYkFaanR4S05aSWpaQWhTIn0=

| | |
|---|---|
| **Script** | None — performed via the BAKTA web interface |
| **Input** | `sequence.fna`; parameter settings saved in `bakta_settings.png` |
| **Output** | `2_bakta` folder containing the `.fna`, `.ffn`, `.faa`, `.gff3`, and other annotation files |

---

## Step 3: Validation of PHB Synthesis Genes

BLAST and domain/structural analysis of the phaA, phaB, and phaC sequences, run via three separate scripts (not invoked by `0_Pipeline.sh`, since ColabFold in particular is normally run interactively/on its own compute).

| Step | Script | Input | Output |
|---|---|---|---|
| 3.1 Information Extraction | `3_validation/step_3.1.R` | `2_bakta` output (`.faa`/`.ffn` sequences for phaA, phaB, phaC) | Extracted sequence/annotation information (e.g., BLAST hits, domain data) for the three genes |
| 3.2 ColabFold | `3_validation/step_3.2 (ColabFold).ipynb` | Sequences from Step 3.1 | Predicted 3D protein structures for phaA, phaB, phaC |
| 3.3 Geometric Analysis | `3_validation/step_3.3.py` | Structures from Step 3.2 | Geometric/structural validation results confirming enzyme identity |

---

## Step 4: KBase

KBase provides a suite of tools for various bioinformatics purposes, such as constructing a GEM. It structures an experiment as a "narrative notebook" with "code cells," each cell representing a function. Performed online and not invoked by `0_Pipeline.sh`.

- **Website:** https://narrative.kbase.us/narratives

| | |
|---|---|
| **Script** | None — performed via the KBase Narrative interface |
| **Input** | `sequence.fna`, `sequence.gff3` |
| **Output** | Draft FBAModel exported as SBML (`.xml`) and TSV |

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

### Step 5.1: Convert KBase Output to CarveMe Format
Since KBase and CarveMe use different notations, the KBase TSV file is converted so CarveMe recognizes it as soft constraints.

| | |
|---|---|
| **Script** | `5_carveme/step_5.1.R` |
| **Input** | KBase FBAModel TSV export |
| **Output** | Reformatted soft-constraints file for CarveMe |

### Step 5.2: Build the Draft Model
CarveMe is run using the reformatted KBase constraints to construct the draft model.

| | |
|---|---|
| **Script** | `5_carveme/step_5.2.sh` |
| **Input** | Soft-constraints file from Step 5.1 |
| **Output** | `model_1.xml`, stored in the `models` folder |

---

## Step 6: Curation

The draft GEM always has errors/mistakes which need to be fixed. After each modification, MEMOTE is used to validate the model and generate a snapshot report.

### Step 6.1: Model Modification
PHB synthesis metabolites and reactions were checked/added, and reactions unsupported by experimental data (e.g., sugar utilization) were removed.

| | |
|---|---|
| **Script** | `6_curation/step_6.1.R` |
| **Input** | `models/model_1.xml` |
| **Output** | `models/model_2.xml`; MEMOTE report `memotes/memote_2.html` |

### Step 6.2: Fix Annotation
Metabolite/reaction annotations are updated against the BiGG database, and gene annotations are updated using the GFF3 file.

| | |
|---|---|
| **Scripts** | `6_curation/step_6.2.1.py` (BiGG metabolite & reaction annotation, log: `step_6.2.1.log`) → `6_curation/step_6.2.3.helper.py` (gene annotation helper, log: `step_6.2.3.helper.log`) → `6_curation/step_6.2.3.py` (gene annotation, log: `step_6.2.3.log`) |
| **Input** | `models/model_2.xml`, GFF3 annotation file |
| **Output** | `models/model_3.xml`; MEMOTE report `memotes/memote_3.html` |

### Step 6.3: Fix Unbound Flux
Unbound reactions are identified and constrained.

| | |
|---|---|
| **Script** | `6_curation/step_6.3.ipynb` |
| **Input** | `models/model_3.xml` |
| **Output** | `models/model_4.xml`; MEMOTE report `memotes/memote_4.html` |

The final curated model, `models/model_5.xml`, is validated with a last MEMOTE snapshot (`memotes/memote_5.html`) and used as the input for all Step 7 analyses.

---

## Step 7: Analyses

Before any analysis, a shared configuration script sets up the environment (e.g., loading the final curated model and common parameters) that every subsequent `7_analyses/` script relies on.

| | |
|---|---|
| **Script** | `7_analyses/step_7_config.py` |
| **Input** | `models/model_5.xml` |
| **Output** | Shared configuration/session used by all Step 7 scripts; log: `step_7_config.log` |

Each analysis below is its own script, run in numerical order by `0_Pipeline.sh`, with console output piped to a matching `.log` file via `tee`.

### Model Characterization
| Step | Script | Input | Output |
|---|---|---|---|
| 7.01 FBA | `step_7.01.py` | `model_5.xml` | Flux distribution results; `step_7.01.log` |
| 7.02 pFBA | `step_7.02.py` | `model_5.xml` | Parsimonious flux distribution; `step_7.02.log` |
| 7.06 Medium Optimization | `step_7.06.py` | `model_5.xml` | Optimized medium composition; `step_7.06.log` |
| 7.09 Aeration and Carbon Source Screen | `step_7.09.py` | `model_5.xml` | Growth/production results across aeration & carbon sources; `step_7.09.log` |

### PHB Pathway Characterization
| Step | Script | Input | Output |
|---|---|---|---|
| 7.03 Production Envelope | `step_7.03.py` | `model_5.xml` | Growth-vs-PHB production envelope; `step_7.03.log` |
| 7.11 PHB Pathway | `step_7.11.py` | `model_5.xml` | PHB pathway flux analysis; `step_7.11.log` |
| 7.13 PHB Pathway Flexibility | `step_7.13.py` | `model_5.xml` | Pathway flexibility/variability results; `step_7.13.log` |

### Parametric Sensitivity and Environmental Perturbations
| Step | Script | Input | Output |
|---|---|---|---|
| 7.07 Robustness Analysis – 1D | `step_7.07.py` | `model_5.xml` | 1D robustness curves; `step_7.07.log` |
| 7.08 Robustness Analysis – 2D | `step_7.08.py` | `model_5.xml` | 2D robustness surfaces; `step_7.08.log` |
| 7.04 PhPP | `step_7.04.py` | `model_5.xml` | Phenotype Phase Plane results; `step_7.04.log` |
| 7.05 dFBA | `step_7.05.py` + `step_7.05_maltose.py` | `model_5.xml` | Dynamic FBA time-course results (glucose and maltose conditions combined); `step_7.05.log`, `step_7.05_maltose.log` |

### System Diagnostics and FVA
| Step | Script | Input | Output |
|---|---|---|---|
| 7.10 Shadow Price and Reduced Cost Analyses | `step_7.10.py` | `model_5.xml` | Shadow price / reduced cost tables; `step_7.10.log` |
| 7.12 Pathway Bottleneck Analysis | `step_7.12.py` | `model_5.xml` | Identified bottleneck reactions; `step_7.12.log` |
| 7.14 Flux Shift Analysis | `step_7.14.py` | `model_5.xml` | Flux shift comparisons; `step_7.14.log` |
| 7.15 Multi-fraction and Loopless FVA | `step_7.15.py` | `model_5.xml` | FVA flux ranges (multi-fraction, loopless); `step_7.15.log` |

### In Silico Strain Design
| Step | Script | Input | Output |
|---|---|---|---|
| 7.16 Essential Reaction Knockout Screening | `step_7.16.py` | `model_5.xml` | List of essential reactions; `step_7.16.log` |
| 7.17 Strain Design (OptKnock) | `step_7.17.py` | `model_5.xml` | OptKnock strain design candidates; `step_7.17.log` |

---

Once all Step 7 scripts have run, `0_Pipeline.sh` deactivates the `thesis_gem` mamba environment, completing the pipeline.
