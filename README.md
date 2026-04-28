# BrSeqTB — A pipeline for antimicrobial resistance inference from Mycobacterium tuberculosis WGS

**BrSeqTB** is a fully modular Nextflow DSL2 pipeline designed for comprehensive Mycobacterium tuberculosis whole-genome sequencing (WGS) analysis. Starting from Illumina paired-end FASTQ files, it produces analysis-ready outputs through an end-to-end workflow that includes:

- **Quality control and reporting**
- **Taxonomic contamination screening (Kaiju)**
- **Reference alignment (BWA-MEM)**
- **Multi-caller variant detection (GATK HaplotypeCaller, LoFreq, Delly)**
- **Variant functional annotation (SnpEff)**
- **Drug-resistance prediction based on the WHO catalogue**
- **Mixed infection inference**
- **Phylogenetic reconstruction and transmission network inference**
- **Variant summary (per sample and cohort-level)**
- **Clinical report**

BrSeqTB supports scalable execution (local, HPC, or cloud), environment isolation via Conda, and modular execution of individual workflow components. The pipeline generates standardized, clinically interpretable outputs, including integrated QC reports, resistance summaries, and cohort-level analyses — enabling robust genomic surveillance and research applications in tuberculosis.

## BrSeqTB Workflow DAG

```
INITIALIZATION
 ├─ KAIJU_DB
 ├─ OMS_CATALOG
 ├─ BWA_REF
 ├─ GATK_DICT
 └─ SNPEFF_DB
        ↓
MAKE_MANIFEST_VALIDATE
        ↓
SAMPLES (fan-out per biosample)

BLOCK 1 — PER BIOSAMPLE

SAMPLE
 ├─ FASTQC
 └─ TRIMMOMATIC
       ├─ KAIJU
       └─ BWA
            ├─ DELLY
            ├─ LOFREQ
            ├─ GATK_GVCF ── TBDR_RCOV
            ├─ GATK_VCF ── NORM ── LINEAGE
            ├─ LOFREQ + GATK_GVCF ── NTM_FILTER
            └─ LOFREQ + GATK_GVCF + GATK_VCF + NORM ── SNPEFF

BLOCK 2 — COHORT LEVEL

COHORT
   ↓
COHORT_FILTER
   ↓
SNP_MATRIX
   ├─ TRANSMISSION
   └─ IQTREE

BLOCK 3 — PER BIOSAMPLE

SAMPLE
 ├─ MIX_INFECTION
 ├─ RESISTANCE_TARGET
 │      ↓
 │  RESISTANCE_REPORT

BLOCK 4 — GLOBAL REPORTS

 ├─ RESISTANCE_SUMMARY
 └─ QC_SUMMARY

FINAL — PER BIOSAMPLE

SAMPLE
 └─ CLINICAL_REPORT

```

# BrSeqTB — Installation and Execution Guide

## Requirements

- **Java (OpenJDK 17)**
- **Nextflow (≥ 25.10.2)**
- **Micromamba** (used by Nextflow to create environments)

Linux or macOS is recommended.

---

# Installation

## 1️⃣ Install Java (OpenJDK 17)

Java is required to run Nextflow.

### Ubuntu/Debian

```bash
sudo apt update
sudo apt install openjdk-17-jdk -y
```

### macOS (Homebrew)

```bash
brew install openjdk@17
```

Verify installation:

```bash
java -version
```

You should see version 17.

---

## 2️⃣ Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

Verify installation:

```bash
nextflow -version
```

---

---

## 3️⃣ Install Conda (Miniconda Recommended)

BrSeqTB uses **Conda** to automatically create the software environment defined in:

```
envs/brseqtb.yml
```

We recommend installing **Miniconda**, a lightweight Conda distribution.

### Install Miniconda

Download and install:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### During Installation

When running the Miniconda installer:

- Accept the license agreement
- Accept the default installation path (`~/miniconda3`)
- When prompted:

  ```
  Proceed with initialization? [yes|no]
  ```

  Type **`yes`**.

  This step is important because it ensures that `conda` is added to your shell environment and is available to Nextflow.


Restart your terminal or run:

```bash
source ~/.bashrc
```

At this point, Conda is correctly configured and available system-wide.

---

### (Optional) Disable Automatic Base Activation

By default, Conda activates the `base` environment every time you open a new terminal.

If you prefer not to auto-activate `base`, run:

```bash
conda config --set auto_activate_base false
```

Then reload your shell:

```bash
source ~/.bashrc
```

This keeps Conda available to Nextflow while preventing automatic activation of the `base` environment.

---

### Verify Installation

```bash
which conda
conda --version
```

You should see the path to `miniconda3` and a Conda version number.

---

## ⚠ Important — Accept Conda Terms of Service

Recent versions of Conda require accepting the Anaconda channel Terms of Service before creating environments in non-interactive mode (such as when running Nextflow).

If you do not accept these terms, environment creation may fail.

Run the following commands **once per machine**:

```bash
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
```

This step enables automated environment creation by Nextflow.

---

## 4️⃣ Clone the Repository

```bash
git clone https://github.com/LaPAM-USP/BrSeqTB.git
cd BrSeqTB
bash install.sh
source ~/.bashrc
```

After installation, the brseqtb command is available in your system.
All executions must be run from the project directory containing your input/ and reads/ folders.

## Input Requirements

BrSeqTB requires two main inputs before execution:

#### 1️ `input/input_table.xlsx`

The input table **must contain a column with biosample IDs** that exactly match the FASTQ file prefixes.

* **Biosample ID** → **Required**
* **Clinical metadata (e.g., patient info, drug resistance, location, etc.)** → Optional

Clinical information is not mandatory for pipeline execution, but if provided, it will be incorporated into downstream reports (e.g., resistance summary and clinical report).

The file must be located at:

```
input/input_table.xlsx
```

---

#### 2️ FASTQ Files in `reads/` Directory

All sequencing reads must be placed inside the `reads/` directory:

```
reads/
```

Files must follow the **strict Illumina paired-end naming convention**:

```
biosample_S1_L001_R1_001.fastq.gz
biosample_S1_L001_R2_001.fastq.gz
```

Naming structure:

```
<BIOSAMPLE>_S<NUM>_L<NNN>_R1_001.fastq.gz
<BIOSAMPLE>_S<NUM>_L<NNN>_R2_001.fastq.gz
```

Explanation:

* `<BIOSAMPLE>` → Must exactly match the biosample ID in `input_table.xlsx`
* `S<NUM>` → Sample number
* `L<NNN>` → Lane number
* `R1/R2` → Read direction (forward/reverse)
* `_001` → File index
* `.fastq.gz` → Gzipped FASTQ format

⚠ The pipeline performs strict validation. Any mismatch between FASTQ filenames and biosample IDs in the input table will result in execution failure.

---

Correct structure example:

```
brseqtb/
├── input/
│   └── input_table.xlsx
├── reads/
│   ├── 1827-22_S1_L001_R1_001.fastq.gz
│   └── 1827-22_S1_L001_R2_001.fastq.gz
```


## Running

### Full Pipeline (Default Execution)
Runs the complete workflow from pre-processing to final clinical reports, including phylogeny and transmission.

```bash
brseqtb
```
### Run a Specific Module
Executes only a single module.

```bash
brseqtb --module <module_name>
```
Example:

```bash
brseqtb --module fastqc
brseqtb --module trimmomatic
```

> ⚠ **Important — Module Dependencies:**
> Although modules can be executed individually, they **respect the logical dependencies defined in the workflow DAG**.
>
> This means required input files must already exist before running a downstream module. For example:
>
> * Running `bwa` requires trimmed reads generated by `trimmomatic`.
> * Running `snpeff` requires completed variant calling steps (`lofreq`, `gatk_gvcf`, `gatk_vcf`, `norm`).
> * Running cohort-level modules requires completed per-sample variant processing.
>
> BrSeqTB does **not automatically execute upstream steps** when using `--module`. The user is responsible for ensuring that prerequisite outputs are available.


### Exclude Optional Modules
Some modules can be excluded during full pipeline execution:

```bash
brseqtb --exclude <module_name>
```
Example:

```bash
brseqtb --exclude kaiju
brseqtb --exclude transmission,iqtree
```
> ⚠ **Important:**
>
> * `--exclude` accepts multiple modules separated by commas **without spaces**.
> * `--module` accepts **only one module at a time** (do not use commas).
>
> ✔ Correct:
> `--exclude transmission,iqtree`
> `--module bwa`
>
> ✘ Incorrect:
> `--exclude transmission, iqtree`
> `--module bwa,lofreq`

---

## Parameters (summary)

| Parameter    |   Description            |
| --------     | ------------------------ |
| --inputTable | Path to the input sample sheet (input_table) used to generate the pipeline manifest. By default, BrSeqTB uses input/input_table.csv, but an equivalent .xlsx version can be provided for easier manual editing. The file must contain biosample identifiers matching the reads directory structure. |
| --readsDir | Directory containing input FASTQ files organized by biosample. Default: reads/. |
| --addKaijuManually     | If true, skips automatic Kaiju database download and expects the database to be already present locally in database/kaiju. |
| --auxCohort         | When activated, BrSeqTB supplements the user dataset with a predefined reference panel of ~10 high-quality TB-DR genomes during the cohort construction and variant filtering stages.     |
| --module          | Executes a single workflow module instead of the full pipeline. Accepts any valid module name (e.g., bwa, lofreq, cohort, clinical_report). Multiple values must be comma-separated.    |
| --exclude          | Excludes optional modules during full pipeline execution. Allowed values: kaiju, transmission, iqtree, clinical_report. Multiple values must be comma-separated.    |
| --profile standard          | Default execution profile. Runs locally using ~65% of available CPUs with dynamic parallelization and full Conda environment isolation.    |
| -profile lowmem          | Safe execution mode for low-resource machines (e.g., 8GB RAM laptops). Reduces CPU usage and memory allocation to prevent crashes.    |
| -profile hpc          | HPC execution profile using SLURM scheduler. Designed for cluster environments with higher CPU and memory allocation.    |

### Available --module Options

| Module Name          | Description                                           |
| -------------------- | ----------------------------------------------------- |
| `fastqc`             | Raw read quality control (FastQC).                    |
| `trimmomatic`        | Adapter trimming and read filtering.                  |
| `kaiju`              | Taxonomic contamination screening (Kaiju).            |
| `bwa`                | Reference alignment using BWA-MEM.                    |
| `delly`              | Structural variant detection (Delly).                 |
| `lofreq`             | Low-frequency variant detection (LoFreq).             |
| `gatk_gvcf`          | Variant calling in GVCF mode (GATK HaplotypeCaller).  |
| `gatk_vcf`           | Variant calling in VCF mode (GATK HaplotypeCaller)    |
| `norm`               | Variant normalization and processing.                 |
| `tbdr_rcov`          | Coverage calculation over TB drug-resistance regions. |
| `lineage`            | Lineage assignment based on variant data.             |
| `ntm_filter`         | Non-tuberculous mycobacteria filtering.               |
| `snpeff`             | Functional annotation of variants (SnpEff).           |
| `cohort`             | Cohort-level variant aggregation.                     |
| `cohort_filter`      | Cohort-level variant filtering.                       |
| `snp_matrix`         | SNP matrix generation for downstream analysis.        |
| `transmission`       | Transmission network inference.                       |
| `iqtree`             | Phylogenetic tree reconstruction (IQ-TREE).           |
| `mix_infection`      | Mixed infection inference per sample.                 |
| `resistance_target`  | Drug-resistance mutation identification per sample.   |
| `resistance_report`  | Per-sample resistance report generation.              |
| `resistance_summary` | Cohort-level resistance summary.                      |
| `qc_summary`         | Global QC summary report.                             |
| `clinical_report`    | Final per-sample clinical report generation.          |


### Available --exclude Options

| Module Name       | Description                                                                 |
| ----------------- | --------------------------------------------------------------------------- |
| `kaiju`           | Skips taxonomic contamination screening (Kaiju) during per-sample analysis. |
| `transmission`    | Skips transmission network inference at the cohort level.                   |
| `iqtree`          | Skips phylogenetic tree reconstruction (IQ-TREE).                           |
| `clinical_report` | Skips final per-sample clinical report generation.                          |

---

## Profiles Summary

| Profile  | Intended Use             | CPU Strategy         |
| -------- | ------------------------ | -------------------- |
| standard (default) | Local laptop/workstation | Auto-scaled (~65%)   |
| lowmem | Low-resource machines (e.g., 8GB RAM) | Conservative, limited parallelism   |
| hpc      | Cluster environments     | Scheduler-controlled |

---

## Environment Creation

On first execution, Nextflow will automatically create the Conda environment defined in `envs/brseqtb.yml`. The environment will be stored in `~/.nextflow_conda_cache/`. Subsequent runs will reuse the cached environment.


## Optional: Clean Environment and Work Directory

If needed:

```bash
nextflow clean -f
rm -rf ~/.nextflow_conda_cache
rm -rf work
```

Then rerun:

```bash
brseqtb
```

