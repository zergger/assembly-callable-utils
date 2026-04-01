# assembly-callable-utils

Portable wrappers for read-backed assembly refinement and callable-constrained assembly materialization.

This repository provides two representative custom scripts developed in the SJF workflow:

- `reapr_patch_nextpolish_ragtag.sh`: a conservative read-backed assembly-refinement wrapper
- `run_callable_assembly.sh`: a utility for materializing callable-constrained assembly views from scaffold FASTA and alignment-derived callable regions

These scripts were developed in the context of a short-read draft assembly workflow, but they are not species-specific and can be adapted to similar projects in other non-model systems.

## Scope

This repository is intentionally narrow.

It documents two reusable utility scripts:

1. one **assembly-refinement wrapper**
2. one **callable-region materialization utility**

It does **not** constitute a full public release of all manuscript-support scripts, intermediate analysis code, figure-generation notebooks, or the complete end-to-end workflow used in the associated study.

## Contents

```text
.
├── LICENSE
├── README.md
├── callable_assembly.env.example
├── reapr_patch_nextpolish_ragtag.sh
├── reapr_patch_nextpolish_ragtag.env.example
└── run_callable_assembly.sh
```

## Script 1: `reapr_patch_nextpolish_ragtag.sh`

This script implements a conservative, read-backed refinement path for a draft assembly.

### Main stages

1. **REAPR** breakpoint detection and breaking
2. **RagTag patch** using a donor assembly
3. **NextPolish** Illumina polishing by default (`POLCA` optional)
4. **RagTag correct** with read validation (optional; requires reference)
5. **RagTag reference-guided scaffolding** (optional; requires reference)

### Intended use

Use this script when you have:

- a primary draft assembly
- a donor assembly for patching
- paired-end short reads
- optionally, a scaffold- or chromosome-level reference assembly

### Key behavior

- **NextPolish** is the default polisher
- `POLCA` is optional
- RagTag correct/scaffold steps are **optional** and are only activated when `REFERENCE_ASM` is provided
- the script uses **resume-style behavior**: completed stages are skipped when their main outputs already exist
- logs are written under `WORKDIR/logs/`

### Typical inputs

- `PRIMARY_ASM`
- `DONOR_ASM`
- `R1`
- `R2`
- optional `REFERENCE_ASM`

### Typical outputs

Depending on which stages are enabled, outputs may include:

- REAPR-broken assembly
- RagTag-patched assembly
- NextPolish-polished assembly
- optional RagTag-correct outputs
- optional RagTag-scaffold outputs
- optional deduplicated outputs when `RUN_SORTNR=1`

### Minimal example

```bash
cp reapr_patch_nextpolish_ragtag.env.example reapr_patch_nextpolish_ragtag.env
$EDITOR reapr_patch_nextpolish_ragtag.env

set -a
source reapr_patch_nextpolish_ragtag.env
set +a

PRIMARY_ASM=/path/to/primary.fasta \
DONOR_ASM=/path/to/donor.fasta \
R1=/path/to/reads_R1.fastq.gz \
R2=/path/to/reads_R2.fastq.gz \
WORKDIR=/path/to/work_reapr_patch_nextpolish_ragtag \
THREADS=10 \
bash reapr_patch_nextpolish_ragtag.sh
```

### Example with reference-guided stages enabled

```bash
PRIMARY_ASM=/path/to/primary.fasta \
DONOR_ASM=/path/to/donor.fasta \
REFERENCE_ASM=/path/to/reference.fasta \
R1=/path/to/reads_R1.fastq.gz \
R2=/path/to/reads_R2.fastq.gz \
WORKDIR=/path/to/work_reapr_patch_nextpolish_ragtag_refA \
THREADS=10 \
bash reapr_patch_nextpolish_ragtag.sh
```

### Help

```bash
bash reapr_patch_nextpolish_ragtag.sh --help
```

---

## Script 2: `run_callable_assembly.sh`

This script materializes a callable-constrained assembly view from a scaffold FASTA plus either:

- a precomputed callable BED
- or a BAM-derived callable-region calculation

### Main modes

- `MODE=scaffold_bam`
- `MODE=corrected_query_agp`

### Intended use

Use this script when you want to:

- derive callable regions from mapped reads
- apply mapping-quality and depth thresholds
- lift corrected-query callable intervals onto scaffold coordinates through AGP
- generate masked and optionally callable-only FASTA views

### Key behavior

- accepts either a **precomputed callable BED** or a **BAM**
- can auto-calculate callable regions from `samtools depth`
- supports explicit or derived depth thresholds
- can emit:
  - callable BED
  - noncallable BED
  - masked scaffold FASTA
  - optional callable-only FASTA intervals

### Mode summary

#### `scaffold_bam`

Use when callable regions should be derived directly from scaffold-coordinate BAM depth.

#### `corrected_query_agp`

Use when callable regions are first defined in corrected-query coordinates and then lifted to scaffold coordinates through an AGP file.

This mode is useful in reference-guided branch workflows where the corrected-query BAM and AGP jointly define the callable scaffold-space product.

### Typical inputs

- `FASTA`
- `BAM` or `CALLABLE_BED`
- optional `AGP` in `corrected_query_agp` mode
- optional `RAGTAG_CORRECT_LOG` for auto-parsing depth center in corrected-query workflows

### Typical outputs

- callable BED
- noncallable BED
- masked FASTA
- optional callable-only FASTA intervals

### Minimal example with a config file

Copy and edit the shipped example config:

```bash
cp callable_assembly.env.example callable_assembly.env
$EDITOR callable_assembly.env
```

Example contents:

```bash
SAMPLE_TAG=sample1_refA
MODE=corrected_query_agp
FASTA=/path/to/ragtag.scaffold.fasta
BAM=/path/to/ragtag.correct.reads.s.bam
AGP=/path/to/ragtag.scaffold.agp
OUTDIR=/path/to/results/callable_assembly/sample1_refA
MQ_CUTOFF=30
DEPTH_CENTER=79
DEPTH_LOW=39
DEPTH_HIGH=158
MASK_CHAR=N
OUTPUT_CALLABLE_ONLY=0
```

Then run:

```bash
bash run_callable_assembly.sh callable_assembly.env
```

### Help

```bash
bash run_callable_assembly.sh --help
```

---

## Requirements

The exact tool set depends on which parts of the scripts you enable.

### Core requirements

- `bash`
- `samtools`
- `python3`

### Common requirements for `reapr_patch_nextpolish_ragtag.sh`

- `reapr`
- `ragtag.py`
- `nextpolish`
- optional `polca.sh`
- `seqtk`
- `samtools`
- one of:
  - `smalt`
  - `bwa-mem2`
  - `strobealign`
- `minimap2`

### Common requirements for `run_callable_assembly.sh`

- `samtools`
- `python3`
- standard Unix text-processing tools
- AGP support files when using `MODE=corrected_query_agp`

## Installation

Clone the repository and make the scripts executable if needed:

```bash
git clone https://github.com/<USER>/assembly-callable-utils.git
cd assembly-callable-utils
chmod +x *.sh
```

Install the external dependencies in your preferred environment manager.

The repository also ships two example config templates:

- `reapr_patch_nextpolish_ragtag.env.example`
- `callable_assembly.env.example`

## Design philosophy

These scripts reflect a conservative workflow style:

- prefer **read-backed** refinement over aggressive graph surgery
- keep **reference-guided** steps optional and explicit
- treat callable-region definition as a **materialized object**, not just a summary statistic
- preserve intermediate structure through stable work directories and logs

## What this repository is not

This repository is **not**:

- a full manuscript-reproduction archive
- a complete assembly benchmarking framework
- a packaged workflow manager
- a species-specific resource release
- a guarantee that the default example paths match your environment

You are expected to provide your own input paths and runtime environment.

## Reproducibility note

These scripts are wrappers around external tools. Reproducibility depends on:

- exact software versions
- consistent input files
- environment setup
- explicit parameter settings
- stable path conventions

For manuscript-grade reproduction, additional intermediate analysis code may be required.

## Citation

If you use these scripts in a publication, please cite:

- the associated manuscript. In preparation.
- the underlying external tools actually used in your run, including as applicable:
  - REAPR
  - RagTag
  - NextPolish
  - POLCA
  - samtools
  - minimap2
  - bwa-mem2 / smalt / strobealign

## License

This repository is released under the **MIT License**.

See `LICENSE`.
