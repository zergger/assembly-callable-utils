# assembly-callable-utils

Portable shell utilities for read-backed short-read assembly refinement, reference-guided projection, read remapping, and callable-constrained assembly materialization.

These scripts were developed in the context of a short-read draft assembly workflow, but they are not species-specific and can be adapted to similar projects in other non-model systems. Paths are supplied through environment variables or `.env` files.

## Scope

This repository provides five reusable utilities:

1. `reapr_patch_nextpolish.sh` — refine an Illumina-only draft assembly before projection.
2. `ragtag_project.sh` — project an already polished assembly onto one reference with RagTag correct/scaffold.
3. `map_for_callable.sh` — remap reads to one final FASTA and create a sorted/indexed BAM.
4. `run_callable_assembly.sh` — derive callable regions and write callable-masked FASTA outputs.
5. `ab_window_core_shell.sh` — classify A/B assembly windows and quantify coordinate core/shell callable space.

It is not a complete manuscript-support archive, workflow manager, figure-generation repository, or turnkey species-specific pipeline.

## Contents

```text
.
├── LICENSE
├── README.md
├── ab_window_core_shell.env.example
├── ab_window_core_shell.sh
├── callable_assembly.env.example
├── map_for_callable.env.example
├── map_for_callable.sh
├── quantify_core_shell.py
├── ragtag_project.env.example
├── ragtag_project.sh
├── reapr_patch_nextpolish.env.example
├── reapr_patch_nextpolish.sh
├── run_callable_assembly.sh
└── window_liftover_stats.py
```

## Recommended workflow

For one sample and two references, the intended order is:

```text
primary assembly + donor assembly + PE reads
  └── reapr_patch_nextpolish.sh
        -> polished unanchored assembly

polished assembly + reference A + PE reads
  └── ragtag_project.sh
        -> reference-A corrected/scaffold branch

polished assembly + reference B + PE reads
  └── ragtag_project.sh
        -> reference-B corrected/scaffold branch

branch FASTA + branch corrected-query BAM + AGP
  └── run_callable_assembly.sh MODE=corrected_query_agp
        -> primary branch callable BED and callable-masked FASTA

branch scaffold FASTA + PE reads
  └── map_for_callable.sh
        -> scaffold-coordinate remap BAM

branch scaffold FASTA + scaffold-coordinate remap BAM
  └── run_callable_assembly.sh MODE=scaffold_bam
        -> branch-specific callable BED and callable-masked FASTA

A/B branch FASTAs + A/B callable BEDs
  └── ab_window_core_shell.sh
        -> A/B window classes, callable core, callable shell, summary tables
```

The critical rule is that a BAM must be aligned to the same coordinate FASTA used by `run_callable_assembly.sh`, unless `MODE=corrected_query_agp` is used with a corrected-query BAM plus the corresponding RagTag AGP file.

`reapr_patch_nextpolish.sh` intentionally stops at the polished unanchored assembly. Reference-guided RagTag correct/scaffold is handled by `ragtag_project.sh`, which should be run separately for each reference branch.

## Script 1: `reapr_patch_nextpolish.sh`

This script implements a conservative read-backed refinement path for a draft assembly.

### Main stages

1. REAPR breakpoint detection and breaking
2. RagTag patch using a donor assembly
3. NextPolish Illumina polishing by default (`POLCA` optional)

This script does not run RagTag correct/scaffold and does not run post-polish redundancy or contamination filtering. Run `ragtag_project.sh` for reference projection, and use a separate local cleanup stage if `sortnr`, FCS-GX, or other post-polish filters are needed.

### Minimal use

```bash
cp reapr_patch_nextpolish.env.example reapr_patch_nextpolish.env
$EDITOR reapr_patch_nextpolish.env
bash reapr_patch_nextpolish.sh
```

### Typical direct invocation

```bash
PRIMARY_ASM=/path/to/primary.fasta \
DONOR_ASM=/path/to/donor.fasta \
R1=/path/to/reads_R1.fastq.gz \
R2=/path/to/reads_R2.fastq.gz \
WORKDIR=/path/to/work_reapr_patch_nextpolish \
THREADS=10 \
bash reapr_patch_nextpolish.sh
```

## Script 2: `ragtag_project.sh`

This script starts from an already polished/unanchored assembly and creates one reference-guided RagTag branch. It intentionally omits REAPR, patching, polishing, `sortnr`, and callable masking.

### Main outputs

- `WORKDIR/correct/ragtag.correct.fasta`
- `WORKDIR/correct/ragtag.correct.reads.s.bam` when read validation is enabled
- `WORKDIR/scaffold/ragtag.scaffold.fasta`
- `WORKDIR/scaffold/ragtag.scaffold.agp`
- `WORKDIR/ragtag_project.outputs.env`

### Minimal use

```bash
cp ragtag_project.env.example ragtag_project.refA.env
$EDITOR ragtag_project.refA.env
bash ragtag_project.sh ragtag_project.refA.env
```

Run once per reference, using a separate `WORKDIR` for each branch.

## Script 3: `map_for_callable.sh`

This script maps paired reads back to one assembly FASTA and writes the sorted/indexed BAM expected by `run_callable_assembly.sh`.

### Main outputs

- `OUTDIR/map/<SAMPLE_TAG>.sorted.bam`
- `OUTDIR/map/<SAMPLE_TAG>.sorted.bam.bai`
- `OUTDIR/stats/flagstat.txt`
- `OUTDIR/stats/samtools.stats.txt`
- `OUTDIR/map_for_callable.outputs.env`

### Minimal use

```bash
cp map_for_callable.env.example map_for_callable.refA.env
$EDITOR map_for_callable.refA.env
bash map_for_callable.sh map_for_callable.refA.env
```

Use this for scaffold-coordinate remap support or for any `MODE=scaffold_bam` callable materialization.

## Script 4: `run_callable_assembly.sh`

This script materializes a callable-constrained assembly view from a scaffold FASTA plus either a precomputed callable BED or a BAM-derived callable-region calculation.

### Modes

- `MODE=scaffold_bam`: derive callable regions directly from a BAM aligned to `FASTA`.
- `MODE=corrected_query_agp`: derive callable regions in corrected-query coordinates, then lift intervals to scaffold coordinates through AGP.

### Minimal use

```bash
cp callable_assembly.env.example callable_assembly.env
$EDITOR callable_assembly.env
bash run_callable_assembly.sh callable_assembly.env
```

### Corrected-query example

```bash
SAMPLE_TAG=sample_refA \
MODE=corrected_query_agp \
FASTA=/path/to/ragtag.scaffold.fasta \
BAM=/path/to/correct/ragtag.correct.reads.s.bam \
AGP=/path/to/scaffold/ragtag.scaffold.agp \
RAGTAG_CORRECT_LOG=/path/to/logs/ragtag_correct.log \
OUTDIR=/path/to/results/callable_assembly/sample_refA \
bash run_callable_assembly.sh
```

### Scaffold-BAM example

```bash
SAMPLE_TAG=sample_refA_scaffold \
MODE=scaffold_bam \
FASTA=/path/to/ragtag.scaffold.fasta \
BAM=/path/to/results/callable_mapping/sample_refA_scaffold/map/sample_refA_scaffold.sorted.bam \
OUTDIR=/path/to/results/callable_assembly/sample_refA_scaffold \
bash run_callable_assembly.sh
```

## Script 5: `ab_window_core_shell.sh`

This script treats the two projected branches as coordinate-specific views of the same sample. It aligns A against B and B against A, classifies fixed windows as one-to-one, one-to-many, or unmapped, and intersects those window classes with branch-specific callable BED files.

### Main outputs

- `OUTDIR/stats/A_to_B_window_classification_primary.tsv`
- `OUTDIR/stats/B_to_A_window_classification_primary.tsv`
- `COMPARE_OUT/tableS_window_primary.tsv`
- `COMPARE_OUT/tableS_callable.tsv`
- `COMPARE_OUT/tableS_core_shell.tsv`
- `OUTDIR/core_shell/bed/<reference>/*.bed`

### Minimal use

```bash
cp ab_window_core_shell.env.example ab_window_core_shell.env
$EDITOR ab_window_core_shell.env
bash ab_window_core_shell.sh ab_window_core_shell.env
```

The primary coordinate space is controlled by `PRIMARY_REGEX_A` and `PRIMARY_REGEX_B`. Use `.*` only when the branch FASTA files contain only the sequences that should contribute to the primary comparison.

## Requirements

The exact tool set depends on which scripts are used.

### Core

- `bash`
- `python3`
- `samtools`
- standard Unix tools (`sort`, `awk`, `grep`, `flock` where locking is used)

### `reapr_patch_nextpolish.sh`

- `reapr`
- `ragtag.py`
- `NextPolish` command-line executable (`nextPolish`)
- optional `polca.sh`
- `seqtk`
- `samtools`
- one REAPR mapper path: `smalt`, `bwa-mem2`, or `minibwa`
- `minimap2` for RagTag read validation when enabled

### `ragtag_project.sh`

- `ragtag.py`
- `samtools`
- `minimap2` when using `RAGTAG_READ_ALIGNER=minimap2`

### `map_for_callable.sh`

- `samtools`
- `bwa-mem2` or `bwa`
- `python3` when optional depth summaries are enabled

### `run_callable_assembly.sh`

- `samtools`
- `bedtools`
- `python3`

### `ab_window_core_shell.sh`

- `samtools`
- `minimap2`
- `python3`

## Installation

```bash
git clone https://github.com/zergger/assembly-callable-utils.git
cd assembly-callable-utils
chmod +x *.sh
```

Install external dependencies in your preferred environment manager. The scripts resolve tools from `PATH` unless an absolute path is supplied in a config file.

## Validation checks

```bash
bash -n reapr_patch_nextpolish.sh
bash -n ragtag_project.sh
bash -n map_for_callable.sh
bash -n run_callable_assembly.sh
bash -n ab_window_core_shell.sh

bash reapr_patch_nextpolish.sh --help
bash ragtag_project.sh --help
bash map_for_callable.sh --help
bash run_callable_assembly.sh --help
bash ab_window_core_shell.sh --help

python3 -m py_compile window_liftover_stats.py quantify_core_shell.py
```

## Design philosophy

- Keep reference-guided projection explicit and per-reference.
- Keep BAM generation separate from callable-mask materialization.
- Treat callable masks as coordinate-specific objects, not portable sample-independent BED files.
- Treat A/B core as the intersection of stable coordinate windows and callable sequence, not as either branch-specific callable BED alone.
- Prefer read-backed correction and conservative callable thresholds over aggressive contiguity gains.
- Use stable work directories and logs so each stage can be audited.

## Citation

A formal citation for this utility repository is being prepared. For now, cite the underlying tools used in the specific workflow you run.

## License

See `LICENSE`.
