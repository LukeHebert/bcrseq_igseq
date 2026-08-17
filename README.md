# BCR Transcript Sequencing and Immunoglobulin Protein Sequencing Pipeline

This repository contains two supported workflows:

- `BCR-seq Transcript Analysis` for processing sequencing reads into clustered BCR annotations and searchable FASTA databases.
- `Ig-seq Bottom-up Proteomic Analysis` for filtering Proteome Discoverer PSM exports, quantifying mapped peptides, and visualizing lineage-level protein evidence.

## Python environment

This project uses a standard Python virtual environment plus `requirements.txt`.

Create and activate a local environment from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Then run scripts with the activated environment's `python`.

## Official workflow files

| Workflow | Official entry files |
| --- | --- |
| BCR-seq Transcript Analysis | `workflows/bcrseq_transcript/bcrseq_pipeline.py`, `workflows/bcrseq_transcript/pipeline_config_template.json` |
| Ig-seq Bottom-up Proteomic Analysis | `workflows/igseq_proteomics/filter_psms.py`, `workflows/igseq_proteomics/quantify_map_peptides.py`, `workflows/igseq_proteomics/plot_lineage_repertoire.py`, `workflows/igseq_proteomics/compare_cdr_lineage_abundance.py` |
| Secondary analysis utilities | `analysis_utils/` |

Repository layout:

- `workflows/bcrseq_transcript/` contains the transcript-side pipeline entrypoints and stage scripts.
- `workflows/igseq_proteomics/` contains the proteomics-side filtering, mapping, and plotting scripts.
- `examples/` contains example configs, example commands, and fixture inputs for onboarding.
- `analysis_utils/` contains secondary helper scripts that are not the main onboarding entrypoints.

## BCR-seq Transcript Analysis

```mermaid
flowchart TD
    A[Sample FASTQ files] --> B[workflows/bcrseq_transcript/trim_merge.py]
    B --> C[Trimmed and merged reads]
    C --> D[workflows/bcrseq_transcript/identify_genes.py]
    D --> E[IgBLAST-annotated BCRseq TSV]
    E --> F[workflows/bcrseq_transcript/filter_collapse.py]
    F --> G[Filtered collapsed BCRseq records]
    G --> H[workflows/bcrseq_transcript/gupta_cluster.py]
    H --> I[Clustered BCRseq annotation TSV]
    I --> J[workflows/bcrseq_transcript/make_searchable.py]
    J --> K[Complementary searchable FASTA database]
```

Start here if you have paired-end transcript sequencing reads and want clustered BCR annotations plus a searchable FASTA database.

Transcript quickstart:

```bash
python workflows/bcrseq_transcript/bcrseq_pipeline.py \
  examples/fixtures/transcript \
  examples/configs/transcript_dry_run_config.json \
  --dry_run
```

This bundled transcript example is a dry-run wiring check. It demonstrates sample discovery, config resolution, and stage command construction.

## Ig-seq Bottom-up Proteomic Analysis

```mermaid
flowchart TD
    A[Proteome Discoverer PSM export] --> B[workflows/igseq_proteomics/filter_psms.py]
    B --> C[Heavy-chain single-lineage PSMs]
    C --> D[workflows/igseq_proteomics/quantify_map_peptides.py]
    E[BCRseq annotation TSV] --> D
    D --> F[Mapped quantified peptides]
    F --> G[workflows/igseq_proteomics/plot_lineage_repertoire.py]
    F --> J[workflows/igseq_proteomics/compare_cdr_lineage_abundance.py]
    E --> G
    G --> H[Lineage abundance plots and TSVs]
    G --> I[Per-lineage logo and coverage plots with TSVs]
    J --> K[CDR-derived lineage-abundance comparison TSVs and plots]
```

Start here if you already have a Proteome Discoverer PSM export and a clustered BCRseq annotation TSV.

Proteomics quickstart:

```bash
python workflows/igseq_proteomics/filter_psms.py \
  examples/fixtures/proteomics/demo_psms.tsv \
  --out-dir scratch/proteomics_example/01_filtered

python workflows/igseq_proteomics/quantify_map_peptides.py \
  scratch/proteomics_example/01_filtered/demo_psms_filtered_heavy_single_lineage.tsv \
  examples/fixtures/proteomics/demo_bcrseq.tsv \
  examples/fixtures/proteomics/demo_suffix.txt \
  --out-dir scratch/proteomics_example/02_mapped

python workflows/igseq_proteomics/plot_lineage_repertoire.py \
  scratch/proteomics_example/02_mapped/demo_psms_filtered_heavy_single_lineage_mapped_peptides.tsv \
  examples/fixtures/proteomics/demo_bcrseq.tsv \
  --out-dir scratch/proteomics_example/03_plots

python workflows/igseq_proteomics/compare_cdr_lineage_abundance.py \
  scratch/proteomics_example/02_mapped/demo_psms_filtered_heavy_single_lineage_mapped_peptides.tsv \
  --out-dir scratch/proteomics_example/04_cdr_comparison
```

### CDR-derived lineage-abundance comparison

`compare_cdr_lineage_abundance.py` compares five lineage-abundance selections from a mapped peptide TSV. A peptide must overlap a CDR by at least three amino acids by default; use `--min-cdr-overlap-aa` to change this threshold and `--top-n` to change the default top-10 comparison.

- CDR1-only: overlaps CDR1 but not CDR2.
- CDR2-only: overlaps CDR2 but not CDR1.
- Combined CDR1/CDR2: the union of CDR1-only, CDR2-only, and peptides overlapping both regions.
- CDR3: overlaps CDR3.
- All peptides: every successfully mapped peptide.

The script writes one five-selection lineage-abundance table, a seven-comparison summary, one lineage and top-N membership table per comparison, and a log-log scatterplot for comparisons with at least two lineages detected by both methods. It uses total precursor abundance and counts each `(ClusterID, peptide_sequence)` only once.

## Examples and fixtures

- [examples/README.md](examples/README.md) describes the bundled example inputs and what each one is intended to validate.
- `examples/fixtures/transcript/` contains a transcript dry-run fixture for sample naming and config validation.
- `examples/fixtures/proteomics/` contains a runnable proteomics example.
