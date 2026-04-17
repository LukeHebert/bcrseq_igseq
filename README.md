# BCR Transcript Sequencing and Immunoglobulin Protein Sequencing Pipeline

## Python environment

This project now uses a standard Python virtual environment plus
`requirements.txt` rather than a conda environment export.

Create and activate a local environment from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Then run scripts with the activated environment's `python`.

## BCR-seq Transcript Analysis

```mermaid
flowchart LR
    A[Sample FASTQ files] --> B[trim_merge.py]
    B --> C[Trimmed and merged reads]
    C --> D[identify_genes.py]
    D --> E[IgBLAST-annotated BCRseq TSV]
    E --> F[filter_collapse.py]
    F --> G[Filtered collapsed BCRseq records]
    G --> H[cluster.py]
    H --> I[Clustered BCRseq annotation TSV]
    I --> J[make_searchable.py]
    J --> K[Complementary searchable FASTA database]
```

## Ig-seq Bottom-up Proteomic Analysis

```mermaid
flowchart LR
    A[Proteome Discoverer PSM export] --> B[filter_psms.py]
    B --> C[Heavy-chain single-lineage PSMs]
    C --> D[quantify_map_peptides.py]
    E[BCRseq annotation TSV] --> D
    D --> F[Mapped quantified peptides]
    F --> G[plot_lineage_repertoire.py]
    E --> G
    G --> H[Lineage abundance plots and TSVs]
    G --> I[Per-lineage logo and coverage plots with TSVs]
```
