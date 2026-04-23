# Examples and fixtures

This directory contains example configs, example commands, and fixture inputs for onboarding.

## Transcript example

- `configs/transcript_dry_run_config.json`
  Minimal config for exercising the transcript pipeline entrypoint in `--dry_run` mode.
- `fixtures/transcript/demo_sample_IgG_ttblood/`
  Tiny paired FASTQ files used to validate sample discovery, filename matching, and command construction.

Run:

```bash
python workflows/bcrseq_transcript/bcrseq_pipeline.py \
  examples/fixtures/transcript \
  examples/configs/transcript_dry_run_config.json \
  --dry_run
```

## Proteomics example

- `fixtures/proteomics/demo_psms.tsv`
  Tiny Proteome Discoverer-style PSM table.
- `fixtures/proteomics/demo_bcrseq.tsv`
  Tiny clustered BCRseq annotation table with one heavy-chain lineage.
- `fixtures/proteomics/demo_suffix.txt`
  Tiny suffix/extension file for peptide-to-BCR mapping.

Run:

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
```
