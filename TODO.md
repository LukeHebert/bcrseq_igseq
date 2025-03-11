## More important

- For the bcrseq-igseq data merging script, do not filter out peptides that merely map to more than 1 bcrseq read. Instead, filter out heavy chain cdr3-covering peptides that map to >1 lineage and, to all peptides, assign them a "map_confidence" score (could use a better name) that shows how confident we are that it belongs to a given bcrseq transcript 

- For instances when user has separate (non natively paired) heavy and light chain sequencing data, provide a way to salvage off-target-but-still-relevant BCRseq data before the clustering step. A simple script to move around any light chain data captured that was accidentally captured with the heavy chain and vice versa might allow for e.g. slightly better mapping of heavy chain CDR3s to peptides.

- Remove all absolute pathways and replace with user arguments to the corresponding software

- Work on the final BCRseq-IgSeq merging script to squeeze all valid mappings from the input datasets (e.g. a peptide "information" score?)

- Add clustering alternatives. At the very least, one should allow for (rare) indels to be clustered together

- Add a VH-VL pair-picking script to help mAb selection & downstream expression

- Add EMPEM data integration scripts?

## Less Important

- Integrate all scripts into a cohesive, user-friendly, bonafide pipeline with a tool like e.g. snakemake

- Include all dependencies in a container or similar

- Create (or at least update old) generally useful plotting scripts

- Make the workflow diagram darker

- Update README.md to guide users through a "pipeline" run.