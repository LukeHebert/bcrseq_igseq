## More important

- For the bcrseq-igseq data merging script, do not filter out peptides that merely map to more than 1 bcrseq read. Instead, filter out heavy chain cdr3-covering peptides that map to >1 lineage and, to all peptides, assign them a "map_confidence" score (could use a better name) that shows how confident we are that it belongs to a given bcrseq transcript 

- For instances when user has separate (non natively paired) heavy and light chain sequencing data, provide a way to salvage off-target-but-still-relevant BCRseq data before the clustering step. A simple script to move around any light chain data captured that was accidentally captured with the heavy chain and vice versa might allow for e.g. slightly better mapping of heavy chain CDR3s to peptides.

- Adding to the above to-do, just make the filter_collapse.py script do the separation, making two output files (one for heavy chain data, one for light). If there is some heavy chain data, it automatically clusters it.

- Remove all absolute pathways and replace with user arguments to the corresponding software

- Work on the final BCRseq-IgSeq merging script to squeeze all valid mappings from the input datasets (e.g. a peptide "information" score?)

- Add clustering alternatives. At the very least, one should allow for (rare) indels to be clustered together

- Add a VH-VL pair-picking script to help mAb selection & downstream expression

- Add EMPEM data integration scripts?

- Add a filter for chimeric PCR products early on in the workflow using e.g. "CHMMAIRRa" (PMID: 40060431 doi: 10.1101/2025.02.21.638809)

- Add a flag for anytime a clone/BCR-seq has an unusually long CDR3 (e.g. >20 AA long) because these are often high binding and broadly neutralizing (e.g. in cows, in norovirus human IgA mAbs)

## Less Important

- Add a caveat in the README.md that says these scripts have been tested primarily on ferret input (for humans e.g. you will need to add C gene database to your IMGT reference datadata files)

- Add a little route in the diagram PNG that shows light chain data needs to skip the clustering step

- Integrate all scripts into a cohesive, user-friendly, bonafide pipeline with a tool like e.g. snakemake

- Include all dependencies in a container or similar

- Create (or at least update old) generally useful plotting scripts

- Make the workflow diagram darker

- Update README.md to guide users through a "pipeline" run.