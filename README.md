# WGBS Snakemake Workflow
This workflow is designed to run the basic steps for a whole-genome bisulfite sequencing experiment. It's intended to automate the workflow for future-use and reproducibility.

## Dependencies
1. trim_galore
2. fastqc
3. bwa-meth
4. samtools
5. Picard Tools

## Workflow
1. Adapter and quality trimming with Trim Galore!
2. Alignment to a reference genome with bwa-meth
3. Marking PCR duplicates with Picard Tools MarkDuplicates

More to come...
