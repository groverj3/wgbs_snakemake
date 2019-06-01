# WGBS Snakemake Workflow
This workflow is designed to run the basic steps for a whole-genome bisulfite sequencing experiment. It's intended to automate the workflow for future-use and reproducibility. If you're looking for something which will do literally everything for an experiment this isn't the workflow you're looking for (yet). However, if you're comfortable running the individual programs and want to save yourself the trouble of running every step separately this will save you trouble.

The advantage of running this through Snakemake is that it intelligently handles threading and replaces completed processes up to the number of cores specified at run-time. Individual options for the steps are mostly hard-coded, as this is intended for reproducibility of our particular workflow. However, options for the thread count for each step are changeable from the .yaml file.

## Getting Started
Edit the .yaml file to include your sample IDs (excluding extensions, pair numbers, lane info, etc.) and an already indexed reference genome. Put your sample .fastq files in the "input_data subdirectory and run with:

`snakemake --cores {cores_here}`

## Dependencies
1. [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
2. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
3. [bwa-meth](https://github.com/brentp/bwa-meth)
4. [samtools](https://www.htslib.org/)
5. [Picard Tools](https://broadinstitute.github.io/picard/)
6. [MethylDackel](https://github.com/dpryan79/MethylDackel)

## Workflow
1. Quality checking, and output of sample information with FastQC
2. Adapter and quality trimming with Trim Galore!
3. Alignment to a reference genome with bwa-meth
4. Marking PCR duplicates with Picard Tools MarkDuplicates
5. Detecting methylation bias per read position with MethylDackel mbias
6. Extracting methylation calls per position into bedGraph and methylKit formats

More to come...
