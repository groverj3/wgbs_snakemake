# Author: Jeffrey Grover
# Purpose: Run the whole-genome bisulfite sequencing workflow
# Created: 2019-05-22

configfile: 'config.yaml'


# Get overall workflow parameters from config.yaml

SAMPLES = config['samples']
REFERENCE_GENOME = config['reference_genome']

rule all:
    input:
        expand('4_methyldackel/{sample}.sorted.markdupes_{context}.{ext}',
               sample=SAMPLES,
               context=[CpG, CHG, CHH],
               ext=[bedGraph, methylKit]
               )


# Combine the individual files from each sample's R1 and R2 files
rule concatenate_reads:
    input:
        'input_data/{sample}/{samples}{lane}R{mate}{id}.fastq.gz'
    output:
        temp('temp_data/{sample}_R{mate}.fastq')
    wildcard_constraints:
        mate = '1|2'
    shell:
        'zcat {input} > {output}'


# Run fastqc and keep the output
rule fastqc_cat:
    input:
        'temp_data/{sample}_R{mate}.fastq'
    output:
        '1_fastqc/{sample}_R{mate}_fastqc.html',
        '1_fastqc/{sample}_R{mate}_fastqc.zip'
    wildcard_constraints:
        mate = '1|2'
    params:
        out_dir = '1_fastqc/'
    shell:
        'fastqc -o {params.out_dir} {input}'


# Trim the concatenated files
rule trim_galore:
    input:
        'temp_data/{sample}_R1.fastq.gz',
        'temp_data/{sample}_R2.fastq.gz'
    output:
        '2_trim_galore/{sample}_R{mate}_val_{mate}.fq.gz',
        '2_trim_galore/{sample}_R{mate}_val_{mate}_fastqc.html'
    wildcard_constraints:
        mate = '1|2'
    params:
        adapter_seq = config[trim_galore]['adapter_seq']
        out_dir = '2_trim_galore'
    shell:
        '''
        trim_galore \
        --a {params.adapter_seq} \
        --gzip \
        --fastqc \
        --trim-n \
        --quality 20 \
        --output_dir {params.out_dir} \
        --paired \
        {input}
        '''


# Align to the reference
rule bwameth_reference:
    input:
        '2_trim_galore/{sample}_R1_val_1.fq.gz',
        '2_trim_galore/{sample}_R2_val_2.fq.gz'
    output:
        temp('temp_data/{sample}.bam')
    threads:
        config[bwameth]['threads']
    params:
        genome = REFERENCE_GENOME
    shell:
        '''
        bwameth.py -t {threads} --reference {params.genome} {input} \
        | samtools view -bhS - > {output}
        '''


# Sort the output files
rule samtools_sort:
    input:
        'temp_data/{sample}.bam'
    output:
        temp('temp_Data/{sample}.sorted.bam')
    threads:
        config[samtools]['threads']
    shell:
        'samtools sort -@ {threads} -O BAM {input} > {output}'


# Mark potential PCR duplicates with Picard Tools
rule mark_dupes:
    input:
        'temp_data/{sample}.sorted.bam'
    output:
        '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bam'
    params:
        picard_path = config[mark_duplicates]['picard_path']
    shell:
        'java -jar {params.picard_path} I={input} O={output} M={output}.log'
        

# Index the sorted and duplicate-marked bam file
rule index_sorted_marked_bam:
    input:
        '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bam'
    output:
        '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bai'
    threads:
        config[samtools_index]['threads']
    shell:
        'samtools index -@ {threads} {input} {output}'


# Run MethylDackel to get the inclusion bounds for methylation calling
rule methyldackel_mbias:
    input:
        bam = '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bam',
        index = '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bai'
    output:
        mbias = '4_methyldackel/{sample}.sorted.markdupes.mbias',
        ob_plot = '4_methyldackel/{sample}.sorted.markdupes.bam_OB.svg',
        ot_plot = '4_methyldackel/{sample}.sorted.markdupes.bam_OT.svg'
    threads:
        config[methyldackel]['threads']
    params:
        out_prefix = '4_methyldackel/{sample}.sorted.markdupes',
        genome = REFERENCE_GENOME
    shell:
        '''
        MethylDackel mbias \
        --CHG \
        --CHH \
        -@ {threads} \
        {params.genome} \
        {input.bam} \
        {params.out_prefix} \
        2> {output.mbias}
        '''


# Run MethylDackel to extract cytosine stats
rule methyldackel_extract:
    input:
        bam = '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bam',
        index = '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bai',
        mbias = '4_methyldackel/{sample}.sorted.markdupes.mbias'
    output:
        '4_methyldackel/{sample}.sorted.markdupes_{context}.bedGraph',
        '4_methyldackel/{sample}.sorted.markdupes_{context}.methylKit'
    wildcard_constraints:
        context = 'CpG|CHG|CHH'
    threads:
        config[methyldackel]['threads']
    params:
        out_prefix = '4_methyldackel/{sample}.sorted.markdupes',
        genome = REFERENCE_GENOME
    shell:
        '''
        # Get bounds for inclusion

        OT=$(cut -d ' ' -f 5 {input.mbias})
        OB=$(cut -d ' ' -f 7 {input.mbias})

        # Get a MethylKit compatible file

        MethylDackel extract \
        --CHG \
        --CHH \
        --OT $OT \
        --OB $OB \
        --methylKit \
        -@ {threads} \
        -o {params.out_prefix} \
        {params.genome} \
        {input.bam}

        # Get the normal bedGraph output file

        MethylDackel extract \
        --CHG \
        --CHH \
        --OT $OT \
        --OB $OB \
        -@ {threads} \
        -o {params.out_prefix} \
        {params.genome} \
        {input.bam}
        '''


# Get the depth for each sample
rule get_depth:
    input:
        '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bam'
    output:
        '5_mosdepth/{sample}.sorted.markdupes.mosdepth.global.dist.txt',
        '5_mosdepth/{sample}.sorted.markdupes.per-base.bed.gz',
        '5_mosdepth/{sample}.sorted.markdupes.per-base.bed.gz.csi'
    threads:
        config[mosdepth]['threads']
    params:
        mapping_quality = config[mosdepth]['mapping_quality'],
        out_prefix = '5_mosdepth/{sample}.sorted.markdupes'
    shell:
        '''
        mosdepth \
        -x \
        -t {threads} \
        -Q {params.mapping_quality} \
        {params.out_prefix} \
        {input}
        '''


# Calculate the coverage from the mosdepth output
rule calc_coverage:
    input:
        '5_mosdepth/{sample}.sorted.markdupes.per-base.bed.gz'
    output:
        '5_mosdepth/{sample}.sorted.markdupes.coverage.txt'
    params:
        genome = REFERENCE_GENOME
    shell:
        '''
        scripts/mosdepth_per-base_to_x_coverage.py \
        -f {params.genome} \
        -m {input} \
        > {output}
        '''


# TODO: the rest of the pipeline below this

