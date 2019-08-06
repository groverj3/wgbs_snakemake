# Author: Jeffrey Grover
# Purpose: Run the whole-genome bisulfite sequencing workflow
# Created: 2019-05-22


# Get overall workflow parameters from config.yaml
configfile: 'config.yaml'

SAMPLES = config['samples']
REFERENCE_GENOME = config['reference_genome']

rule all:
    input:
        expand('6_mosdepth/{sample}.sorted.markdupes.coverage.txt',
               sample=SAMPLES),
        expand('5_methyldackel_extract/{sample}.sorted.markdupes_{context}.{ext}',
               sample=SAMPLES, context=['CpG', 'CHG', 'CHH'],
               ext=['bedGraph', 'methylKit']),
        expand('2_trim_galore/{sample}_R{mate}_val_{mate}_fastqc.{ext}',
               sample=SAMPLES, mate=[1, 2], ext=['html', 'zip'])


# Index the reference genome
# ancient() will assume the reference is older than output files if they exist
rule bwameth_index:
    input:
        ancient(REFERENCE_GENOME)
    output:
        REFERENCE_GENOME + '.bwameth.c2t',
        REFERENCE_GENOME + '.bwameth.c2t.amb',
        REFERENCE_GENOME + '.bwameth.c2t.ann',
        REFERENCE_GENOME + '.bwameth.c2t.bwt',
        REFERENCE_GENOME + '.bwameth.c2t.pac',
        REFERENCE_GENOME + '.bwameth.c2t.sa'
    params:
        bwameth_path = config['paths']['bwameth_path'],
    shell:
        '{params.bwameth_path} index {input}'


# Index the reference genome with faidx
# ancient() will assume the reference is older than output files if they exist
rule samtools_faidx:
    input:
        ancient(REFERENCE_GENOME)
    output:
        REFERENCE_GENOME + '.fai'
    params:
        samtools_path = config['paths']['samtools_path']
    shell:
        '{params.samtools_path} faidx {input}'


# Run fastqc on the raw .fastq.gz files
rule fastqc_raw:
    input:
        'input_data/{sample}_R{mate}.fastq.gz'
    output:
        '1_fastqc_raw/{sample}_R{mate}_fastqc.html',
        '1_fastqc_raw/{sample}_R{mate}_fastqc.zip'
    params:
        fastqc_path = config['paths']['fastqc_path'],
        out_dir = '1_fastqc_raw/'
    shell:
        '{params.fastqc_path} -o {params.out_dir} {input}'


# Trim the read pairs
rule trim_galore:
    input:
        '1_fastqc_raw/{sample}_R1_fastqc.html',
        '1_fastqc_raw/{sample}_R1_fastqc.zip',
        '1_fastqc_raw/{sample}_R2_fastqc.html',
        '1_fastqc_raw/{sample}_R2_fastqc.zip',
        R1 = 'input_data/{sample}_R1.fastq.gz',
        R2 = 'input_data/{sample}_R2.fastq.gz'
    output:
        '2_trim_galore/{sample}_R1_val_1.fq.gz',
        '2_trim_galore/{sample}_R1.fastq.gz_trimming_report.txt',
        '2_trim_galore/{sample}_R2_val_2.fq.gz',
        '2_trim_galore/{sample}_R2.fastq.gz_trimming_report.txt'
    params:
        adapter_seq = config['trim_galore']['adapter_seq'],
        out_dir = '2_trim_galore',
        trim_galore_path = config['paths']['trim_galore_path']
    shell:
        '''
        {params.trim_galore_path} \
        --a {params.adapter_seq} \
        --gzip \
        --trim-n \
        --quality 20 \
        --output_dir {params.out_dir} \
        --paired \
        {input.R1} {input.R2} \
        '''


# Run fastqc on the trimmmed reads
rule fastqc_trimmmed:
    input:
        '2_trim_galore/{sample}_R{mate}.fastq.gz_trimming_report.txt',
        fq_gz = '2_trim_galore/{sample}_R{mate}_val_{mate}.fq.gz'
    output:
        '2_trim_galore/{sample}_R{mate}_val_{mate}_fastqc.html',
        '2_trim_galore/{sample}_R{mate}_val_{mate}_fastqc.zip'
    params:
        fastqc_path = config['paths']['fastqc_path'],
        out_dir = '2_trim_galore/'
    shell:
        '{params.fastqc_path} -o {params.out_dir} {input.fq_gz}'


# Align to the reference
rule bwameth_align:
    input:
        {rules.bwameth_index.output},
        R1 = '2_trim_galore/{sample}_R1_val_1.fq.gz',
        R2 = '2_trim_galore/{sample}_R2_val_2.fq.gz'
    output:
        temp('temp_data/{sample}.bam')
    log:
        '3_aligned_sorted_markdupes/{sample}.bwameth.log'
    threads:
        config['bwameth']['threads']
    params:
        bwameth_path = config['paths']['bwameth_path'],
        genome = REFERENCE_GENOME
    shell:
        '''
        {params.bwameth_path} \
        -t {threads} \
        --reference {params.genome} \
        {input.R1} {input.R2} \
        2> {log} \
        | samtools view -b - \
        > {output}
        '''


# Sort the output files
rule samtools_sort:
    input:
        'temp_data/{sample}.bam'
    output:
        temp('temp_data/{sample}.sorted.bam')
    threads:
        config['samtools_sort']['threads']
    params:
        samtools_path = config['paths']['samtools_path'],
        mem = config['samtools_sort']['mem']
    shell:
        '''
        {params.samtools_path} sort \
        -@ {threads} \
        -m {params.mem} \
        -O BAM \
        -T {input}.samtools_sort.tmp \
        -o {output} \
        {input} \
        '''


# Mark potential PCR duplicates with Picard Tools
rule mark_dupes:
    input:
        'temp_data/{sample}.sorted.bam'
    output:
        '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bam'
    log:
        '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.log'
    params:
        picard_path = config['paths']['picard_path']
    shell:
        '''
        {params.picard_path} MarkDuplicates \
        I={input} \
        O={output} \
        M={log}
        '''


rule samtools_index:
    input:
        '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bam'
    output:
        '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bai'
    threads:
        config['samtools_index']['threads']
    params:
        samtools_path = config['paths']['samtools_path']
    shell:
        '''
        {params.samtools_path} index \
        -@ {threads} \
        -b \
        {input} \
        {output}
        '''


# Run MethylDackel to get the inclusion bounds for methylation calling
rule methyldackel_mbias:
    input:
        {rules.samtools_faidx.output},
        '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bai',
        bam = '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bam'
    output:
        '4_methyldackel_mbias/{sample}.sorted.markdupes_OB.svg',
        '4_methyldackel_mbias/{sample}.sorted.markdupes_OT.svg',
        mbias = '4_methyldackel_mbias/{sample}.sorted.markdupes.mbias'
    threads:
        config['methyldackel']['threads']
    params:
        methyldackel_path = config['paths']['methyldackel_path'],
        out_prefix = '4_methyldackel_mbias/{sample}.sorted.markdupes',
        genome = REFERENCE_GENOME
    shell:
        '''
        {params.methyldackel_path} mbias \
        --CHG \
        --CHH \
        -@ {threads} \
        {params.genome} \
        ./{input.bam} \
        ./{params.out_prefix} \
        2> ./{output.mbias}
        '''


# Run MethylDackel to extract cytosine stats
rule methyldackel_extract:
    input:
        {rules.samtools_faidx.output},
        '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bai',
        bam = '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bam',
        mbias = '4_methyldackel_mbias/{sample}.sorted.markdupes.mbias'
    output:
        '5_methyldackel_extract/{sample}.sorted.markdupes_CpG.bedGraph',
        '5_methyldackel_extract/{sample}.sorted.markdupes_CHG.bedGraph',
        '5_methyldackel_extract/{sample}.sorted.markdupes_CHH.bedGraph',
        '5_methyldackel_extract/{sample}.sorted.markdupes_CpG.methylKit',
        '5_methyldackel_extract/{sample}.sorted.markdupes_CHG.methylKit',
        '5_methyldackel_extract/{sample}.sorted.markdupes_CHH.methylKit'
    threads:
        config['methyldackel']['threads']
    params:
        methyldackel_path = config['paths']['methyldackel_path'],
        out_prefix = './5_methyldackel_extract/{sample}.sorted.markdupes',
        genome = REFERENCE_GENOME
    shell:
        '''
        # Get bounds for inclusion from the .mbias files

        OT=$(cut -d ' ' -f 5 {input.mbias})
        OB=$(cut -d ' ' -f 7 {input.mbias})

        # Get a MethylKit compatible file

        {params.methyldackel_path} extract \
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

        {params.methyldackel_path} extract \
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
rule mosdepth:
    input:
        bam = '3_aligned_sorted_markdupes/{sample}.sorted.markdupes.bam'
    output:
        '6_mosdepth/{sample}.sorted.markdupes.mosdepth.global.dist.txt',
        '6_mosdepth/{sample}.sorted.markdupes.mosdepth.summary.txt',
        '6_mosdepth/{sample}.sorted.markdupes.per-base.bed.gz',
        '6_mosdepth/{sample}.sorted.markdupes.per-base.bed.gz.csi'
    threads:
        config['mosdepth']['threads']
    params:
        mapping_quality = config['mosdepth']['mapping_quality'],
        mosdepth_path = config['paths']['mosdepth_path'],
        out_prefix = '6_mosdepth/{sample}.sorted.markdupes'
    shell:
        '''
        {params.mosdepth_path} \
        -x \
        -t {threads} \
        -Q {params.mapping_quality} \
        {params.out_prefix} \
        {input.bam}
        '''


# Calculate the coverage from the mosdepth output
rule calc_coverage:
    input:
        bed = '6_mosdepth/{sample}.sorted.markdupes.per-base.bed.gz'
    output:
        '6_mosdepth/{sample}.sorted.markdupes.coverage.txt'
    params:
        genome = REFERENCE_GENOME
    shell:
        '''
        scripts/mosdepth_to_x_coverage.py \
        -f {params.genome} \
        -m {input.bed} \
        > {output}
        '''
