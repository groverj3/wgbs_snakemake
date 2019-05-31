# Author: Jeffrey Grover
# Purpose: Run the whole-genome bisulfite sequencing workflow
# Created: 2019-05-22

configfile: 'config.yaml'


# Get overall workflow parameters from config.yaml

samples = config['samples']
reference_genome = config['reference_genome']

# rule all:
#    input:


# Combine the individual files from each sample's R1 and R2 files
rule concatenate_R1:
    input:
        '{samples}'
    output:
        temp('temp_data/{samples}_1.fastq')
    shell:
        'zcat input_data/{input}/{input}*R1*.fastq.gz > {output}'

rule concatenate_R2:
    input:
        '{samples}'
    output:
        temp('temp_data/{samples}_2.fastq')
    shell:
        'zcat input_data/{input}/{input}*R2*.fastq.gz > {output}'


# Run fastqc and keep the output
rule fastqc_cat:
    input:
        'temp_data/{samples}_1.fastq',
        'temp_data/{samples}_2.fastq'
    output:
        '1_fastqc/{samples}_1_fastqc.html',
        '1_fastqc/{samples}_2_fastqc.html',
        '1_fastqc/{samples}_1_fastqc.zip',
        '1_fastqc/{samples}_2_fastqc.zip'
    params:
        out_dir = '1_fastqc/'
    shell:
        'fastqc -o {params.out_dir} {input}'


# Trim the concatenated files
rule trim_galore:
    input:
        'temp_data/{samples}_1.fastq.gz',
        'temp_data/{samples}_2.fastq.gz'
    output:
        '2_trim_galore/{samples}_1_val_1.fq.gz',
        '2_trim_galore/{samples}_2_val_2.fq.gz',
        '2_trim_galore/{samples}_1_val_1_fastqc.html',
        '2_trim_galore/{samples}_1_val_2_fastqc.html'
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


# Align to the reference, I'll handle conversion rate calculation some other time
rule bwameth_reference:
    input:
        '2_trim_galore/{samples}_1_val_1.fq.gz',
        '2_trim_galore/{samples}_2_val_2.fq.gz',
    output:
        temp('temp_data/{samples}.bam')
    threads:
        config[bwameth]['threads']
    shell:
        '''
        bwameth.py \
        -t {threads} \
        --reference {reference_genome} \
        {input} \
        | samtools view -bhS - > {output}
        '''


# Sort the output files
rule samtools_sort:
    input:
        'temp_data/{samples}.bam'
    output:
        temp('temp_Data/{samples}.sorted.bam')
    threads:
        config[samtools]['threads']
    shell:
        'samtools sort -@ {threads} -O BAM {input} > {output}'


# Mark potential PCR duplicates with Picard Tools
rule mark_dupes:
    input:
        'temp_data/{samples}.sorted.bam'
    output:
        '3_aligned_sorted_markdupes/{samples}.sorted.markdupes.bam'
    params:
        picard_path = config[mark_duplicates]['picard_path']
    shell:
        '''
        java -jar \
        {params.picard_path} \
        I={input} \
        O={output} \
        M={output}.log
        '''


# Index the sorted and duplicate-marked bam file
rule index_sorted_marked_bam:
    input:
        '3_aligned_sorted_markdupes/{samples}.sorted.markdupes.bam'
    output:
        '3_aligned_sorted_markdupes/{samples}.sorted.markdupes.bai'
    threads:
        config[samtools_index]['threads']
    shell:
        'samtools index -@ {threads} {input} {output}'


# Run MethylDackel to get the inclusion bounds for methylation calling
rule methyldackel_mbias:
    input:
        bam = '3_aligned_sorted_markdupes/{samples}.sorted.markdupes.bam'
        index = '3_aligned_sorted_markdupes/{samples}.sorted.markdupes.bai'
    output:
        mbias = '4_methyldackel/{samples}.sorted.markdupes.mbias'
        ob_plot = '4_methyldackel/{samples}.sorted.markdupes.bam_OB.svg'
        ot_plot = '4_methyldackel/{samples}.sorted.markdupes.bam_OT.svg'
    threads:
        config[methyldackel]['threads']
    params:
        out_prefix = '4_methyldackel/{samples}.sorted.markdupes'
    shell:
        '''
        MethylDackel mbias \
        --CHG \
        --CHH \
        -@ {threads} \
        {genome} \
        {input.bam} \
        {params.out_prefix} \
        2> {output.mbias}
        '''


# Run MethylDackel to extract cytosine stats
rule methyldackel_extract:
    input:
        bam = '3_aligned_sorted_markdupes/{samples}.sorted.markdupes.bam'
        index = '3_aligned_sorted_markdupes/{samples}.sorted.markdupes.bai'
        mbias = '4_methyldackel/{samples}.sorted.markdupes.mbias'
    output:
        '4_methyldackel/{samples}.sorted.markdupes_CpG.bedGraph',
        '4_methyldackel/{samples}.sorted.markdupes_CHG.bedGraph',
        '4_methyldackel/{samples}.sorted.markdupes_CHH.bedGraph',
        '4_methyldackel/{samples}.sorted.markdupes_CpG.methylKit',
        '4_methyldackel/{samples}.sorted.markdupes_CHG.methylKit',
        '4_methyldackel/{samples}.sorted.markdupes_CHH.methylKit'
    threads:
        config[methyldackel]['threads']
    params:
        out_prefix = '4_methyldackel/{samples}.sorted.markdupes'
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
        {genome} \
        {input.bam}

        # Get the normal bedGraph output file

        MethylDackel extract \
        --CHG \
        --CHH \
        --OT $OT \
        --OB $OB \
        -@ {threads} \
        -o {params.out_prefix} \
        {genome} \
        {input.bam}
        '''


# TODO: the rest of the pipeline below this

