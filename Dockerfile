# Filename: Dockerfile
FROM ubuntu:16.04
MAINTAINER Jeffrey Grover <jeffrey.w.grover@gmail.com>


# Install system software dependencies
RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y autoconf automake curl default-jre bsdtar make gcc git \
    tar wget python3 python3-dev python3-pip python python-dev python-pip \
    libtbb2 libtbb-dev libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-gnutls-dev libssl-dev


# Install FastQC 0.11.8
WORKDIR /opt
RUN wget -qO- https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip \
    | bsdtar -xvf- \
    && chmod +x /opt/FastQC/fastqc \
    && ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc


# Install cutadapt v2.9 with pip
RUN pip3 install 'cutadapt=2.9'


# Install trim_galore 0.6.4
RUN curl -L https://github.com/FelixKrueger/TrimGalore/archive/0.6.4.tar.gz \
    | tar xz \
    && ln -s /opt/TrimGalore-0.6.4/trim_galore /usr/local/bin/trim_galore


# Install Samtools 1.9
RUN curl -L https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    | tar xj
WORKDIR ./samtools-1.9
RUN ./configure
RUN make && make install
WORKDIR /opt


# Install bwa 0.7.17
RUN curl -L https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 \
    | tar xj
WORKDIR /opt/bwa-0.7.17
RUN make && ln -s /opt/bwa-0.7.17/bwa /usr/local/bin
WORKDIR /opt


# Install bwa-meth, which due to stupid ubuntu 16 defaults to python2
# Require Brent Pedersen's toolshed package (Python 2). No updates since 2016, 0.4.6 is the most recent.
RUN pip install 'toolshed=0.4.6'
RUN curl -L https://github.com/brentp/bwa-meth/archive/v0.2.2.tar.gz \
    | tar xz \
    && ln -s /opt/bwa-meth-0.2.2/bwameth.py /usr/local/bin


# Install htslib 1.9
RUN curl -L https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 \
    | tar xj
WORKDIR /opt/htslib-1.9
RUN ./configure && make && make install
WORKDIR /opt


# Install libBigWig 0.4.4
RUN curl -L https://github.com/dpryan79/libBigWig/archive/0.4.4.tar.gz \
    | tar xz
WORKDIR /opt/libBigWig-0.4.4
RUN make && make install
WORKDIR /opt


# Make sure htslib and libBigWig can be found
ENV LD_LIBRARY_PATH /usr/local/lib/:$LD_LIBRARY_PATH


# Install MethylDackel 0.4
RUN curl -L https://github.com/dpryan79/MethylDackel/archive/0.4.0.tar.gz \
    | tar xz
WORKDIR /opt/MethylDackel-0.4.0
RUN make && make install
WORKDIR /opt


# Install picard tools 2.22.3
RUN mkdir /opt/picard_tools-2.22.3
WORKDIR /opt/picard_tools-2.22.3
RUN wget https://github.com/broadinstitute/picard/releases/download/2.22.3/picard.jar \
    && printf '#!/bin/sh\nexec /usr/bin/java $JVM_OPTS -jar /opt/picard_tools-2.22.3/picard.jar "$@"' > /usr/local/bin/picard \
    && chmod +x /usr/local/bin/picard
WORKDIR /opt


# Install mosdepth 0.2.9
RUN mkdir /opt/mosdepth-0.2.9
WORKDIR /opt/mosdepth-0.2.9
RUN wget https://github.com/brentp/mosdepth/releases/download/v0.2.9/mosdepth \
    && chmod +x mosdepth \
    && ln -s /opt/mosdepth-0.2.9/mosdepth /usr/local/bin


# Install Snakemake v5.14.0
RUN pip3 install 'snakemake=5.14.0'


# Set final WORKDIR back to root
WORKDIR /
