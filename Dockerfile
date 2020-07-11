FROM ubuntu:16.04

# Update necessary packages
RUN apt-get update && apt-get install -qqy \
    awscli \
    build-essential \
    git \
    libbz2-dev \
    liblzma-dev \
    make \
    pkg-config \
    python \
    python3 \
    python3-dev \
    python3-pip \
    wget \
    unzip \
    zlib1g-dev

RUN pip3 install --upgrade pip

RUN pip3 install pybedtools pyvcf scipy numpy

# Install samtools (needed to index reference fasta files)
RUN wget -O samtools-1.9.tar.bz2 https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
RUN tar -xjf samtools-1.9.tar.bz2
WORKDIR samtools-1.9
RUN ./configure --without-curses && make && make install
WORKDIR ..

# Install bedtools (needed for DumpSTR)
RUN wget -O bedtools-2.27.1.tar.gz https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz
RUN tar -xzvf bedtools-2.27.1.tar.gz
WORKDIR bedtools2
RUN make && make install
WORKDIR ..

# Download, compile, and install GangSTR
RUN wget -O GangSTR-2.4.tar.gz https://github.com/gymreklab/GangSTR/releases/download/v2.4/GangSTR-2.4.tar.gz
RUN tar -xzvf GangSTR-2.4.tar.gz
WORKDIR GangSTR-2.4
RUN ./install-gangstr.sh
RUN ldconfig
WORKDIR ..

# Download and install TRTools
RUN git clone https://github.com/ileenamitra/TRTools
WORKDIR TRTools
RUN git checkout xchrom
RUN python3 setup.py install
WORKDIR ..
