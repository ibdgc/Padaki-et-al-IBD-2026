FROM ubuntu:20.04

ARG GITHUB_REPO
LABEL org.opencontainers.image.title="bge-cutler-qc"
LABEL org.opencontainers.image.description="Container for BGE genome Cutler QC workflows"
LABEL org.opencontainers.image.source="${GITHUB_REPO}"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install basic dependencies
RUN apt-get update && apt-get install -y \
    perl \
    wget \
    gzip \
    unzip \
    tar \
    curl \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install PLINK 1.9
RUN cd /tmp && \
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip && \
    unzip plink_linux_x86_64_20241022.zip && \
    mv plink /usr/local/bin/plink1.9 && \
    chmod +x /usr/local/bin/plink1.9 && \
    rm -f /tmp/plink_linux_x86_64_20241022.zip && \
    rm -f /tmp/LICENSE /tmp/README.txt

# Install PLINK 2
RUN cd /tmp && \
    wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20250129.zip && \
    unzip plink2_linux_avx2_20250129.zip && \
    mv plink2 /usr/local/bin/ && \
    chmod +x /usr/local/bin/plink2 && \
    rm -f /tmp/plink2_linux_avx2_20250129.zip

# Create working directory
WORKDIR /app

# Copy workflow scripts
COPY workflows/add_rs_to_pvar.pl /usr/local/bin/
COPY workflows/make_bge_exonic_list.pl /usr/local/bin/
COPY workflows/unpacking.sh /usr/local/bin/

# Make scripts executable
RUN chmod +x /usr/local/bin/add_rs_to_pvar.pl \
    /usr/local/bin/make_bge_exonic_list.pl \
    /usr/local/bin/unpacking.sh

# Set PATH
ENV PATH="/usr/local/bin:${PATH}"
