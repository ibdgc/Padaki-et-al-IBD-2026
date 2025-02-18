FROM ubuntu:20.04

ARG GITHUB_REPO
LABEL org.opencontainers.image.title="bge-analysis"
LABEL org.opencontainers.image.description="Container for BGE genome analysis workflows"
LABEL org.opencontainers.image.source="${GITHUB_REPO}"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install basic dependencies
RUN apt-get update && apt-get install -y \
    perl \
    wget \
    gzip \
    tar \
    curl \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install PLINK 2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20250129.zip && \
    unzip plink2_linux_avx2_20250129.zip && \
    mv plink2 /usr/local/bin/ && \
    rm plink2_linux_avx2_20250129.zip


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
