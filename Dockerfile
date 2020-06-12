FROM continuumio/miniconda3:4.7.12

LABEL maintainer="Wen-Wei Liang <liang.w@wustl.edu>"

# Configure locale and timezone
RUN echo "America/Chicago" > /etc/timezone && \
    rm /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata

RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -y python=3.7 \
        snakemake-minimal=5.10.0 \
        pandas=1.0.1 \
        scipy \
    && conda clean -y --all
