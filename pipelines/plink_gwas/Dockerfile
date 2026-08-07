FROM mambaorg/micromamba:1.5.10

USER root
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      awscli=1.* htslib=1.20 plink=1.90b6.21 plink2=2.00a5.10 \
      python=3.11 matplotlib \
    && micromamba clean --all --yes
USER $MAMBA_USER
