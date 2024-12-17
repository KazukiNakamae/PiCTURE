############################################################
# Dockerfile to build STAR Mapping in human (GATK resource bundle hg38)
############################################################

FROM mambaorg/micromamba:0.13.1

# File Author / Maintainer
MAINTAINER Kazuki Nakamae
RUN conda install -c bioconda star=2.7.4a --yes \
  && conda clean --all --yes \
  && mkdir hg38_index

# Copy reference genome
COPY resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta .

RUN STAR --runMode genomeGenerate --genomeDir hg38_index --genomeFastaFiles resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta --limitGenomeGenerateRAM 200000000000