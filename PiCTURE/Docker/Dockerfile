############################################################
# Dockerfile to build STAR Mapping in human (GATK resource bundle hg38)
############################################################

FROM mambaorg/micromamba:0.13.1

RUN conda install -c bioconda star=2.7.4a --yes \
  && conda clean --all --yes \
  && mkdir hg38_index

# Copy reference genome
COPY Homo_sapiens_assembly38.fasta resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta

RUN STAR --runThreadN 16 --genomeSAsparseD 3 --runMode genomeGenerate --genomeDir hg38_index --genomeFastaFiles resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta --limitGenomeGenerateRAM 100000000000