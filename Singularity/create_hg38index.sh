#!/bin/bash
output=$1

echo "Indexing hg38 genome using STAR"
if [ ! -f ./hg38_index/genomeParameters.txt ]; then
     echo "Create STAR index in docker image"
     singularity exec star.sif STAR --runMode genomeGenerate --genomeDir $1/hg38_index --genomeFastaFiles $1/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta --limitGenomeGenerateRAM 200000000000
else
     echo "The docker image already existed."
fi
