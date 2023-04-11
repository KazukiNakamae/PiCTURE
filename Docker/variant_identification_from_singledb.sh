#!/bin/bash
input_name=$1
output=$2

if [[ ! -e $output ]]; then
      echo ${output}" cannot be found.";
      exit 1;
fi
cd $output;

# variant call
if [[ ! -f 6_haplotypecaller/${input_name}.hg38.vcf.gz ]]; then
  mkdir 6_haplotypecaller/tmp;
  docker run \
      -u "$(id -u $USER):$(id -g $USER)" \
      -v /etc/passwd:/etc/passwd:ro \
      -v /etc/group:/etc/group:ro \
      --name nakamae_haplotypecaller --memory 120g -itv $PWD:/data -w /data --rm broadinstitute/gatk:4.3.0.0 \
      gatk HaplotypeCaller \
      -R 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
      -I 5_recal_data/${input_name}.dbsnp_only.BQSR.bam \
      -O 6_haplotypecaller/${input_name}.hg38.vcf.gz \
      -ERC GVCF --tmp-dir 6_haplotypecaller/tmp \
      --sample-name ${input_name};
  rm -rf 6_haplotypecaller/tmp;
fi

if [[ ! -f 6_haplotypecaller/${input_name}.hg38.vcf.gz ]]; then
    echo "Error."
    exit 1;
fi

# make DB from sigle data
if [[ ! -f single_db/gvcfs_single_db_${input_name}/vidmap.json ]]; then
  docker run \
      -u "$(id -u $USER):$(id -g $USER)" \
      -v /etc/passwd:/etc/passwd:ro \
      -v /etc/group:/etc/group:ro \
      --name nakamae_genomicsdbimport --memory 120g -itv $PWD:/data -w /data --rm broadinstitute/gatk:4.3.0.0 \
      gatk GenomicsDBImport -R 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
      -V 6_haplotypecaller/${input_name}.hg38.vcf.gz \
      -L 7a_single_genomics_dbImport/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list \
      --genomicsdb-workspace-path single_db/gvcfs_single_db_${input_name};
fi

if [[ ! -f single_db/gvcfs_single_db_${input_name}/vidmap.json ]]; then
    echo "Error."
    exit 1;
fi

# バリアント判定
if [[ ! -f 8_vcf_identification/${input_name}.hg38.identified.vcf ]]; then
  docker run \
      -u "$(id -u $USER):$(id -g $USER)" \
      -v /etc/passwd:/etc/passwd:ro \
      -v /etc/group:/etc/group:ro \
      --name nakamae_genotypegvcfs --memory 120g -itv $PWD:/data -w /data --rm broadinstitute/gatk:4.3.0.0 \
      gatk GenotypeGVCFs -R 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
      -V gendb://single_db/gvcfs_single_db_${input_name} \
      -O 8_vcf_identification/${input_name}.hg38.identified.vcf;
fi

if [[ ! -f 8_vcf_identification/${input_name}.hg38.identified.vcf ]]; then
    echo "Error."
    exit 1;
fi

echo "Done."
exit 0;