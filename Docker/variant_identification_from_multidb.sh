#!/bin/bash
if [ $# -lt 3 ]; then
     echo "3 arguments are required.";
     exit 1;
fi

output=$1
group_name=$2
# サンプル名を配列で取得
args=("$@")
start_index=2
end_index=${#args[@]}
input_name_arr=("${args[@]:$start_index:$end_index}")

if [[ ! -e $output ]]; then
      echo ${output}" cannot be found.";
      exit 1;
fi
cd $output;

# variant call
for i in "${!input_name_arr[@]}"; do 
  if [[ ! -f 6_haplotypecaller/${group_name}_${i}.hg38.vcf.gz ]]; then
    mkdir 6_haplotypecaller/tmp;
    docker run \
        -u "$(id -u $USER):$(id -g $USER)" \
        -v /etc/passwd:/etc/passwd:ro \
        -v /etc/group:/etc/group:ro \
        --name nakamae_haplotypecaller --memory 120g -itv $PWD:/data -w /data --rm broadinstitute/gatk:4.3.0.0 \
        gatk HaplotypeCaller \
        -R 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
        -I 5_recal_data/${input_name_arr[$i]}.dbsnp_only.BQSR.bam \
        -O 6_haplotypecaller/${group_name}_${i}.hg38.vcf.gz \
        -ERC GVCF --tmp-dir 6_haplotypecaller/tmp \
        --sample-name ${input_name_arr[$i]};
    rm -rf 6_haplotypecaller/tmp;
  fi
  if [[ ! -f 6_haplotypecaller/${group_name}_${i}.hg38.vcf.gz ]]; then
    echo "Error."
    exit 1;
  fi
  echo -e "${group_name}_${i}\t6_haplotypecaller/${group_name}_${i}.hg38.vcf.gz" >> "7b_joint_genomics_dbImport/${group_name}.sample_map"
done

# make DB from sigle data
if [[ ! -f joint_db/gvcfs_joint_db_${group_name}/vidmap.json ]]; then
  docker run \
      -u "$(id -u $USER):$(id -g $USER)" \
      -v /etc/passwd:/etc/passwd:ro \
      -v /etc/group:/etc/group:ro \
      --name nakamae_genomicsdbimport --memory 120g -itv $PWD:/data -w /data --rm broadinstitute/gatk:4.3.0.0 \
      gatk GenomicsDBImport -R 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
      --sample-name-map 7b_joint_genomics_dbImport/${group_name}.sample_map \
      -L 7b_joint_genomics_dbImport/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list \
      --genomicsdb-workspace-path joint_db/gvcfs_joint_db_${group_name};
fi

if [[ ! -f joint_db/gvcfs_joint_db_${group_name}/vidmap.json ]]; then
    echo "Error."
    exit 1;
fi

# バリアント判定
if [[ ! -f 8_vcf_identification/${group_name}.hg38.identified.vcf; ]]; then
  docker run \
      -u "$(id -u $USER):$(id -g $USER)" \
      -v /etc/passwd:/etc/passwd:ro \
      -v /etc/group:/etc/group:ro \
      --name nakamae_genotypegvcfs --memory 120g -itv $PWD:/data -w /data --rm broadinstitute/gatk:4.3.0.0 \
      gatk GenotypeGVCFs -R 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
      -V gendb://joint_db/gvcfs_joint_db_${group_name} \
      -O 8_vcf_identification/${group_name}.hg38.identified.vcf;
fi

if [[ ! -f 8_vcf_identification/${group_name}.hg38.identified.vcf; ]]; then
    echo "Error."
    exit 1;
fi

echo "Done."
exit 0;