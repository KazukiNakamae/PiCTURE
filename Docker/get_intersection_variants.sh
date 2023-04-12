#!/bin/bash
if [ ${#args[@]} -lt 4 ]; then
     echo "3 arguments are required.";
     exit 1;
fi

output=$1
set_name=$2
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

# get intersection varitnts set

for element in "${input_name_arr[@]}"; do
  file_name_arr+=("8_vcf_identification/${element}.hg38.identified.vcf")
done

if [[ ! -f 6_haplotypecaller/${group_name}_${i}.hg38.vcf.gz ]]; then
    echo "Add VAF..."
    docker run \
        -u "$(id -u $USER):$(id -g $USER)" \
        -v /etc/passwd:/etc/passwd:ro \
        -v /etc/group:/etc/group:ro \
        --name nakamae_vaf_calculation --memory 120g -itv $PWD:/data -w /data --rm staphb/bcftools:1.16 \
        bcftools isec -n=${#file_name_arr[@]} file_name_arr[@] -o 8_vcf_identification/${set_name}.hg38.identified.vcf;
fi
if [[ ! -f 8_vcf_identification/${set_name}.hg38.identified.vcf ]]; then
    echo "Error."
    exit 1;
fi

echo "Done."
exit 0;