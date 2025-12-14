#!/bin/bash
if [ $# -lt 4 ]; then
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

file_name_arr=();
for element in "${input_name_arr[@]}"; do
  # file_name="--variant 8_vcf_identification/${element}.hg38.identified.vcf"
  file_name="8_vcf_identification/${element}.hg38.identified.vcf.gz"
  file_name_arr+=("$file_name")
done

for element in "${input_name_arr[@]}"; do
    if [[ ! -f 8_vcf_identification/${element}.hg38.identified.vcf.idx ]]; then
        echo "Indexing VCF..."
        docker run \
            -u "$(id -u $USER):$(id -g $USER)" \
            -v /etc/passwd:/etc/passwd:ro \
            -v /etc/group:/etc/group:ro \
            --name nakamae_haplotypecaller --memory 120g -itv $PWD:/data -w /data --rm broadinstitute/gatk:4.3.0.0 \
            gatk IndexFeatureFile \
            -I 8_vcf_identification/${element}.hg38.identified.vcf;
    fi
done

# if [[ ! -f 8_vcf_identification/${set_name}.merge.vcf ]]; then
#     echo "Merge VCF..."
#     docker run \
#         -u "$(id -u $USER):$(id -g $USER)" \
#         -v /etc/passwd:/etc/passwd:ro \
#         -v /etc/group:/etc/group:ro \
#         --name nakamae_haplotypecaller --memory 120g -itv $PWD:/data -w /data --rm broadinstitute/gatk:4.3.0.0 \
#         gatk CombineGVCFs \
#         ${file_name_arr[@]} \
#         -O 8_vcf_identification/${set_name}.merge.vcf \
#         -R 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta;
# fi
# if [[ ! -f 8_vcf_identification/${set_name}.merge.vcf ]]; then
#     echo "Error."
#     exit 1;
# fi






for element in "${input_name_arr[@]}"; do
    if [[ ! -f 8_vcf_identification/${element}.hg38.identified.vcf.gz ]]; then
        echo "Compressing VCF..."
        docker run \
            -u "$(id -u $USER):$(id -g $USER)" \
            -v /etc/passwd:/etc/passwd:ro \
            -v /etc/group:/etc/group:ro \
            --name nakamae_vaf_calculation --memory 120g -itv $PWD:/data -w /data --rm staphb/bcftools:1.16 \
            bcftools view -Oz 8_vcf_identification/${element}.hg38.identified.vcf \
            -o 8_vcf_identification/${element}.hg38.identified.vcf.gz;
        echo "Indexing VCF..."
        docker run \
            -u "$(id -u $USER):$(id -g $USER)" \
            -v /etc/passwd:/etc/passwd:ro \
            -v /etc/group:/etc/group:ro \
            --name nakamae_vaf_calculation --memory 120g -itv $PWD:/data -w /data --rm biocontainers/tabix:v1.9-11-deb_cv1 \
            tabix -p vcf 8_vcf_identification/${element}.hg38.identified.vcf.gz;
    fi
done
if [[ ! -f 8_vcf_identification/${set_name}.merge.vcf ]]; then
    echo "Merge VCF..."
    docker run \
        -u "$(id -u $USER):$(id -g $USER)" \
        -v /etc/passwd:/etc/passwd:ro \
        -v /etc/group:/etc/group:ro \
        --name nakamae_vaf_calculation --memory 120g -itv $PWD:/data -w /data --rm staphb/bcftools:1.16 \
        bcftools merge --merge all \
        ${file_name_arr[@]} \
        -O v -o 8_vcf_identification/${set_name}.merge.vcf;
fi
if [[ ! -f 8_vcf_identification/${set_name}.merge.vcf ]]; then
    echo "Error."
    exit 1;
fi

if [[ ! -f 8_vcf_identification/${set_name}.hg38.identified.vcf ]]; then
    echo "Intersect VCF..."

    # Build JEXL: all samples must be non-reference (not HomRef)
    jexl=""
    for s in "${input_name_arr[@]}"; do
        cond="!vc.getGenotype('${s}').isHomRef()"
        if [[ -z "${jexl}" ]]; then
            jexl="${cond}"
        else
            jexl="${jexl} && ${cond}"
        fi
    done
    echo "[INFO] JEXL = ${jexl}"

    docker run \
    -u "$(id -u $USER):$(id -g $USER)" \
    -v /etc/passwd:/etc/passwd:ro \
    -v /etc/group:/etc/group:ro \
    --memory 120g -itv "$PWD:/data" -w /data --rm broadinstitute/gatk:4.3.0.0 \
    gatk SelectVariants \
        -V "8_vcf_identification/${set_name}.merge.vcf" \
        -O "8_vcf_identification/${set_name}.hg38.identified.vcf" \
        --max-nocall-number 0 \
        -select "${jexl}" \
        -R "5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta"
fi
if [[ ! -f 8_vcf_identification/${set_name}.hg38.identified.vcf ]]; then
    echo "Error."
    exit 1;
fi

echo "Done."
exit 0;