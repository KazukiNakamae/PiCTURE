#!/bin/bash
set -euo pipefail

if [ $# -lt 4 ]; then
  echo "Usage: $0 <working_directory> <set_name> <sample_label_1> <sample_label_2> [sample_label_3 ...]"
  exit 1
fi

output=$1
set_name=$2

# sample labels (used for filenames and (after reheader) sample names in VCF)
args=("$@")
input_name_arr=("${args[@]:2}")

if [[ ! -e "${output}" ]]; then
  echo "[ERROR] ${output} cannot be found."
  exit 1
fi

cd "${output}"

# -------------------------------
# 1) Ensure per-sample VCFs are indexed/compressed
#    Expected inputs:
#      8_vcf_identification/<label>.hg38.identified.vcf
# -------------------------------
for element in "${input_name_arr[@]}"; do
  if [[ ! -f "8_vcf_identification/${element}.hg38.identified.vcf" ]]; then
    echo "[ERROR] Missing input VCF: 8_vcf_identification/${element}.hg38.identified.vcf"
    echo "        (This script assumes per-sample identified VCFs already exist.)"
    exit 1
  fi

  if [[ ! -f "8_vcf_identification/${element}.hg38.identified.vcf.idx" ]]; then
    echo "Indexing VCF (GATK IndexFeatureFile): ${element}"
    docker run --rm -it       -u "$(id -u "${USER}")":"$(id -g "${USER}")"       -v /etc/passwd:/etc/passwd:ro       -v /etc/group:/etc/group:ro       -v "${PWD}":/data -w /data       broadinstitute/gatk:4.3.0.0       gatk IndexFeatureFile         -I "8_vcf_identification/${element}.hg38.identified.vcf"
  fi

  if [[ ! -f "8_vcf_identification/${element}.hg38.identified.vcf.gz" ]]; then
    echo "Compressing VCF (bgzip via bcftools): ${element}"
    docker run --rm -it       -u "$(id -u "${USER}")":"$(id -g "${USER}")"       -v /etc/passwd:/etc/passwd:ro       -v /etc/group:/etc/group:ro       -v "${PWD}":/data -w /data       staphb/bcftools:1.16       bcftools view -Oz         "8_vcf_identification/${element}.hg38.identified.vcf"         -o "8_vcf_identification/${element}.hg38.identified.vcf.gz"

    echo "Indexing VCF (tabix): ${element}"
    docker run --rm -it       -u "$(id -u "${USER}")":"$(id -g "${USER}")"       -v /etc/passwd:/etc/passwd:ro       -v /etc/group:/etc/group:ro       -v "${PWD}":/data -w /data       biocontainers/tabix:v1.9-11-deb_cv1       tabix -p vcf "8_vcf_identification/${element}.hg38.identified.vcf.gz"
  fi
done

# -------------------------------
# 2) Reheader (rename) samples so that:
#    - bcftools merge does NOT hit "Duplicate sample names"
#    - JEXL that refers to sample labels is valid
#    This is critical when the VCF header sample name does not match your label.
# -------------------------------
merge_inputs=()

for element in "${input_name_arr[@]}"; do
  in_vcf="8_vcf_identification/${element}.hg38.identified.vcf.gz"

  # Query sample names in the VCF
  vcf_samples=$(docker run --rm -i     -u "$(id -u "${USER}")":"$(id -g "${USER}")"     -v /etc/passwd:/etc/passwd:ro     -v /etc/group:/etc/group:ro     -v "${PWD}":/data -w /data     staphb/bcftools:1.16     bcftools query -l "${in_vcf}" | tr -d '\r')

  n_samples=$(printf "%s\n" "${vcf_samples}" | sed '/^$/d' | wc -l | awk '{print $1}')

  if [[ "${n_samples}" -ne 1 ]]; then
    echo "[WARN] ${in_vcf} contains ${n_samples} samples (expected 1). Skipping reheader."
    merge_inputs+=("${in_vcf}")
    continue
  fi

  old_name=$(printf "%s" "${vcf_samples}" | head -n 1)

  if [[ "${old_name}" == "${element}" ]]; then
    merge_inputs+=("${in_vcf}")
    continue
  fi

  echo "[INFO] Reheader: ${old_name} -> ${element} (${in_vcf})"
  name_file="8_vcf_identification/${element}.new_sample_name.txt"
  printf "%s\n" "${element}" > "${name_file}"

  out_vcf="8_vcf_identification/${element}.hg38.identified.renamed.vcf.gz"
  if [[ ! -f "${out_vcf}" ]]; then
    docker run --rm -it       -u "$(id -u "${USER}")":"$(id -g "${USER}")"       -v /etc/passwd:/etc/passwd:ro       -v /etc/group:/etc/group:ro       -v "${PWD}":/data -w /data       staphb/bcftools:1.16       bcftools reheader -s "${name_file}" "${in_vcf}" -Oz -o "${out_vcf}"

    docker run --rm -it       -u "$(id -u "${USER}")":"$(id -g "${USER}")"       -v /etc/passwd:/etc/passwd:ro       -v /etc/group:/etc/group:ro       -v "${PWD}":/data -w /data       biocontainers/tabix:v1.9-11-deb_cv1       tabix -f -p vcf "${out_vcf}"
  fi

  merge_inputs+=("${out_vcf}")
done

# -------------------------------
# 3) Merge into multi-sample VCF
# -------------------------------
if [[ ! -f "8_vcf_identification/${set_name}.merge.vcf" ]]; then
  echo "Merge VCF..."
  docker run --rm -it     -u "$(id -u "${USER}")":"$(id -g "${USER}")"     -v /etc/passwd:/etc/passwd:ro     -v /etc/group:/etc/group:ro     -v "${PWD}":/data -w /data     staphb/bcftools:1.16     bcftools merge --merge all       "${merge_inputs[@]}"       -O v -o "8_vcf_identification/${set_name}.merge.vcf"
fi

if [[ ! -f "8_vcf_identification/${set_name}.merge.vcf" ]]; then
  echo "[ERROR] Merge failed."
  exit 1
fi

# -------------------------------
# 4) Intersection using Method C (JEXL on merged multi-sample VCF)
#    Recommended robust expression:
#      - vc.hasGenotypes(): sanity check
#      - getNoCallCount()==0 : all samples called
#      - getHomRefCount()==0 : no sample is HomRef => all samples are non-reference (Het/HomVar)
# -------------------------------
if [[ ! -f "8_vcf_identification/${set_name}.hg38.identified.vcf" ]]; then
  echo "Intersect VCF..."

  jexl="vc.hasGenotypes() && vc.getNoCallCount() == 0 && vc.getHomRefCount() == 0"
  echo "[INFO] JEXL = ${jexl}"

  docker run --rm -it     -u "$(id -u "${USER}")":"$(id -g "${USER}")"     -v /etc/passwd:/etc/passwd:ro     -v /etc/group:/etc/group:ro     -v "${PWD}":/data -w /data     broadinstitute/gatk:4.3.0.0     gatk SelectVariants       -V "8_vcf_identification/${set_name}.merge.vcf"       -O "8_vcf_identification/${set_name}.hg38.identified.vcf"       -select "${jexl}"       -R "5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta"
fi

if [[ ! -f "8_vcf_identification/${set_name}.hg38.identified.vcf" ]]; then
  echo "[ERROR] Intersect step failed."
  exit 1
fi

echo "Done."
