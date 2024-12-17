#!/bin/bash
input_name=$1
output=$2
result_name=$3

if [[ ! -e $output ]]; then
      echo ${output}" cannot be found.";
      exit 1;
fi
cd $output;
[ ! -d report ] && mkdir report && echo "Create report"};

mkdir report/${result_name};
# All
mkdir report/${result_name}/all;
mkdir report/${result_name}/all/vcf;
cp 10_snp_hard_filter/${input_name}.hg38.identified.snp.fltr.vcf report/${result_name}/all/vcf;
cp -r 11_snp_classification report/${result_name}/all/sequence;
cp -r 12_motif_extraction report/${result_name}/all/motif;
cp 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed/stats report/${result_name}/all/summary.txt;
cp -r 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed/multiqc_data report/${result_name}/all/figure;
# Each
mkdir report/${result_name}/each;
for condition in ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.*vaf
  do
  mkdir report/${result_name}/each/${condition};
  mkdir report/${result_name}/each/${condition}/vcf;
  cp 14_vaf_classification/${condition}.vcf report/${result_name}/each/condition/vcf/;
  mkdir report/${result_name}/each/${condition}/sequence;
  cp ${condition}/*.fa report/${result_name}/each/sequence/;
  mkdir report/${result_name}/each/${condition}/motif;
  cp ${condition}/*.png report/${result_name}/each/motif/;
  cp ${condition}/*.txt report/${result_name}/each/motif/;
  cp 14_vaf_classification/${condition}/stats report/${result_name}/each/summary.txt;
  cp -r 14_vaf_classification/${condition}/multiqc_data report/${result_name}/each/figure;
done
tar -zcvf report/${result_name}.tar.gz report/${result_name};
rm -rf report/${result_name};

exit 0;