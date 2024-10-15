#!/bin/bash
input_name=$1
output=$2
vaf_threshold=$3

if [[ ! -e $output ]]; then
      echo ${output}" cannot be found.";
      exit 1;
fi
cd $output;


if (( $(echo "${vaf_threshold} > 1.000" | bc -l) )); then
     echo "VAF threshold should be 0.0-1.0.";
     exit 1;
fi
if (( $(echo "${vaf_threshold} < 0.000" | bc -l) )); then
     echo "VAF threshold should be 0.0-1.0.";
     exit 1;
fi

###

# SNP抽出
if [[ ! -f 9_snp_detection/${input_name}.hg38.identified.snp.vcf ]]; then
     echo "SNP detection..."
     singularity exec --bind `pwd`:/data -W /data ../gatk4.sif \
          gatk SelectVariants -R 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
          -V 8_vcf_identification/${input_name}.hg38.identified.vcf \
          --select-type-to-include SNP \
          -O 9_snp_detection/${input_name}.hg38.identified.snp.vcf;
fi

if [[ ! -f 9_snp_detection/${input_name}.hg38.identified.snp.vcf ]]; then
    echo "Error."
    exit 1;
fi
echo "Done"

# SNPフィルタリング
# QD: Variant Confidence/Quality by Depth.
# QUAL: -10*log10 (posterior genotype probability of a homozygous-reference genotype (GT=0/0))
# SOR: Symmetric Odds Ratio of 2x2 contingency table to detect strand bias.
# FS: Phred-scaled p-value using Fisher's exact test to detect strand bias.
# MQ: RMS Mapping Quality.
# MQRankSum: Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities.
# ReadPosRankSum: Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias.
if [[ ! -f 10_snp_hard_filter/${input_name}.hg38.identified.snp.fltr.vcf ]]; then
     echo "SNP filtering..."
     singularity exec --bind `pwd`:/data -W /data ../gatk4.sif \
          gatk VariantFiltration -R 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
          -V 9_snp_detection/${input_name}.hg38.identified.snp.vcf \
          -O 10_snp_hard_filter/${input_name}.hg38.identified.snp.fltr.vcf \
          -filter "QD < 2.0" --filter-name "QD2" \
          -filter "QUAL < 30.0" --filter-name "QUAL30" \
          -filter "SOR > 3.0" --filter-name "SOR3" \
          -filter "FS > 60.0" --filter-name "FS60" \
          -filter "MQ < 40.0" --filter-name "MQ40" \
          -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
          -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8";
fi

if [[ ! -f 10_snp_hard_filter/${input_name}.hg38.identified.snp.fltr.vcf ]]; then
    echo "Error."
    exit 1;
fi
echo "Done"


if [ $(grep -v "#" 10_snp_hard_filter/${input_name}.hg38.identified.snp.fltr.vcf | wc -l) -eq 0 ]; then
     echo "There is no valid SNV in "${input_name};
     exit 0;
fi

if [[ ! -f 11_snp_classification/${input_name}.hg38.identified.snp.fltr.all.fa ]]; then
     echo "SNP classification..."
     singularity exec --bind `pwd`:/data -W / ../pysam.sif \
          python /extract_flanking_seq_from_vcf.py \
          -i 10_snp_hard_filter/${input_name}.hg38.identified.snp.fltr.vcf \
          -r 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
          -o 11_snp_classification/${input_name}.hg38.identified.snp.fltr;
fi

if [[ ! -f 11_snp_classification/${input_name}.hg38.identified.snp.fltr.all.fa ]]; then
    echo "Error."
    exit 1;
fi
echo "Done"

mut_seq_arr=(${input_name}".hg38.identified.snp.fltr.CtoG.fa" \
${input_name}".hg38.identified.snp.fltr.CtoT.fa" \
${input_name}".hg38.identified.snp.fltr.CtoA.fa" \
${input_name}".hg38.identified.snp.fltr.AtoG.fa" \
${input_name}".hg38.identified.snp.fltr.AtoC.fa" \
${input_name}".hg38.identified.snp.fltr.AtoT.fa" \
${input_name}".hg38.identified.snp.fltr.TtoC.fa" \
${input_name}".hg38.identified.snp.fltr.TtoA.fa" \
${input_name}".hg38.identified.snp.fltr.TtoG.fa" \
${input_name}".hg38.identified.snp.fltr.GtoA.fa" \
${input_name}".hg38.identified.snp.fltr.GtoC.fa" \
${input_name}".hg38.identified.snp.fltr.GtoT.fa" \
${input_name}".hg38.identified.snp.fltr.all.fa");


if [[ ! -f 12_motif_extraction/${file}.txt ]]; then
     echo "motif estimation..."
     for file in ${mut_seq_arr[@]}
     do
          echo ${file};
          singularity exec --bind `pwd`:/data -W /data ../motif_extraction.sif \
               weblogo \
               -f 11_snp_classification/${file} \
               -o 12_motif_extraction/${file}.png \
               -F png -A dna -U bits --resolution 600;
          singularity exec --bind `pwd`:/data -W /data ../motif_extraction.sif \
               weblogo \
               -f 11_snp_classification/${file} \
               -o 12_motif_extraction/${file}.txt \
               -F logodata -A dna -U bits;
     done
fi

if [[ ! -f 12_motif_extraction/${file}.txt ]]; then
    echo "Error."
    exit 1;
fi
echo "Done"



# Add VAF
if [[ ! -f 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.vcf ]]; then
     echo "Add VAF..."
     singularity exec --bind `pwd`:/data -W /data ../bcftools.sif \
          bcftools +fill-tags \
          10_snp_hard_filter/${input_name}.hg38.identified.snp.fltr.vcf \
          -Ov -o 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.vcf \
          -- -t FORMAT/VAF;
fi

if [[ ! -f 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.vcf ]]; then
    echo "Error."
    exit 1;
fi
echo "Done"

# Fix header
if [[ ! -f 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.vcf ]]; then
     echo "Fix VCF header..."
     singularity exec --bind $PWD:/data -W /opt/picard ../picard_cwl.sif \
          java -jar /opt/picard/picard.jar FixVcfHeader \
          I=/data/13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.vcf \
          O=/data/13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.vcf;
fi

if [[ ! -f 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.vcf ]]; then
    echo "Error."
    exit 1;
fi
echo "Done"

# Export Summary
[ ! -d 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed ] \
&& mkdir 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed \
&& echo "Create 13_vaf_calculation/"${input_name}".hg38.identified.snp.fltr.vaf.headerfixed";
if [[ ! -f 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed/stats ]]; then
     echo "Export Summary of all SNPs..."
     singularity exec --bind `pwd`:/data -W /data ../bcftools.sif \
          bcftools stats 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.vcf \
          > 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed/stats;
fi

if [[ ! -f 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed/stats ]]; then
    echo "Error."
    exit 1;
fi
echo "Done"

# Export figures
if [[ ! -d 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed/multiqc_data ]]; then
     echo "Visualization of all SNPs..."
     singularity exec --bind `pwd`:/data -W /data ../multiqc.sif \
          multiqc 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed \
          --outdir 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed/multiqc_data;
fi

if [[ ! -d 13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed/multiqc_data ]]; then
    echo "Error."
    exit 1;
fi
echo "Done"

# Extract SNP with VAF>=${vaf_threshold}
if [[ ! -f 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf.vcf ]]; then
     echo "Extract SNP with VAF>="${vaf_threshold}"...";
     singularity exec --bind `pwd`:/data -W /data ../bcftools.sif \
          bcftools view -O v\
          -o 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf.vcf \
          -e 'FORMAT/VAF<'${vaf_threshold} \
          13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.vcf;
fi

if [[ ! -f 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf.vcf ]]; then
    echo "Error."
    exit 1;
fi
echo "Done"


# Export Summary
[ ! -d 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf ] \
&& mkdir 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf \
&& echo "14_vaf_classification/"${input_name}".hg38.identified.snp.fltr.vaf.headerfixed."${vaf_threshold}"_1.0vaf";
if [[ ! -f 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf/stats ]]; then
     echo "Export Summary of SNPs with VAF>="${vaf_threshold}"..."
     singularity exec --bind `pwd`:/data -W /data ../bcftools.sif \
          bcftools stats 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf.vcf \
          > 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf/stats;
fi

if [[ ! -f 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf/stats ]]; then
    echo "Error."
    exit 1;
fi
echo "Done"

# Export figures
if [[ ! -d 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf/multiqc_data ]]; then
     echo "Visualization of SNPs with VAF>="${vaf_threshold}"..."
     singularity exec --bind `pwd`:/data -W /data ../multiqc.sif \
          multiqc 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf \
          --outdir 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf/multiqc_data;
fi

if [[ ! -d 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf/multiqc_data ]]; then
    echo "Error."
    exit 1;
fi
echo "Done"


[ ! -d ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf ] \
&& mkdir ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf \
&& echo "Create "${input_name}".hg38.identified.snp.fltr.vaf.headerfixed."${vaf_threshold}"_1.0vaf";
echo "Extract flanking sequences around SNPs with VAF>="${vaf_threshold}"..."
singularity exec --bind `pwd`:/data -W / ../pysam.sif \
     python /extract_flanking_seq_from_vcf.py \
     -i 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf.vcf \
     -r 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
     -o ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf/${input_name}.hg38.identified.snp.fltr;
echo "Done"

echo "motif estimation for SNPs with VAF>="${vaf_threshold}"..."
for file in ${mut_seq_arr[@]}
do
     echo ${file};
     singularity exec --bind `pwd`:/data -W /data ../motif_extraction.sif \
          weblogo \
          -f ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf/${file} \
          -o ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf/${file}.png \
          -F png -A dna -U bits --resolution 600;
     singularity exec --bind `pwd`:/data -W /data ../motif_extraction.sif \
          weblogo \
          -f ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf/${file} \
          -o ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.${vaf_threshold}_1.0vaf/${file}.txt \
          -F logodata -A dna -U bits;
done
echo "Done"

# When VAF threshold is 0, the blow process is slipped, because the processed data must be empty.
if (( $(echo "${vaf_threshold} > 0" | bc -l) )); then
     # Extract SNP with VAF<${vaf_threshold}
     if [[ ! -f 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf.vcf ]]; then
          singularity exec --bind `pwd`:/data -W /data ../bcftools.sif \
               bcftools view -O v\
               -o 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf.vcf \
               -e 'FORMAT/VAF>='${vaf_threshold} \
               13_vaf_calculation/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.vcf;
     fi

     if [[ ! -f 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf.vcf ]]; then
     echo "Error."
     exit 1;
     fi
     echo "Done"

     # Export Summary
     [ ! -d 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf ] \
     && mkdir 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf \
     && echo "14_vaf_classification/"${input_name}".hg38.identified.snp.fltr.vaf.headerfixed.0.0_"${vaf_threshold}"vaf";
     if [[ ! -f 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf/stats ]]; then
          echo "Export Summary of SNPs with VAF<"${vaf_threshold}"..."
          singularity exec --bind `pwd`:/data -W /data ../bcftools.sif \
               bcftools stats 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf.vcf \
               > 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf/stats;
     fi

     if [[ ! -f 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf/stats ]]; then
     echo "Error."
     exit 1;
     fi
     echo "Done"

     # Export figures
     if [[ ! -d 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf/multiqc_data ]]; then
          echo "Visualization of SNPs with VAF<"${vaf_threshold}"..."
          singularity exec --bind `pwd`:/data -W /data ../multiqc.sif \
               multiqc 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf \
               --outdir 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf/multiqc_data;
     fi

     if [[ ! -d 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf/multiqc_data ]]; then
     echo "Error."
     exit 1;
     fi
     echo "Done"

     [ ! -d ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf ] \
     && mkdir ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf \
     && echo "Create "${input_name}".hg38.identified.snp.fltr.vaf.headerfixed.0.0_"${vaf_threshold}"vaf";
     echo "Extract flanking sequences around SNPs with VAF<"${vaf_threshold}"..."
     singularity exec --bind `pwd`:/data -W / ../pysam.sif \
          python /extract_flanking_seq_from_vcf.py \
          -i 14_vaf_classification/${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf.vcf \
          -r 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
          -o ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf/${input_name}.hg38.identified.snp.fltr;

     echo "motif estimation for SNPs with VAF<"${vaf_threshold}"..."
     for file in ${mut_seq_arr[@]}
     do
          echo ${file};
          singularity exec --bind `pwd`:/data -W /data ../motif_extraction.sif \
               weblogo \
               -f ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf/${file} \
               -o ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf/${file}.png \
               -F png -A dna -U bits --resolution 600;
          singularity exec --bind `pwd`:/data -W /data ../motif_extraction.sif \
               weblogo \
               -f ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf/${file} \
               -o ${input_name}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_${vaf_threshold}vaf/${file}.txt \
               -F logodata -A dna -U bits;
     done
     echo "Done"
fi

exit 0;