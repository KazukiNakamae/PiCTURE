#!/bin/bash
output=$1

echo "Make directories"
[ ! -d ${output} ] && mkdir ${output} && echo "Create "${output};
[ ! -d ${output}/0_rawdata ] && mkdir ${output}/0_rawdata && echo "Create "${output}"/0_rawdata";
[ ! -d ${output}/1_trim_galore ] && mkdir ${output}/1_trim_galore && echo "Create "${output}"/1_trim_galore";
[ ! -d ${output}/2_star ] && mkdir ${output}/2_star && echo "Create "${output}"/2_star";
[ ! -d ${output}/3_twopass_align ] && mkdir ${output}/3_twopass_align && echo "Create "${output}"/3_twopass_align";
[ ! -d ${output}/4_bam_preparation ] && mkdir ${output}/4_bam_preparation && echo "Create "${output}"/4_bam_preparation";
[ ! -d ${output}/5_recal_data ] && mkdir ${output}/5_recal_data && echo "Create "${output}"/5_recal_data";
[ ! -d ${output}/6_haplotypecaller ] && mkdir ${output}/6_haplotypecaller && echo "Create "${output}"/6_haplotypecaller";

[ ! -d ${output}/7a_single_genomics_dbImport ] && mkdir ${output}/7a_single_genomics_dbImport && echo "Create "${output}"/7a_single_genomics_dbImport";
[ ! -d ${output}/7b_joint_genomics_dbImport ] && mkdir ${output}/7b_joint_genomics_dbImport && echo "Create "${output}"/7b_joint_genomics_dbImport";
[ ! -d ${output}/single_db ] && mkdir ${output}/single_db && echo "Create "${output}"/single_db";
[ ! -d ${output}/joint_db ] && mkdir ${output}/joint_db && echo "Create "${output}"/joint_db";
[ ! -d ${output}/8_vcf_identification ] && mkdir ${output}/8_vcf_identification && echo "Create "${output}"/8_vcf_identification";

[ ! -d ${output}/9_snp_detection ] && mkdir ${output}/9_snp_detection && echo "Create "${output}"/9_snp_detection";
[ ! -d ${output}/10_snp_hard_filter ] && mkdir ${output}/10_snp_hard_filter && echo "Create "${output}"/10_snp_hard_filter";
[ ! -d ${output}/11_snp_classification ] && mkdir ${output}/11_snp_classification && echo "Create "${output}"/11_snp_classification";
[ ! -d ${output}/12_motif_extraction ] && mkdir ${output}/12_motif_extraction && echo "Create "${output}"/12_motif_extraction";
[ ! -d ${output}/13_vaf_calculation ] && mkdir ${output}/13_vaf_calculation && echo "Create "${output}"/13_vaf_calculation";
[ ! -d ${output}/14_vaf_classification ] && mkdir ${output}/14_vaf_classification && echo "Create "${output}"/14_vaf_classification";

# Download GATK bundle
if [ ! -f ./${output}/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.dict ]; then
     curl -OL https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict;
     mv Homo_sapiens_assembly38.dict ./${output}/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.dict && \
     echo "Locate "${output}"/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.dict";
fi
if [ ! -f ./Homo_sapiens_assembly38.fasta ]; then
     curl -OL https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta;
     echo "Locate Homo_sapiens_assembly38.fasta";
fi
if [ ! -f ./${output}/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta ]; then
     cp Homo_sapiens_assembly38.fasta ./${output}/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta && \
     echo "Locate "${output}"/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta";
fi
if [ ! -f ./${output}/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta.fai ]; then
     curl -OL https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai;
     mv Homo_sapiens_assembly38.fasta.fai ./${output}/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta.fai && \
     echo "Locate "${output}"/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta.fai";
fi
if [ ! -f ./${output}/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf ]; then
     curl -OL https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf;
     mv Homo_sapiens_assembly38.dbsnp138.vcf ./${output}/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf && \
     echo "Locate "${output}"/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf";
fi
if [ ! -f ./${output}/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf.idx ]; then
     curl -OL https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx;
     mv Homo_sapiens_assembly38.dbsnp138.vcf.idx ./${output}/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf.idx && \
     echo "Locate "${output}"/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf.idx";
fi
if [ ! -f ./${output}/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.dict ]; then
     cp ./${output}/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.dict ./${output}/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.dict && \
     echo "Locate "${output}"/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.dict";
fi
if [ ! -f ./${output}/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta.fai ]; then
     cp ./${output}/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta.fai ./${output}/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta.fai && \
     echo "Locate "${output}"/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta.fai";
fi
if [ ! -f ./${output}/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta ]; then
     cp ./${output}/4_bam_preparation/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta ./${output}/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta && \
     echo "Locate "${output}"/5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta";
fi
if [ ! -f ./${output}/7a_single_genomics_dbImport/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list ]; then
     curl -OL https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list;
     mv wgs_calling_regions.hg38.interval_list ./${output}/7a_single_genomics_dbImport/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list && \
     echo "Locate "${output}"/7a_single_genomics_dbImport/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list";
fi
if [ ! -f ./${output}/7b_joint_genomics_dbImport/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list ]; then
     cp ./${output}/7a_single_genomics_dbImport/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list ./${output}/7b_joint_genomics_dbImport/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list && \
     echo "Locate "${output}"/7b_joint_genomics_dbImport/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list";
fi

echo "Build singularity container"
if [ ! -f ./trim_galore.sif ]; then
     singularity build trim_galore.sif docker://clinicalgenomics/trim_galore:0.6.7;
fi
if [ ! -f ./gatk4.sif ]; then
     singularity build gatk4.sif docker://broadinstitute/gatk:4.3.0.0;
fi
if [ ! -f ./motif_extraction.sif ]; then
     singularity build motif_extraction.sif docker://kazukinakamae/motif_extraction:amd64_1.0;
fi
if [ ! -f ./bcftools.sif ]; then
     singularity build bcftools.sif docker://staphb/bcftools:1.16;
fi
if [ ! -f ./picard.sif ]; then
     singularity build picard.sif docker://biocontainers/picard:2.3.0;
fi
if [ ! -f ./picard_cwl.sif ]; then
     singularity build picard_cwl.sif docker://mgibio/picard-cwl:2.18.1;
fi
if [ ! -f ./multiqc.sif ]; then
     singularity build multiqc.sif docker://ewels/multiqc:v1.14;
fi
if [ ! -f ./pysam.sif ]; then
     singularity build pysam.sif docker://kazukinakamae/pysam:amd64_0.19.1;
fi
if [ ! -f ./star.sif ]; then
     singularity build star.sif docker://kazukinakamae/conda_star:2.7.4a;
fi

if [ ! -f ${output}/4_bam_preparation/bam_preparation_v2.sh ]; then
cat << EOT > ${output}/4_bam_preparation/bam_preparation_v2.sh
if [[ ! -f "\$2".addRG.bam ]]; then
	# bam_preparation.sh <reference.fa> <mapped.bam> <readname> <platform> <output.bam>
	# update read group
	echo "Update Read Group..."
	singularity exec --bind \$PWD:/data -W /data ../../picard.sif\
	     picard AddOrReplaceReadGroups \
	     INPUT="\$2" \
	     OUTPUT="\$2".addRG.bam \
	     SO=coordinate \
	     RGID=FLOWCELLID \
	     RGLB="\$3"_library1 \
	     RGPU=unknown \
	     RGPL="\$4" \
	     RGSM="\$3" \
	     CREATE_INDEX=true;
	echo "DONE"
fi
if [[ ! -f "\$2".addRG.bam ]]; then
    echo "Error."
    exit 1;
fi

if [[ ! -f "\$2".addRG.duprm.bam ]]; then
	# mark duplication
	echo "Mark duplication..."
	singularity exec --bind \$PWD:/data -W /data ../../gatk4.sif\
	     gatk MarkDuplicatesSpark \
	     -I "\$2".addRG.bam \
	     -O "\$2".addRG.duprm.bam \
	     -M MarkDuplicatesSpark_output.metrics.txt;
	echo "DONE"
fi
if [[ ! -f "\$2".addRG.duprm.bam ]]; then
    echo "Error."
    exit 1;
fi
if [[ ! -f \$5 ]]; then
	# Splits reads that contain Ns in their cigar string
	echo "Split reads that contain Ns in their cigar string..."
	singularity exec --bind \$PWD:/data -W /data ../../gatk4.sif\
	     gatk SplitNCigarReads \
	     -R \$1 \
	     -I "\$2".addRG.duprm.bam \
	     -O \$5;
	echo "DONE"
fi
if [[ ! -f \$5 ]]; then
    echo "Error."
    exit 1;
fi
EOT
fi
chmod +x ${output}/4_bam_preparation/bam_preparation_v2.sh;


