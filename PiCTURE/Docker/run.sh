#!/bin/bash
input_fore=$1
input_rev=$2
input_name=$3
output=$4
memory=$5

if [[ ! -e $output ]]; then
      echo ${output}" cannot be found.";
      exit 1;
fi
cd $output;

### Copy rawdata

echo "Copy raw fastq data"
if [[ ../${input_fore} == *.fq ]]; then
      cp ../${input_fore} 0_rawdata/${input_name}_F.fastq;
      cp ../${input_rev} 0_rawdata/${input_name}_R.fastq;
elif [[ ../${input_fore} == *.FASTQ ]]; then
      cp ../${input_fore} 0_rawdata/${input_name}_F.fastq;
      cp ../${input_rev} 0_rawdata/${input_name}_R.fastq;
elif [[ ../${input_fore} == *.fastq ]]; then
      cp ../${input_fore} 0_rawdata/${input_name}_F.fastq;
      cp ../${input_rev} 0_rawdata/${input_name}_R.fastq;
else
      echo "Input files must be .fastq format."
      exit 1;
fi

### Adapter Trimming and QC filtering

if [[ ! -f 1_trim_galore/${input_name}_F_val_1.fq ]]; then
    echo "Adapter Trimming and QC filtering"
    docker run \
        -u "$(id -u $USER):$(id -g $USER)" \
        -v /etc/passwd:/etc/passwd:ro \
        -v /etc/group:/etc/group:ro \
        --memory ${memory}g --rm -v `pwd`:/DATA -w /DATA -i clinicalgenomics/trim_galore:0.6.7 \
        trim_galore \
        -j 14 --fastqc --trim1 \
        --output_dir 1_trim_galore --paired 0_rawdata/${input_name}_F.fastq 0_rawdata/${input_name}_R.fastq;
fi

if [[ ! -f 1_trim_galore/${input_name}_F_val_1.fq ]]; then
    echo "Error."
    exit 1;
fi

### STAR mapping

if [[ ! -f 2_star/star_${input_name}_SJ.out.tab ]]; then
    echo "STAR Mapping (1/2)"
    docker run \
        -u "$(id -u $USER):$(id -g $USER)" \
        -v /etc/passwd:/etc/passwd:ro \
        -v /etc/group:/etc/group:ro \
        --memory ${memory}g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/star_for_human_gatk:1.0 \
        STAR \
        --genomeDir ../hg38_index --readFilesIn 1_trim_galore/${input_name}_F_val_1.fq 1_trim_galore/${input_name}_R_val_2.fq \
        --runThreadN 14 --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix 2_star/star_${input_name}_;
fi

if [[ ! -f 2_star/star_${input_name}_SJ.out.tab ]]; then
    echo "Error."
    exit 1;
fi

if [[ ! -f 4_bam_preparation/${input_name}.original.bam ]]; then
    echo "STAR Mapping (2/2)"
    docker run \
        -u "$(id -u $USER):$(id -g $USER)" \
        -v /etc/passwd:/etc/passwd:ro \
        -v /etc/group:/etc/group:ro \
        --memory ${memory}g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/star_for_human_gatk:1.0 \
        STAR \
        --genomeDir ../hg38_index \
        --sjdbFileChrStartEnd 2_star/star_${input_name}_SJ.out.tab \
        --readFilesIn 1_trim_galore/${input_name}_F_val_1.fq 1_trim_galore/${input_name}_R_val_2.fq \
        --runThreadN 14 --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix 3_twopass_align/star_${input_name}_;
    cp 3_twopass_align/star_${input_name}_Aligned.sortedByCoord.out.bam 4_bam_preparation/${input_name}.original.bam;
    rm 1_trim_galore/${input_name}_F_val_1.fq;
    rm 1_trim_galore/${input_name}_R_val_2.fq;
fi

if [[ ! -f 4_bam_preparation/${input_name}.original.bam ]]; then
    echo "Error."
    exit 1;
fi

if [[ ! -f 4_bam_preparation/${input_name}.complete.bam ]]; then
    echo "Formating BAM file"
    cd 4_bam_preparation;
    sudo ./bam_preparation_v2.sh resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
        ${input_name}.original.bam \
        ${input_name} \
        illumina \
        ${input_name}.complete.bam \
        ${memory};
    cd ..
fi

if [[ ! -f 4_bam_preparation/${input_name}.complete.bam ]]; then
    echo "Error."
    exit 1;
fi

if [[ ! -f 5_recal_data/${input_name}_dbsnp_only_recal_data.table ]]; then
    echo "Make recall data"
    docker run \
        -u "$(id -u $USER):$(id -g $USER)" \
        -v /etc/passwd:/etc/passwd:ro \
        -v /etc/group:/etc/group:ro \
        --memory ${memory}g -itv $PWD:/data -w /data --rm broadinstitute/gatk:4.3.0.0 \
        gatk BaseRecalibrator \
        -R 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
        -I 4_bam_preparation/${input_name}.complete.bam \
        --known-sites 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf \
        -O 5_recal_data/${input_name}_dbsnp_only_recal_data.table;
fi

if [[ ! -f 5_recal_data/${input_name}_dbsnp_only_recal_data.table ]]; then
    echo "Error."
    exit 1;
fi

if [[ ! -f 5_recal_data/${input_name}.dbsnp_only.BQSR.bam ]]; then
    echo "BQSR"
    docker run \
        -u "$(id -u $USER):$(id -g $USER)" \
        -v /etc/passwd:/etc/passwd:ro \
        -v /etc/group:/etc/group:ro \
        --memory ${memory}g -itv $PWD:/data -w /data --rm broadinstitute/gatk:4.3.0.0 \
        gatk ApplyBQSR \
        -R 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
        -I 4_bam_preparation/${input_name}.complete.bam \
        --bqsr-recal-file 5_recal_data/${input_name}_dbsnp_only_recal_data.table \
        -O 5_recal_data/${input_name}.dbsnp_only.BQSR.bam;
fi

if [[ ! -f 5_recal_data/${input_name}.dbsnp_only.BQSR.bam ]]; then
    echo "Error."
    exit 1;
fi

if [[ ! -f 5_recal_data/${input_name}_dbsnp_only_post_recal_data.table ]]; then
    # BQSR後のrecal_data.table作成
    docker run \
        -u "$(id -u $USER):$(id -g $USER)" \
        -v /etc/passwd:/etc/passwd:ro \
        -v /etc/group:/etc/group:ro \
        --memory ${memory}g -itv $PWD:/data -w /data --rm broadinstitute/gatk:4.3.0.0 \
        gatk BaseRecalibrator \
        -R 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta \
        -I 5_recal_data/${input_name}.dbsnp_only.BQSR.bam \
        --known-sites 5_recal_data/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf \
        -O 5_recal_data/${input_name}_dbsnp_only_post_recal_data.table;
fi

if [[ ! -f 5_recal_data/${input_name}_dbsnp_only_post_recal_data.table ]]; then
    echo "Error."
    exit 1;
fi

if [[ ! -f 5_recal_data/${input_name}_recalibration_dbsnp_only_plots.csv ]]; then
    # 較正レポート
    docker run \
        -u "$(id -u $USER):$(id -g $USER)" \
        -v /etc/passwd:/etc/passwd:ro \
        -v /etc/group:/etc/group:ro \
        --memory ${memory}g -itv $PWD:/data -w /data --rm broadinstitute/gatk:4.3.0.0 \
        gatk AnalyzeCovariates \
        -before 5_recal_data/${input_name}_dbsnp_only_recal_data.table \
        -after 5_recal_data/${input_name}_dbsnp_only_post_recal_data.table \
        -csv 5_recal_data/${input_name}_recalibration_dbsnp_only_plots.csv \
        -plots 5_recal_data/${input_name}_recalibration_dbsnp_only_plots.pdf;
fi

if [[ ! -f 5_recal_data/${input_name}_recalibration_dbsnp_only_plots.csv ]]; then
    echo "Error."
    exit 1;
fi

echo "Clean files..."
rm 0_rawdata/${input_name}_F.fastq;
rm 0_rawdata/${input_name}_R.fastq;
rm 1_trim_galore/${input_name}_F_val_1.fq;
rm 1_trim_galore/${input_name}_R_val_2.fq;
rm 4_bam_preparation/${input_name}.original.bam;
rm 4_bam_preparation/${input_name}.original.bai;
rm 4_bam_preparation/${input_name}.original.bam.addRG.bam;
rm 4_bam_preparation/${input_name}.original.bam.addRG.bai;
rm 4_bam_preparation/${input_name}.original.bam.addRG.duprm.bam;
rm 4_bam_preparation/${input_name}.original.bam.addRG.duprm.bai;

echo "Done."
exit 0;
