output=$1

echo "Check files"
if [[ ! -f resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta ]]; then
    echo "resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta chould be in the current directory."
    exit 1;
fi

echo "Make directories"
mkdir ${output}
mkdir ${output}/0_rawdata;
mkdir ${output}/1_trim_galore;
mkdir ${output}/2_star;
mkdir ${output}/3_twopass_align;
mkdir ${output}/4_bam_preparation;
mkdir ${output}/5_recal_data;
mkdir ${output}/6_haplotypecaller;

echo "Download Docker images"
docker pull clinicalgenomics/trim_galore:0.6.7;
docker pull broadinstitute/gatk:4.3.0.0;

echo "Create STAR index in docker image"
docker build -t kazukinakamae/star_for_human_gatk:1.0 .;


cat << EOT >> ${output}/4_bam_preparation/bam_preparation_v2.sh
# bam_preparation.sh <reference.fa> <mapped.bam> <readname> <platform> <output.bam>
# update read group
echo "Update Read Group..."
docker run \
     -u "\$(id -u \$USER):\$(id -g \$USER)" \
     -v /etc/passwd:/etc/passwd:ro \
     -v /etc/group:/etc/group:ro \
     --name add_rg --memory 120g -itv \$PWD:/data -w /data --rm biocontainers/picard:2.3.0 \
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
# mark duplication
echo "Mark duplication..."
docker run \
     -u "\$(id -u \$USER):\$(id -g \$USER)" \
     -v /etc/passwd:/etc/passwd:ro \
     -v /etc/group:/etc/group:ro \--name mark_dup --memory 120g -itv \$PWD:/data -w /data --rm broadinstitute/gatk:4.3.0.0 \
     gatk MarkDuplicatesSpark \
     -I "\$2".addRG.bam \
     -O "\$2".addRG.duprm.bam \
     -M MarkDuplicatesSpark_output.metrics.txt;

# Splits reads that contain Ns in their cigar string
echo "Split reads that contain Ns in their cigar string..."
docker run \
     -u "\$(id -u \$USER):\$(id -g \$USER)" \
     -v /etc/passwd:/etc/passwd:ro \
     -v /etc/group:/etc/group:ro \
     --name mark_dup --memory 120g -itv \$PWD:/data -w /data --rm broadinstitute/gatk:4.3.0.0 \
     gatk SplitNCigarReads \
     -R \$1 \
     -I "\$2".addRG.duprm.bam \
     -O \$5;
echo "DONE"
EOT
chmod +x ${output}/4_bam_preparation/bam_preparation_v2.sh;

