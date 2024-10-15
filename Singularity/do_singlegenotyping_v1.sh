INPUT_1=$1
INPUT_2=$2
NAME=$3
OUTPUT_DIRECTORY=$4
THRESHOLD=$5
RESULT_DIRECTORY=$6

if [ $# != 6 ]; then
    echo "Number of argumanents must be 6."
    exit 1
else
    echo "Processing..."
fi

if [[ ! ${INPUT_1} ]]; then
    echo "There is no paired fastq files."
    exit 1;
fi

if [[ ! ${INPUT_2} ]]; then
    echo "There is no paired fastq files."
    exit 1;
fi

chmod +x run.sh;
bash run.sh \
${INPUT_1} \
${INPUT_2} \
${NAME} \
${OUTPUT_DIRECTORY};

if [[ ! -f ${OUTPUT_DIRECTORY}/5_recal_data/${NAME}_recalibration_dbsnp_only_plots.csv ]]; then
    echo "run.sh was failed"
    exit 1;
fi

chmod +x variant_identification_from_singledb.sh;
bash variant_identification_from_singledb.sh ${NAME} ${OUTPUT_DIRECTORY};

if [[ ! -f ${OUTPUT_DIRECTORY}/8_vcf_identification/${NAME}.hg38.identified.vcf ]]; then
    echo "variant_identification_from_singledb.sh was failed"
    exit 1;
fi

chmod +x motif_estimation.sh;
bash motif_estimation.sh ${NAME} ${OUTPUT_DIRECTORY} ${THRESHOLD};

chmod +x get_result.sh;
bash get_result.sh ${NAME} ${OUTPUT_DIRECTORY} ${RESULT_DIRECTORY};
