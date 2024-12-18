#!/bin/bash

#$ -cwd
#$ -V 
#$ -l gpu
#$ -l d_rt=100:00:00
#$ -l s_rt=100:00:00
#$ -l s_vmem=80G 
#$ -l mem_req=80G
#$ -N ft_dnabert2_20240119_2
#$ -S /bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate dna
cd /home/m216692/DNABERT_2/finetune/

export DATA_PATH=/home/m216692/dnabert2_ft_data
export MAX_LENGTH=10

export LR=3e-5
python train.py \
    --model_name_or_path zhihan1996/DNABERT-2-117M \
    --data_path  ${DATA_PATH} \
    --kmer -1 \
    --run_name DNABERT2_${DATA_PATH} \
    --model_max_length ${MAX_LENGTH} \
    --per_device_train_batch_size 10 \
    --per_device_eval_batch_size 16 \
    --gradient_accumulation_steps 1 \
    --learning_rate ${LR} \
    --num_train_epochs 10 \
    --fp16 \
    --save_steps 200 \
    --output_dir ft_output_20240119_2 \
    --evaluation_strategy steps \
    --eval_steps 200 \
    --warmup_steps 55 \
    --logging_steps 100 \
    --overwrite_output_dir True \
    --log_level info \
    --find_unused_parameters False \
    --save_model True
