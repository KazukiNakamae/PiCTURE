#!/bin/bash

#$ -cwd
#$ -V 
#$ -l gpu
#$ -l cuda=4
#$ -l d_rt=100:00:00
#$ -l s_rt=100:00:00
#$ -l s_vmem=384G
#$ -l mem_req=384G
#$ -N rAPOBEC1_dataset_v1_union_40bp
#$ -S /bin/bash

source /home/nakamae_rois_ds/miniconda3/etc/profile.d/conda.sh
conda activate dna
cd /home/nakamae_rois_ds/DNABERT_2/finetune/

export DATA_PATH=/home/nakamae_rois_ds/rAPOBEC1_dataset_v1_union_40bp_balance_kx_tb8_eb16_gas8_lr2e5_te8
export MAX_LENGTH=10
export NAKAMAE_NAME_TAG=balance_kx_tb8_eb16_gas8_lr2e5_te8

export LR=2e-5
python train.py \
    --model_name_or_path zhihan1996/DNABERT-2-117M \
    --data_path  ${DATA_PATH} \
    --kmer -1 \
    --run_name DNABERT2_${DATA_PATH} \
    --model_max_length ${MAX_LENGTH} \
    --per_device_train_batch_size 8 \
    --per_device_eval_batch_size 16 \
    --gradient_accumulation_steps 8 \
    --learning_rate ${LR} \
    --num_train_epochs 8 \
    --fp16 True \
    --save_steps 200 \
    --output_dir ft_output_rAPOBEC1_dataset_v1_union_40bp_${NAKAMAE_NAME_TAG} \
    --evaluation_strategy steps \
    --eval_steps 200 \
    --warmup_steps 50 \
    --logging_steps 100 \
    --overwrite_output_dir True \
    --log_level info \
    --find_unused_parameters False \
    --save_model True
