num_gpus=8
lr=1e-4
dtype=bf16
per_device_train_batch_size=8
per_device_eval_batch_size=8
gradient_accumulation_steps=8
num_train_epochs=10
save_total_limit=10
eval_steps=1000
save_steps=2000
lora_rank=64
lora_alpha=128
lora_dropout=0.05
lora_trainable="q_proj,v_proj,k_proj,o_proj,gate_proj,down_proj,up_proj"
modules_to_save="embed_tokens,lm_head"
load_best_model_at_end=False

deepspeed_config_file=config/ds_zero1_no_offload.json

pretrained_model=path/to/llama2-7b-chat
tokenizer_path=path/to/llama2-7b-chat/tokenizer

dataset_dir=path/to/molopt-instructions
train_file=${dataset_dir}/train.json
validation_file=${dataset_dir}/val.json
data_cache_dir=${dataset_dir}/data_cache

output_dir=drugassist_output

torchrun --nnodes 1 --nproc_per_node ${num_gpus} run_sft_lora.py \
    --deepspeed ${deepspeed_config_file} \
    --model_name_or_path ${pretrained_model} \
    --tokenizer_name_or_path ${tokenizer_path} \
    --train_file ${train_file} \
    --validation_file ${validation_file} \
    --data_cache_dir ${data_cache_dir} \
    --per_device_train_batch_size ${per_device_train_batch_size} \
    --per_device_eval_batch_size ${per_device_eval_batch_size} \
    --do_train \
    --do_eval \
    --seed 42 \
    --${dtype} \
    --num_train_epochs ${num_train_epochs} \
    --lr_scheduler_type cosine \
    --learning_rate ${lr} \
    --warmup_ratio 0.03 \
    --weight_decay 0 \
    --logging_strategy steps \
    --logging_steps 10 \
    --save_strategy steps \
    --save_total_limit ${save_total_limit} \
    --evaluation_strategy steps \
    --eval_steps ${eval_steps} \
    --save_steps ${save_steps} \
    --gradient_accumulation_steps ${gradient_accumulation_steps} \
    --preprocessing_num_workers 8 \
    --max_seq_length 1024 \
    --output_dir ${output_dir} \
    --overwrite_output_dir \
    --ddp_timeout 30000 \
    --logging_first_step True \
    --lora_rank ${lora_rank} \
    --lora_alpha ${lora_alpha} \
    --trainable ${lora_trainable} \
    --modules_to_save ${modules_to_save} \
    --lora_dropout ${lora_dropout} \
    --torch_dtype float16 \
    --gradient_checkpointing \
    --ddp_find_unused_parameters False \
    --load_best_model_at_end ${load_best_model_at_end}