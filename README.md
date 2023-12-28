<h1 align="center"> 🐹 DrugAssist  </h1>
<h3 align="center"> A Large Language Model for Molecule Optimization </h3>

<p align="center">
  📃 <a href="https://arxiv.org" target="_blank">Paper</a><br>
</p>

<div align="center">
  <img src="fig/logo.png" width="200">
</div>

## 📌 Contents
- [Install](#install)
- [Dataset](#dataset)
- [Train](#train)
- [Demo](#demo)
- [About](#about)

## 🛠️ Install
1. Clone this repository and navigate to DrugAssist folder
```bash
git clone https://github.com/blazerye/DrugAssist.git
cd DrugAssist
```

2. Install Package
```Shell
conda create -n drugassist python=3.8 -y
conda activate drugassist
pip install -r requirements.txt
```

## 🤗 Dataset
* TODO: coming soon

## 🚆 Train
You can use LoRA to finetune `Llama2-7B-Chat` model on the `MolOpt-Instructions` dataset, the running command is as follows:
```Shell
sh run_sft_lora.sh
```

## 👀 Demo
#### Step 1: Merge model weights
You can merge LoRA weights to generate full model weights using the following command:
```Shell
python merge_model.py \
    --base_model $BASE_MODEL_PATH \
    --lora_model $LORA_MODEL_PATH \
    --output_dir $OUTPUT_DIR \
    --output_type huggingface \
    --verbose
```

#### Step 2: Launch web demo
You can use gradio to launch web demo by running the following command:
```Shell
python gradio_service.py \
    --base_model $FULL_MODEL_PATH \
    --ip $IP \
    --port $PORT
```

## 📝 About
### Citation
If you find DrugAssist useful for your research and applications, please cite using this BibTeX:
```bibtex
@misc{ye2024drugassist,
      title={DrugAssist: A Large Language Model for Molecule Optimiaztion}, 
      author={Ye, Geyan and Cai, Xibao and Lai, Houtim and Wang, Xing and Wang, Longyue and Liu, Wei and Zeng, Xiangxiang},
      publisher={arXiv:2410.03744},
      year={2024},
}
```
### Acknowledgements
We appreciate [LLaMA](https://github.com/facebookresearch/llama), [Chinese-LLaMA-Alpaca-2](https://github.com/ymcui/Chinese-LLaMA-Alpaca-2), [Alpaca](https://crfm.stanford.edu/2023/03/13/alpaca.html), [iDrug](https://drug.ai.tencent.com) and many other related works for their open-source contributions.