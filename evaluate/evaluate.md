# Evaluation Script Usage

## Precautions:
Our evaluation process adopts the workflow of [ChatDrug](https://github.com/chao1224/ChatDrug) for evaluating small molecules. Therefore, the evaluation code is based on the script provided by this work. Please note that, unlike this work, our success rate calculation method is `the number of successful optimizations` / `the total number of molecules`, while this work is `the number of successful optimizations` / `the number of valid molecules`.  
In addition, in our optimization tasks involving "esol", "bbbp", and "herg" properties, the property values are given by the idrug model. Currently, the API call of the idrug model is not open to the public, but the [web query](https://drug.ai.tencent.com/cn) is available for everyone. You can query the aforementioned property values of molecules through this method. If you want to evaluate the above tasks in batches, please implement your own property prediction model at the code entry we reserved.

## Script Parameter Explanation:
`--task`: Task name, currently supported values are: `qed+`, `acceptor+`, `donor+`, `esol+`, `bbbp+`, `herg-`, `esol+acceptor+`, `qed+bbbp+`.  
`--log_file`: Log file to record the experiment output.  
`--record_file`: JSON file to record the experiment output.  
`--constraint`: Optimization standard is loose or strict, optional values are `loose` and `strict`, default is `loose`.  
`--num_round`: Maximum number of rounds the model is allowed to modify the result based on feedback, default is `2`.  
`--client_add`: IP and port when the model is deployed through gradio.  
`--database`: The molecular dataset used to prompt the model, default is the file `mol_DB.csv` which we've provided with the script.  
`--test_mol`: Test set molecules, default is the file `testset.csv` which we've provided with the script.

## Usage Examples:
Evaluate the model's performance in optimizing QED values, with strict optimization criteria, and allowing the model to modify the result up to 3 times:
```
python main_eval.py --task qed+ --log_file DrugAssist_qed+_strict.log --record_file DrugAssist_qed+_strict.json --constraint strict  --num_round 3 --client_add 127.0.0.1:19324
```

Evaluate the model's performance in optimizing esol and acceptor simultaneously, with loose optimization criteria, and allowing the model to modify the result up to 2 times:
```
python main_eval.py --task esol+acceptor+ --log_file DrugAssist_esol+acc+_loose.log --record_file DrugAssist_esol+acc+_loose.json --client_add 127.0.0.1:19324
```

# Evaluation of Model Performance under Different Quantization Precisions

To understand the quality loss brought by quantization, we tested the performance of DrugAssist in its 8-bit and 4-bit quantized versions. The table below shows the results. In each data cell, the left side represents the success rate, and the right side represents the validity rate.

| Precision         | qed+_l | qed+_s | donor+_l | donor+_s | acceptor+_l | acceptor+_s | solubility+_l | solubility+_s | bbbp+_l | bbbp+_s | herg-_l | herg-_s | sol+ & acc+_l | sol+ & acc+_s | sol+ & acc+_l | qed+ & bbbp+_s |
|-------------------|--------|--------|----------|----------|-------------|-------------|---------------|---------------|---------|---------|---------|---------|---------------|---------------|---------------|----------------|
| full              | 0.76/0.99 | 0.63/0.97 | 0.72/0.98 | 0.76/0.95 | 0.71/0.97 | 0.67/0.96 | 0.80/0.98 | 0.41/0.98 | 0.82/0.99 | 0.61/0.98 | 0.71/0.99 | 0.67/0.98 | 0.50/0.95 | 0.27/0.95 | 0.65/0.99 | 0.41/0.98 |
| 8-bit quantized   | 0.75/0.98 | 0.61/0.93 | 0.69/0.96 | 0.74/0.93 | 0.69/0.96 | 0.65/0.94 | 0.78/0.96 | 0.36/0.93 | 0.81/0.98 | 0.60/0.96 | 0.69/0.97 | 0.64/0.94 | 0.47/0.93 | 0.23/0.93 | 0.62/0.96 | 0.38/0.96 |
| 4-bit quantized   | 0.73/0.98 | 0.60/0.93 | 0.68/0.95 | 0.70/0.90 | 0.68/0.94 | 0.61/0.91 | 0.77/0.96 | 0.35/0.91 | 0.81/0.98 | 0.56/0.93 | 0.68/0.96 | 0.62/0.92 | 0.43/0.90 | 0.22/0.92 | 0.60/0.95 | 0.35/0.94 |