# Evaluation Script Usage

## Precautions:
Our evaluation process adopts the workflow of [ChatDrug](https://github.com/chao1224/ChatDrug) for evaluating small molecules. Therefore, the evaluation code is based on the script provided by this work. Please note that, unlike this work, our success rate calculation method is `the number of successful optimizations` / `the total number of molecules`, while this work is `the number of successful optimizations` / `the number of valid molecules`.  
In addition, in our optimization tasks involving "esol", "bbbp", and "herg" properties, the property values are given by the idrug model. Currently, the API call of the idrug model is not open to the public, but the [web query](https://drug.ai.tencent.com/cn) is free for everyone. You can query the aforementioned property values of molecules through this method. If you want to evaluate the above tasks in batches, please implement your own property prediction model at the code entry we reserved.

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
