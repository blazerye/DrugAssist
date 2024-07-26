import pandas as pd
import csv
import time
from rdkit.Chem import AllChem
from rdkit import Chem, DataStructs
from small_molecule_editing import task_specification_dict_molecule_strict, evaluate_molecule, evaluate_molecule_local, task_specification_dict_molecule, task2threshold_list, parse_molecule


def construct_PDDS_prompt(task_specification_dict, input_drug, task, mode, file_name):
    if (mode == 'strict'):
        if (task == "esol+" or task == "esol+acceptor+"):
            data = pd.read_csv(file_name)
            filtered_data = data[data['smiles'] == input_drug]

            bound_l_inc = str(round(float(filtered_data['esol'].iloc[0]) + 0.5, 2))
            bound_h_inc = str(round(float(filtered_data['esol'].iloc[0]) + 1.5, 2))
            bound_l_dec = str(round(float(filtered_data['esol'].iloc[0]) - 1.5, 2))
            bound_h_dec = str(round(float(filtered_data['esol'].iloc[0]) - 0.5, 2))

            if task == "esol+":
                task_prompt_template = "Can you give me an optimized version of the molecule [SMILES] with a water solubility value ranging from [BOUNDL] to [BOUNDH] (logarithm of mol/L) while maintaining similarity to the original molecule?"
                prompt = task_prompt_template.replace('[SMILES]', input_drug)
                prompt = prompt.replace('[BOUNDL]', bound_l_inc)
                prompt = prompt.replace('[BOUNDH]', bound_h_inc)
            elif task == "esol+acceptor+":
                task_prompt_template = "Can you give me an optimized version of the molecule [SMILES] with a water solubility value ranging from [BOUNDL] to [BOUNDH] (logarithm of mol/L), and with at least 2 more hydrogen bond acceptors while maintaining similarity to the original molecule?"
                prompt = task_prompt_template.replace('[SMILES]', input_drug)
                prompt = prompt.replace('[BOUNDL]', bound_l_inc)
                prompt = prompt.replace('[BOUNDH]', bound_h_inc)

            prompt = prompt + " Give me only the SMILES string of the output molecule. No explanation is needed."
            return prompt

    task_prompt_template = task_specification_dict[task]
    prompt = task_prompt_template.replace('SMILES_PLACEHOLDER', input_drug)
    prompt = prompt + " Give me only the SMILES string of the output molecule. No explanation is needed."
        
    return prompt


def retrieve_and_feedback(query_type, task, query_db, DB, input_drug, generated_drug, constraint):
    candidate_mol = []
    gen_mol = Chem.MolFromSmiles(generated_drug)
    gen_fp = AllChem.GetMorganFingerprintAsBitVect(gen_mol, 2, nBits=2048, useChirality=False)

    tic = time.time()
    for db_smi, fp in DB.items():
        answer = evaluate(query_db, query_type, input_drug, db_smi, task, constraint)
        if answer:
            candidate_mol.append([db_smi, fp])
    toc = time.time()
    print(f"success retrieval in {toc - tic:.2f}s")


    if len(candidate_mol) == 0:
        raise Exception("Sorry, Cannot fined a good one")
    else:
        max_sim = 0
        res_smi = candidate_mol[0][0]
        for mol in candidate_mol:
            cur_sim = cal_sim_withfp(gen_fp, mol[1])
            if cur_sim > max_sim:
                max_sim = cur_sim
                res_smi = mol[0]
        print(f"successful retrieval with similarity {max_sim:.2f}")
        return res_smi, max_sim
    

def complete(clt, cur_message, history, api_name="/drugassist"):
    try:
        ans = clt.predict(cur_message, history, api_name=api_name)
        print("User: " + cur_message)
        print("Assistant: " + ans[-1][1])

        return ans[-1][0], ans[-1][1], ans

    except:
        raise Exception('failed to retrieve answer')
    

def evaluate(query_db, query_type, input_drug, generated_drug, task, constraint):
    if len(generated_drug) > 200:
        return -1
    
    if constraint == 'loose':
        threshold_list = task2threshold_list[task][0]
    else:
        threshold_list = task2threshold_list[task][1]

    if query_type == 'local':
        _, _, answer = evaluate_molecule_local(query_db, input_drug, generated_drug, task,
                                                threshold_list=threshold_list)
    else:
        _, _, answer = evaluate_molecule(input_drug, generated_drug, task, threshold_list=threshold_list)
    
    return answer


def get_task_specification_dict(mode):
    if mode == 'loose':
        return task_specification_dict_molecule
    else:
        return task_specification_dict_molecule_strict


def load_dataset(file_path):
    smiles_list = []
    with open(file_path, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            smiles_list.append(row['smiles'])
    return smiles_list


def similarity(a, b):
    if a is None or b is None:
        return 0.0
    amol = Chem.MolFromSmiles(a)
    bmol = Chem.MolFromSmiles(b)
    if amol is None or bmol is None:
        return 0.0

    fp1 = AllChem.GetMorganFingerprintAsBitVect(amol, 2, nBits=2048, useChirality=False)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(bmol, 2, nBits=2048, useChirality=False)
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def cal_sim_withfp(fp1, fp2):
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def load_local_data(file_name):
    data = {}
    with open(file_name, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            smiles = row['smiles']
            data[smiles] = row
    return data


def cal_mol_fp(file_name):
    data = {}
    with open(file_name, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            smiles = row['smiles']
            mol = Chem.MolFromSmiles(smiles)
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048, useChirality=False)
            data[smiles] = fp
    return data


def parse(generated_text, addition_drug=None):
    return parse_molecule(generated_text, addition_drug)

