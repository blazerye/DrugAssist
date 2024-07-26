import os
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
# from idrug_utils import grpc_client

def cal_none_rdkit_prop(mol, task_name):
    # Please implement your own property prediction model here
    # We are using the iDrug model here

    # val = cal_prop_idrug(mol, task_name)
    # return val
    return None

def cal_prop_idrug(mol, task_name):
    try:
        client = grpc_client.Client(ip='ip address', port='port')
    except:
        print("connection failed")
        return -999

    query_mol = [mol]
    task = [task_name]
    res = client.get_mole_prop(tasks=task, smiles=query_mol)

    for task in res.score:
        try:
            cur_val = res.score[task].task_score[query_mol[0]].val
            cur_val = round(cur_val[0], 4)
            print(f"query result: {task} {cur_val} {mol}")
            return cur_val
        except:
            print("failed query")
            return -999       


props = ["MolLogP", "qed", "TPSA", "NumHAcceptors", "NumHDonors"]
prop_pred = [(n, func) for n, func in Descriptors.descList if n.split("_")[-1] in props]

prop2func = {}
for prop, func in prop_pred:
    prop2func[prop] = func

task_specification_dict_molecule = {
    "qed+": "Help me make molecule with the SMILES string SMILES_PLACEHOLDER more like a drug. The output molecule should be similar to the input molecule.",
    "acceptor+": "Can you make molecule SMILES_PLACEHOLDER with more hydrogen bond acceptors? The output molecule should be similar to the input molecule.",
    "donor+": "Can you make molecule SMILES_PLACEHOLDER with more hydrogen bond donors? The output molecule should be similar to the input molecule.",
    "esol+": "How can we modify the molecule SMILES_PLACEHOLDER to increase its water solubility value while keeping it similar to the input molecule?",
    "bbbp+": "How can we modify the molecule SMILES_PLACEHOLDER to increase its blood-brain barrier penetration (BBBP) value while keeping it similar to the input molecule?",
    "herg-": "How can we modify the molecule SMILES_PLACEHOLDER to decrease its hERG inhibition value while keeping it similar to the input molecule?",
    "esol+acceptor+": "How can we modify the molecule SMILES_PLACEHOLDER to increase its water solubility value and to have more hydrogen bond acceptors? The output molecule should be similar to the input molecule.",
    "qed+bbbp+": "How can we modify the molecule SMILES_PLACEHOLDER to increase its blood-brain barrier penetration (BBBP) value and make it more like a drug? The output molecule should be similar to the input molecule."
}

task_specification_dict_molecule_strict = {
    "qed+": "How can we modify the molecule SMILES_PLACEHOLDER to increase its QED value by at least 0.1 compared to the pre-optimized value to make it more drug-like while keeping it similar to the input molecule?",
    "acceptor+": "Help me increase the number of hydrogen bond acceptors in the molecule SMILES_PLACEHOLDER by at least 2 compared to the pre-optimized value. The output molecule should be similar to the input molecule.",
    "donor+": "Help me increase the number of hydrogen bond donors in the molecule SMILES_PLACEHOLDER by at least 2 compared to the pre-optimized value. The output molecule should be similar to the input molecule.",
    "esol+": "How can we modify the molecule SMILES_PLACEHOLDER to increase its water solubility value while keeping it similar to the input molecule?",
    "bbbp+": "How can we modify the molecule SMILES_PLACEHOLDER to increase its blood-brain barrier penetration (BBBP) value by at least 0.1 compared to the pre-optimized value while keeping it similar to the input molecule?",
    "herg-": "How can we modify the molecule SMILES_PLACEHOLDER to decrease its hERG inhibition value by at least 0.1 compared to the pre-optimized value while keeping it similar to the input molecule?",
    "esol+acceptor+": "How can we modify the molecule SMILES_PLACEHOLDER to increase its water solubility value and to have more hydrogen bond acceptors? The output molecule should be similar to the input molecule.",
    "qed+bbbp+": "How can we modify the molecule SMILES_PLACEHOLDER to increase its blood-brain barrier penetration (BBBP) value by at least 0.1 and increase its QED value by at least 0.1 compared to the pre-optimized value to make it more drug-like? The output molecule should be similar to the input molecule."
}

task2threshold_list = { 
    "qed+": [[0], [0.1]], 
    "acceptor+": [[0], [1]], 
    "donor+": [[0], [1]], 
    "esol+": [[0], [999]], 
    "bbbp+": [[0], [0.1]], 
    "herg-": [[0], [0.1]],  
    "esol+acceptor+": [[0, 0], [999, 1]], 
    "qed": [[0, 0], [0.1, 0.1]], 
}


def parse_molecule(raw_text, retrieval_sequence):
    record = raw_text.strip()
    record = record.replace('\r\n', '\n')
    record = record.replace('\n', ' ')
    if record.endswith('.'):
        record = record[:-1]
    try:
        record = record.split(":")[-1]
    except:
        pass

    if '"' in record and '->' not in record:
        start = record.find('"') + 1
        end = record.find('"', start)
        tmp = record[start:end]
    else:
        if ' ' not in record:
            tmp = record
        else:
            words = record.split(' ')
            try:
                index = words.index('->')
                tmp = words[index + 1]
                tmp = tmp.replace('"', '')
            except ValueError:
                longest_word = max(words, key=len)
                tmp = longest_word

    if  retrieval_sequence != None and tmp == retrieval_sequence:
        print("The generated molecule is the same as the prompt molecule, skip ReDF!")
        return []
    else:
        return [tmp]


def query_local_data(data, mol, prop):
    row = data.get(mol)
    return float(row[prop]) if row else None


def evaluate_molecule(input_SMILES, output_SMILES, task_id, threshold_list=[0]):
    input_mol = Chem.MolFromSmiles(input_SMILES)
    Chem.Kekulize(input_mol)

    try:
        output_mol = Chem.MolFromSmiles(output_SMILES)
        Chem.Kekulize(output_mol)
    except:
        return None, None, -1

    try:
        if output_mol is None:
            return None, None, -1

        elif task_id == "qed+":
            prop = "qed"
            threshold = threshold_list[0]
            input_value = prop2func[prop](input_mol)
            output_value = prop2func[prop](output_mol)
            if input_value > 0.9 and output_value > input_value:
                return input_value, output_value, True
            return input_value, output_value, output_value > input_value + threshold

        elif task_id == "acceptor+":
            prop = "NumHAcceptors"
            threshold = threshold_list[0]
            input_value = prop2func[prop](input_mol)
            output_value = prop2func[prop](output_mol)
            return input_value, output_value, output_value > input_value + threshold

        elif task_id == "donor+":
            prop = "NumHDonors"
            threshold = threshold_list[0]
            input_value = prop2func[prop](input_mol)
            output_value = prop2func[prop](output_mol)
            return input_value, output_value, output_value > input_value + threshold

        elif task_id == "esol+":
            prop = "esol"
            threshold = threshold_list[0]
            input_value = cal_none_rdkit_prop(input_SMILES, prop)
            output_value = cal_none_rdkit_prop(output_SMILES, prop)
            if input_value is None or output_value is None:
                return None, None, -2
            if threshold != 999:  #loose
                return input_value, output_value, output_value > input_value + threshold
            else:  #strict
                if output_value - input_value > 0.5 and output_value - input_value < 1.5:
                    return input_value, output_value, True
                else:
                    return input_value, output_value, False

        elif task_id == "bbbp+":
            prop = "bbbp"
            threshold = threshold_list[0]
            input_value = cal_none_rdkit_prop(input_SMILES, prop)
            output_value = cal_none_rdkit_prop(output_SMILES, prop)
            if input_value is None or output_value is None:
                return None, None, -2
            if input_value > 0.9 and output_value > input_value:
                return input_value, output_value, True
            return input_value, output_value, output_value > input_value + threshold

        elif task_id == "herg-":
            prop = "hergc10"
            threshold = threshold_list[0]
            input_value = cal_none_rdkit_prop(input_SMILES, prop)
            output_value = cal_none_rdkit_prop(output_SMILES, prop)
            if input_value is None or output_value is None:
                return None, None, -2
            if input_value < 0.1 and output_value < input_value:
                return input_value, output_value, True
            return input_value, output_value, output_value + threshold < input_value

        elif task_id == "esol+acceptor+":
            input_value_01, output_value_01, result_01 = evaluate_molecule(input_SMILES, output_SMILES, "esol+", [threshold_list[0]])
            input_value_02, output_value_02, result_02 = evaluate_molecule(input_SMILES, output_SMILES, "acceptor+", [threshold_list[1]])
            if result_01 == -2 or result_02 == -2:
                return None, None, -2
            return (input_value_01, input_value_02), (output_value_01, output_value_02), result_01 and result_02

        elif task_id == "qed+bbbp+":
            input_value_01, output_value_01, result_01 = evaluate_molecule(input_SMILES, output_SMILES, "qed+", [threshold_list[0]])
            input_value_02, output_value_02, result_02 = evaluate_molecule(input_SMILES, output_SMILES, "bbbp+", [threshold_list[1]])
            if result_01 == -2 or result_02 == -2:
                return None, None, -2
            return (input_value_01, input_value_02), (output_value_01, output_value_02), result_01 and result_02

    except:
        return None, None, -1


def evaluate_molecule_local(local_data, input_SMILES, output_SMILES, task_id, threshold_list=[0]):
    if task_id == "qed+":
        prop = "qed"
        threshold = threshold_list[0]
        input_value = query_local_data(local_data, input_SMILES, prop)
        output_value = query_local_data(local_data, output_SMILES, prop)
        return input_value, output_value, output_value > input_value + threshold

    elif task_id == "acceptor+":
        prop = "NumHAcceptors"
        threshold = threshold_list[0]
        input_value = query_local_data(local_data, input_SMILES, prop)
        output_value = query_local_data(local_data, output_SMILES, prop)
        return input_value, output_value, output_value > input_value + threshold

    elif task_id == "donor+":
        prop = "NumHDonors"
        threshold = threshold_list[0]
        input_value = query_local_data(local_data, input_SMILES, prop)
        output_value = query_local_data(local_data, output_SMILES, prop)
        return input_value, output_value, output_value > input_value + threshold

    elif task_id == "esol+":
        prop = "esol"
        threshold = threshold_list[0]
        input_value = query_local_data(local_data, input_SMILES, prop)
        output_value = query_local_data(local_data, output_SMILES, prop)
        if input_value is None or output_value is None:
            return None, None, -2
        if threshold != 999:
            return input_value, output_value, output_value > input_value + threshold
        else:
            if output_value - input_value > 0.5 and output_value - input_value < 1.5:
                return input_value, output_value, True
            else:
                return input_value, output_value, False

    elif task_id == "bbbp+":
        prop = "bbbp"
        threshold = threshold_list[0]
        input_value = query_local_data(local_data, input_SMILES, prop)
        output_value = query_local_data(local_data, output_SMILES, prop)
        if input_value is None or output_value is None:
            return None, None, -2
        if input_value > 0.9 and output_value > input_value:
            return input_value, output_value, True
        return input_value, output_value, output_value > input_value + threshold

    elif task_id == "herg-":
        prop = "hergc10"
        threshold = threshold_list[0]
        input_value = query_local_data(local_data, input_SMILES, prop)
        output_value = query_local_data(local_data, output_SMILES, prop)
        if input_value is None or output_value is None:
            return None, None, -2
        if input_value < 0.1 and output_value < input_value:
            return input_value, output_value, True
        return input_value, output_value, output_value + threshold < input_value

    elif task_id == "esol+acceptor+":
        input_value_01, output_value_01, result_01 = evaluate_molecule_local(local_data, input_SMILES, output_SMILES, "esol+",
                                                                       [threshold_list[0]])
        input_value_02, output_value_02, result_02 = evaluate_molecule_local(local_data, input_SMILES, output_SMILES, "acceptor+",
                                                                       [threshold_list[1]])
        return (input_value_01, input_value_02), (output_value_01, output_value_02), result_01 and result_02

    elif task_id == "qed+bbbp+":
        input_value_01, output_value_01, result_01 = evaluate_molecule_local(local_data, input_SMILES, output_SMILES, "qed+",
                                                                             [threshold_list[0]])
        input_value_02, output_value_02, result_02 = evaluate_molecule_local(local_data, input_SMILES, output_SMILES, "bbbp+",
                                                                             [threshold_list[1]])
        return (input_value_01, input_value_02), (output_value_01, output_value_02), result_01 and result_02
