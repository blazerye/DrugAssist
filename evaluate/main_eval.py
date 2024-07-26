import json
import argparse
import sys
from gradio_client import Client
from rdkit import RDLogger

from utils import (
    construct_PDDS_prompt,retrieve_and_feedback, load_local_data,
    cal_mol_fp, load_dataset, complete, get_task_specification_dict, evaluate, parse)

RDLogger.DisableLog('rdApp.*')

def conversation(clt, query_db, history, messages, C, round_index, task,
                 input_drug, retrieval_DB, record, logfile, constraint):
    try:
        usr_in, generated_text, history = complete(clt, messages[-1], history)
    except:
        return -1, None, history

    print("----------------", file=logfile)
    print("User:" + usr_in, file=logfile)
    print("Assistant:" + generated_text, file=logfile)

    record[input_drug]['retrieval_conversation'][round_index]['user'] = usr_in
    record[input_drug]['retrieval_conversation'][round_index]['assist'] = generated_text

    if round_index >= 1:
        try:
            closest_drug = record[input_drug]['retrieval_conversation'][round_index - 1]['retrieval_drug']
        except:
            closest_drug = None
    else:
        closest_drug = None

    generated_drug_list = parse(generated_text, closest_drug)

    if generated_drug_list == None:
        record[input_drug]['skip_round'] = round_index
        return -1, None, history
    elif len(generated_drug_list) == 0:
        record[input_drug]['retrieval_conversation'][round_index]['answer'] = 'False'
        return 0, None, history
    else:  
        generated_drug = generated_drug_list[0]
        print(f"Generated molecule: {generated_drug}")
        print(f"Generated Result: {generated_drug}", file=logfile)
        record[input_drug]['retrieval_conversation'][round_index]['generated_drug'] = generated_drug
    
    answer = evaluate(None, 'remote', input_drug, generated_drug, task, constraint)

    if answer == -1:
        print('Invalid SMILES generated, skip ReDF', file=logfile)
        print('Invalid SMILES generated, skip ReDF')
        record[input_drug]['skip_round'] = round_index
        return -1, None, history
    elif answer == -2:
        print('Failed to retrieve molecule property value, skip', file=logfile)
        print('Failed to retrieve molecule property value, skip')
        record[input_drug]['skip_round'] = round_index
        return -2, None, history

    print(f"Evaluation result is: {answer}", file=logfile)
    print(f"Evaluation result is: {answer}")
    record[input_drug]['retrieval_conversation'][round_index]['answer'] = str(answer)

    if answer:
        return 1, generated_drug, history
    else:
        if round_index < C:
            answer, generated_drug = ReDF(query_db, messages, round_index, task, input_drug, generated_drug,
                                          retrieval_DB, record, logfile, constraint)
        return answer, generated_drug, history


def ReDF(query_db, messages, round_index, task, input_drug, generated_drug,
         retrieval_DB, record, logfile, constraint):
    print(f'Start Retrieval {round_index + 1}', file=logfile)

    try:
        closest_drug, _= retrieve_and_feedback('local', task, query_db, retrieval_DB, input_drug, generated_drug, constraint)
    except:
        error = sys.exc_info()
        print(error)
        if error[0] == Exception:
            print('Cannot find a retrieval result.', file=logfile)
            return 0, None
        else:
            print('Invalid drug. Failed to evaluate. Skipped.', file=logfile)
            record[input_drug]['skip_round'] = round_index
            return -1, None

    print(f"Retrieval Result: {closest_drug}", file=logfile)
    record[input_drug]['retrieval_conversation'][round_index]['retrieval_drug'] = closest_drug
    prompt_ReDF = f'Your provided molecule {generated_drug} is not correct. We find a molecule {closest_drug} which is correct and similar to the molecule you provided. Give me only the SMILES string of a new molecule satisfying the requirements. No explanation is needed. '
    messages.append(prompt_ReDF)
    return 0, generated_drug


def main_assit(args):
    f = open(args['log_file'], 'w')
    record = {}

    client = Client(args['client_add'], serialize=False)
    print("connection OK")

    local_data = load_local_data(args['database'])
    print("molecule database loaded")

    prompt_db = cal_mol_fp(args['database'])
    print("finish calculating molecule fingerprint")

    task_specification_dict = get_task_specification_dict(args['constraint'])
    input_drug_list = load_dataset(args['test_mol'])

    print("all preprocess procedure finished")
    num_valid = num_correct = num_all = num_invalid_eval = 0

    for index, input_drug in enumerate(input_drug_list):
        print(f">>Sample {index}", file=f)
        print(f"No. {index}")
        messages = []
        history = []
        answer = 0

        record[input_drug] = {}
        record[input_drug]['skip_conversation_round'] = -1
        record[input_drug]['retrieval_conversation'] = [{'result': i} for i in range((args['num_round'] + 1))]

        PDDS_prompt = construct_PDDS_prompt(task_specification_dict, input_drug, args['task'], args['constraint'], args['test_mol'])
        messages.append(PDDS_prompt)

        for round_index in range((args['num_round'] + 1)):
            answer, output_drug, history = conversation(
                clt=client, query_db=local_data, messages=messages, history=history, C=args['num_round'],
                round_index=round_index, task=args['task'], input_drug=input_drug, retrieval_DB=prompt_db,
                record=record, logfile=f, constraint=args['constraint'])

            if answer != 0 or output_drug == None:
                break

        num_all += 1
        if answer == -1:  # invalid smiles
            pass
        elif answer == -2:  # invalid evaluation process
            num_invalid_eval += 1
        elif answer == 1:  # successful optimization
            num_correct += 1
            num_valid += 1
        else:  # failed optimization
            num_valid += 1

        print(f'Acc = {num_correct}/{num_all}', file=f)
        print("----------------", file=f)
        print(f'Acc = {num_correct}/{num_all}')
        print("----------------")

    print("--------Final Acc--------", file=f)
    print(f'Acc = {num_correct}/{num_all-num_invalid_eval}', file=f)
    print("----------------", file=f)

    print("--------Final Validity--------", file=f)
    print(f'Acc = {num_valid}/{num_all-num_invalid_eval}', file=f)
    print("----------------", file=f)

    with open(args['record_file'], 'w') as rf:
        json.dump(record, rf, ensure_ascii=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--task', action='store', required=True, type=str, choices=["qed+", "acceptor+", "donor+", "esol+", "bbbp+", "herg-", "esol+acceptor+", "qed+bbbp+"], help='task name')
    parser.add_argument('--log_file', action='store', required=False, type=str, default='results/DrugAssist.log', help='saved log file name')
    parser.add_argument('--record_file', action='store', required=False, type=str, default='results/DrugAssist.json', help='saved record file name')
    parser.add_argument('--constraint', choices=['loose', 'strict'], required=False, type=str, default='loose', help='Choose between loose and strict mode')
    parser.add_argument('--num_round', required=False, type=int, default=2, help='number of conversation round')
    parser.add_argument('--client_add', required=True, type=str, help='DrugAssist client address set by gradio service')
    parser.add_argument('--database', required=False, type=str, default='mol_DB.csv',help='molecule database')
    parser.add_argument('--test_mol', required=False, type=str, default='testset.csv', help='testset')
    args = parser.parse_args()
    args = vars(args)


    if 'esol' in args['task'] or 'bbbp' in args['task'] or 'herg' in args['task']:
        choice = input("The property value involved in this task will be given by a specific prediction model, make sure you have implemented it yourself!\
                       If you have any questions, please refer to https://github.com/blazerye/DrugAssist/blob/main/evaluate/evaluate.md. Enter y to confirm, enter n to exit\n")
        if choice == 'y':
            main_assit(args)
    else:
        main_assit(args)

