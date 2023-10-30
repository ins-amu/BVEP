import subprocess
import re
import random
import string
import os
import optuna
import xml.etree.ElementTree as ET
import argparse
from absl import logging
from scipy.stats import gmean

def generate_random_name():
    vowels = 'aeiou'
    consonants = ''.join(set(string.ascii_lowercase) - set(vowels))
    name_length = random.randint(4, 10)  # Generate a random length between 4 and 10

    name = ''
    for i in range(name_length):
        if i % 2 == 0:
             name += random.choice(consonants)
        else:
             name += random.choice(vowels)
    return name.capitalize()

def replace_value_in_xml(xml_file_path, nstuck, evalmax, minimum_num_sendSol, reception_threshold, evals_threshold, mult_num_sendSol):
    # Parse the XML file
    tree = ET.parse(xml_file_path)
    root = tree.getroot()

    element = root.find(".//minimum_num_sendSol")
    if element is not None:
        element.text = str(minimum_num_sendSol)
    
    element = root.find(".//n_stuck")
    if element is not None:
        element.text = str(nstuck)

    element = root.find(".//evalmax")
    if element is not None:
        element.text = str(evalmax)

    element = root.find(".//reception_threshold")
    if element is not None:
        element.text = str(reception_threshold)

    element = root.find(".//evals_threshold")
    if element is not None:
        element.text = str(evals_threshold)        

    element = root.find(".//mult_num_sendSol")
    if element is not None:
        element.text = str(mult_num_sendSol)  

    tree.write(xml_file_path)

def opt_func(trial):
#include
    filename, file_extension = os.path.splitext(input_file)
    new_filename = filename + str(generate_random_name())
    new_xml_file_path = new_filename + file_extension
    command = 'cp ' +  input_file + ' ' + new_xml_file_path
    output = subprocess.check_output(command, shell=True)

    nstuck = trial.suggest_categorical('user_options_nstuck', [10 , 20 , 30, 40, 50, 60, 70, 80, 90, 100 ])
    evalmax = trial.suggest_int('local_options_evalmax', 1, 1000) * 85
    evals_threshold = trial.suggest_categorical('cooperation_evals_threshold', list(range(500, 50001, 500)))
    mult_num_sendSol = trial.suggest_int('cooperation_mult_num_sendSol', 5, 20)
    minimum_num_sendSol = trial.suggest_int('cooperation_minimum_num_sendSol', 10, 40)
    reception_threshold = trial.suggest_float('cooperation_threshold_adapt', 1, 20)

    print(new_xml_file_path)
    print("nstuck "+ str(nstuck))
    print("evalmax " +  str(evalmax))
    print("minimum_num_sendSol " + str(minimum_num_sendSol))
    print("reception_threshold " + str(reception_threshold))
    print("evals_threshold " + str(evals_threshold))
    print("mult_num_sendSol "+ str(mult_num_sendSol))

    replace_value_in_xml(new_xml_file_path, nstuck, evalmax, minimum_num_sendSol, reception_threshold, evals_threshold, mult_num_sendSol)

    name_output = 'output/' + 'optuna_' +  str(generate_random_name());

    list_values = []
    for i in range(5):
        print(command)
        command = 'mpirun -np 12 ' +  './bin/paralleltestbed ' + new_xml_file_path + ' ' +  name_output 
        print(command)

        subprocess.check_output(command, shell=True)

        pattern = r"fx = (\d+\.\d+)"
        command = 'cat ' + name_output + "/logfile_id0 | grep 'fx ='" 
        print(command)

        output = subprocess.check_output(command, shell=True)
        match = re.search(pattern, output.decode())

        if match:
            performance_metric = float(match.group(1))
            print(performance_metric)
        else:
            print("No match found.")

        list_values.append(performance_metric)

        command = 'rm -Rf ' +  name_output                            
        subprocess.check_output(command, shell=True)
    command = 'rm ' + new_xml_file_path
    subprocess.check_output(command, shell=True)
    geometric_mean = gmean(list_values)

    return geometric_mean
    
    

# -max-iters
# -evals_th 
# -min_n_sendSol
# -mult_n_sendSol
# -nstuck
#search_space = {
#    'max_iters_ls': ('integer', [1000, 1000000]),
#    'evals_threshold': ('integer', [1000, 100000]),
#    'nstuck': ('integer', [5, 100]),
#    'minimum_num_sendSol': ('integer', [10, 100]),
#    'mult_num_sendSol': ('integer', [5, 100]),
#}

parser = argparse.ArgumentParser(description="Process an XML file.")
parser.add_argument('input_file', type=str, help='The path to the XML file')

args = parser.parse_args()
input_file = args.input_file
print(input_file)

study = optuna.create_study(direction='minimize')
study.optimize(opt_func, n_trials=300)

print("Best trial:")
trial = study.best_trial
print("  Value: ", trial.value)
print("  Params: ")
for key, value in trial.params.items():
    print("    {}: {}".format(key, value))
