# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 19:50:15 2025

@author: drgarza
"""

import os, warnings
from pathlib import Path
os.environ['GRB_LICENSE_FILE'] = os.path.join(Path(os.getcwd()), 'licenses', 'gurobi.lic')

import csv
import numpy as np

import cobra
#cobra.Configuration().solver = "glpk"
warnings.filterwarnings("ignore", category=UserWarning)

from dfba import *
import gc

def load_media_sample(sample_id: int, media_file="media.tsv"):
    """
    Load a media sample and return:
    - env_dict: {exchange_id: concentration}
    - env_names: list of metabolite names (in order)
    - env_ids: list of exchange reaction IDs (in order)
    """
    with open(media_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # First line: IDs (e.g., EX_glc(e), EX_o2(e), ...)
    ids = lines[0].rstrip('\n').split('\t')[1:]  # skip 'sample' column

    # Second line: metabolite names, must match same number of columns
    name_line = lines[1].rstrip('\n').split('\t')[1:]
    if len(name_line) != len(ids):
        raise ValueError(f"Mismatch: {len(name_line)} names but {len(ids)} IDs")

    env_names = name_line
    line_idx = sample_id + 1
    if line_idx >= len(lines):
        raise IndexError(f"Sample {sample_id} out of range. Only {len(lines)-2} samples in file.")

    # Sample row
    values = lines[line_idx].rstrip('\n').split('\t')[1:]  # skip sample name
    concentrations = list(map(float, values))

    if len(concentrations) != len(ids):
        raise ValueError(f"Mismatch in data row: {len(concentrations)} values but {len(ids)} IDs")

    env_dict = dict(zip(ids, concentrations))
    return env_dict, env_names, ids


def store_simulation_to_tsv(filepath, 
                            sol, 
                            model_names, 
                            metabolite_names):
    """
    Store simulation results in a TSV file.

    Parameters:
    - filepath: str, path to the output .tsv file
    - sol: object returned by solve_ivp (must have .t and .y)
    - model_names: list of strings, names for each biomass trajectory
    - metabolite_names: list of strings, names for each metabolite
    """
    n_models = len(model_names)
    n_mets = len(metabolite_names)
    timepoints = sol.t

    # Sanity check
    if sol.y.shape[0] != n_models + n_mets:
        raise ValueError("Mismatch: sol.y rows do not match model + metabolite counts.")

    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')

        # Write header
        header = ['time'] + model_names + metabolite_names
        writer.writerow(header)

        # Write data
        for i, t in enumerate(timepoints):
            row = [t]
            # Biomasses
            row += [sol.y[j, i] for j in range(n_models)]
            # Metabolites
            row += [sol.y[n_models + k, i] for k in range(n_mets)]
            writer.writerow(row)

    print(f"Saved simulation results to: {filepath}")



model_folder = os.path.join(Path(os.getcwd()).parents[0], 'AgoraModels', 'Akkermansia')


for env_num in range(1,26):    
    env, env_names, env_ids = load_media_sample(env_num)
    
    
    
    models = [cobra.io.read_sbml_model(os.path.join(model_folder, i)) for i in os.listdir(model_folder) if '.xml' in i]
    
    biomassNames = []
    
    for mod in models:
        for react in mod.reactions:
            if react.objective_coefficient == 1.0:
                biomassNames.append(react.id)
    np.random.seed(666)
    monod_c = np.random.lognormal(1, size = len(env))          
    sol = simulate(models, 
                   env_ids,
                   env_names, 
                   env, 
                   monod_c=monod_c,
                   show_plots=True, 
                   biomassNames=biomassNames)
    
    model_names = [i.id.split('_')[1] + '_' + i.id.split('_')[2] for i in models]
    
    filePath = os.path.join(Path(os.getcwd()).parents[0], 'files', 'akkermansia', str(env_num)+'.tsv')
    store_simulation_to_tsv(filePath, 
                                sol, 
                                model_names, 
                                env_names)
    gc.collect()








model_names = []