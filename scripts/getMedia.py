# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 10:36:14 2025

@author: drgarza
"""

import os
from pathlib import Path
os.environ['GRB_LICENSE_FILE'] = os.path.join(Path(os.getcwd()), 'licenses', 'gurobi.lic')

import csv
import numpy as np

import cobra
#cobra.Configuration().solver = "glpk"

from dfba import *

model_folder = os.path.join(Path(os.getcwd()).parents[0], 'AgoraModels')




NUM_SAMPLES = 100
WATER_V = 10
OXYGEN_V = 0


def get_exchange(modelPath):
    model = cobra.io.read_sbml_model(modelPath)

    exchanges_ids = [model.reactions.get_by_id(i.id).id for i in model.reactions if 'EX_' in i.id]
    
    fva = cobra.flux_analysis.flux_variability_analysis(model, reaction_list = exchanges_ids, processes=1)
    selected =  list(fva[fva['minimum']<0].index)
    
    return {i: model.reactions.get_by_id(i).reactants[0].name for i in selected}
    

exchanges = {}

for i in os.listdir(model_folder):
    exchanges = exchanges | get_exchange(os.path.join(model_folder, i))
    
exchanges.pop('EX_h2o(e)')
exchanges.pop('EX_o2(e)')

metabolites = list(exchanges.keys())


environments = np.random.dirichlet(np.ones(len(metabolites)), size=NUM_SAMPLES) * 1000


# Append water and oxygen to metabolites and their names
full_ids = metabolites + ['EX_h2o(e)', 'EX_o2(e)']
full_names = [exchanges[mid] for mid in metabolites] + ['Water', 'Oxygen']

# Add water and oxygen to each sample
media_with_fixed = [list(env) + [WATER_V, OXYGEN_V] for env in environments]

# Write TSV
output_path = os.path.join(os.getcwd(), "media.tsv")
with open(output_path, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')

    # First row: metabolite IDs
    writer.writerow(["sample"] + full_ids)

    # Second row: metabolite names
    writer.writerow([""] + full_names)

    # Data rows: media samples
    for i, row in enumerate(media_with_fixed):
        writer.writerow([f"sample_{i+1}"] + row)

print(f"Wrote {NUM_SAMPLES} media compositions to {output_path}")

