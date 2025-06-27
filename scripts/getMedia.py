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


environments = np.random.dirichlet(np.ones(len(metabolites)))*10


