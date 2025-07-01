# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 18:08:08 2025

@author: drgarza
"""

import os
from pathlib import Path

import numpy as np

from scipy.spatial.distance import braycurtis

from scipy.spatial.distance import cosine

def loadKinetics(filePath):
    TIME_POINTS = 100
    time = np.zeros(TIME_POINTS)
    
    with open(filePath) as f:
        variables = f.readline().rstrip('\n').split('\t')[1::]
        data = {var:np.zeros(TIME_POINTS) for var in variables}
        counter = 0
        for line in f:
            a = line.rstrip('\n').split('\t')
            for num in range(len(variables)):
                data[variables[num]][counter] = float(a[num + 1])
            
            time[counter] = float(a[0])
            counter+=1
    return time, data


def getTransformedDict(folder, max_v = 16):
    rootPath = os.path.join(Path(os.getcwd()).parents[0], 'files', folder)

    all_kinetics = {}
    counter = 0

    for i in range(1,max_v):
        all_kinetics[counter] = {}
        all_kinetics[counter]['time'],all_kinetics[counter]['data'] = loadKinetics(os.path.join(rootPath, str(i) + '.tsv'))
        counter+=1
        
    all_transformed_dict = {}


    for i in all_kinetics:
        all_transformed_dict[i] = {}
        for z in all_kinetics[i]['data']:
            if z != 'Oxygen':
                all_transformed_dict[i][z] = all_kinetics[i]['data'][z]/np.max(all_kinetics[i]['data'][z])
            
    return all_transformed_dict



def compareTransformedDict(dict_1, dict_2, threshold = 0.1):
    'always make the leave one out the second one'
    v1, v2 = [],[]
    
    used = []
    for i in dict_2:
        trace = dict_2[i]
        if (trace.max() - trace.min()) < threshold:
            trace2 = dict_1[i]
            if (trace2.max() - trace2.min()) < threshold:
                print(i)
                pass
            else:
                used.append(i)
        else:
            used.append(i)
    
    for i in used:
        v1.append(dict_1[i][-1])
        v2.append(dict_2[i][-1])
    return braycurtis(v1,v2)
            

all_transformed_dict = getTransformedDict('100_env_all_models')

akkermansia_dict = getTransformedDict('akkermansia')

anaerostipes_dict = getTransformedDict('anaerostipes')

faecalibacterium_dict = getTransformedDict('Faecalibacterium')

lactobacillus_dict = getTransformedDict('Lactobacillus')


blautia_dict = getTransformedDict('Blautia', max_v = 10)

akkermansia_bc = [compareTransformedDict(all_transformed_dict[i], akkermansia_dict[i]) for i in range(15)]

anaerostipes_bc = [compareTransformedDict(all_transformed_dict[i], anaerostipes_dict[i]) for i in range(15)]

faecalibacterium_bc = [compareTransformedDict(all_transformed_dict[i], faecalibacterium_dict[i]) for i in range(15)]

lactobacillus_bc = [compareTransformedDict(all_transformed_dict[i], lactobacillus_dict[i]) for i in range(15)]

blautia_bc = [compareTransformedDict(all_transformed_dict[i], blautia_dict[i]) for i in range(9)]


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Create a DataFrame for plotting
df = pd.DataFrame({
    'Akkermansia': akkermansia_bc,
    'Anaerostipes': anaerostipes_bc,
    'Blautia': blautia_bc + [np.nan]*(15-9),  # pad with NaN to align length
    'Faecalibacterium':faecalibacterium_bc,
    'Lactobacillus': lactobacillus_bc
})

# Melt the DataFrame to long format for seaborn
df_melted = df.melt(var_name='Species', value_name='Bray-Curtis Distance')

# Make the violin plot
plt.figure(figsize=(6, 4))
sns.boxplot(data=df_melted, x='Species', y='Bray-Curtis Distance')
sns.swarmplot(data=df_melted, x='Species', y='Bray-Curtis Distance', color='k', alpha=0.5, size=3)

plt.title('Bray-Curtis Dissimilarities by Species Removed')
plt.tight_layout()
plt.show()
