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

from aquarel import load_theme


theme = load_theme("boxy_light")
theme.apply()


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



def compareTransformedDict(dict_1, dict_2, threshold = 0.05):
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

bacteroidesfragilis_dict = getTransformedDict('BacteroidesFragilis')

bacteroidestheta_dict = getTransformedDict('BacteroidesTheta')

bifidobacterium_dict = getTransformedDict('Bifidobacterium')

blautia_dict = getTransformedDict('Blautia')

clostridum_dict = getTransformedDict('Clostridum')

faecalibacterium_dict = getTransformedDict('Faecalibacterium')

lactobacillus_dict = getTransformedDict('Lactobacillus')

prevotella_dict = getTransformedDict('Prevotella' )

roseburia_dict = getTransformedDict('Roseburia')




akkermansia_bc = [compareTransformedDict(all_transformed_dict[i], akkermansia_dict[i]) for i in range(15)]


anaerostipes_bc = [compareTransformedDict(all_transformed_dict[i], anaerostipes_dict[i]) for i in range(15)]

bacteroidesfragilis_bc = [compareTransformedDict(all_transformed_dict[i], bacteroidesfragilis_dict[i]) for i in range(15)]

bacteroidestheta_bc = [compareTransformedDict(all_transformed_dict[i], bacteroidestheta_dict[i]) for i in range(15)]

bifidobacterium_bc = [compareTransformedDict(all_transformed_dict[i], bifidobacterium_dict[i]) for i in range(15)]

blautia_bc = [compareTransformedDict(all_transformed_dict[i], blautia_dict[i]) for i in range(15)]

clostridum_bc = [compareTransformedDict(all_transformed_dict[i], clostridum_dict[i]) for i in range(15)]

faecalibacterium_bc = [compareTransformedDict(all_transformed_dict[i], faecalibacterium_dict[i]) for i in range(15)]

lactobacillus_bc = [compareTransformedDict(all_transformed_dict[i], lactobacillus_dict[i]) for i in range(15)]

prevotella_bc = [compareTransformedDict(all_transformed_dict[i], prevotella_dict[i]) for i in range(15)]

roseburia_bc = [compareTransformedDict(all_transformed_dict[i], roseburia_dict[i]) for i in range(15)]



import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Create a DataFrame for plotting
df = pd.DataFrame({
    'Akkermansia': akkermansia_bc,
    'Anaerostipes': anaerostipes_bc,
    'BacteroidesFragilis': bacteroidesfragilis_bc,
    'BacteroidesTheta': bacteroidestheta_bc,
    'Bifidobacterium' : bifidobacterium_bc,
    'Blautia': blautia_bc,
    'Clostridium': clostridum_bc,
    'Faecalibacterium':faecalibacterium_bc,
    'Lactobacillus': lactobacillus_bc,
    'Prevotella' : prevotella_bc,
    'Roseburia' : roseburia_bc
})

# Melt the DataFrame to long format for seaborn
df_melted = df.melt(var_name='Species', value_name='Bray-Curtis Distance')

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Define custom colors
custom_colors = [
    '#000000',  # Black
    '#E69F00',  # Orange
    '#56B4E9',  # Sky Blue
    '#009E73',  # Bluish Green
    '#F0E442',  # Yellow
    '#0072B2',  # Blue
    '#D55E00',  # Vermilion
    '#CC79A7',  # Reddish Purple
    '#8B4513',  # Saddle Brown
    '#32CD32',  # Lime Green
    '#FF1493'   # Deep Pink
]

plt.figure(figsize=(12, 4.5))

ax = sns.boxplot(
    data=df_melted, x='Species', y='Bray-Curtis Distance',
    showcaps=True,
    boxprops=dict(facecolor='white', edgecolor='black', linewidth=2),
    whiskerprops=dict(linewidth=2),
    capprops=dict(linewidth=2),
    medianprops=dict(color='black', linewidth=2),
    flierprops=dict(marker='o', markersize=5, linestyle='none')
)

# Apply custom box facecolors
for i, patch in enumerate(ax.patches[:len(custom_colors)]):
    patch.set_facecolor(custom_colors[i])
    #patch.set_alpha(0)
    patch.set_edgecolor('black')
    patch.set_linewidth(2)

# Overlay swarm points
sns.swarmplot(
    data=df_melted, x='Species', y='Bray-Curtis Distance',
    color='k', size=5
)

# Hide x-axis labels
ax.set_xticklabels([''] * len(df.columns), rotation=45)
ax.tick_params(axis='y', labelsize=16)

# Remove axis labels
ax.set_xlabel('')
ax.set_ylabel('')

# Optional: also remove axis label padding (if needed)
ax.xaxis.labelpad = 0
ax.yaxis.labelpad = 0

#plt.title('Bray-Curtis Dissimilarities by Species Removed', fontsize=12)
plt.tight_layout()
plt.grid(axis='y', linestyle='--', linewidth=0.5, alpha=0.6)
plt.savefig(os.path.join(Path(os.getcwd()).parents[0], 'files', 'Figures', 'braycurtis.png'), transparent =True, dpi =600)
plt.show()

