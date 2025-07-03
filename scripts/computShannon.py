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

def hill_diversity(abundances, q=1):
    abundances = np.asarray(abundances, dtype=np.float64)
    total = abundances.sum()
    if total == 0:
        return 0.0
    p = abundances / total

    if q == 1:
        # Shannon diversity
        p = p[p > 0]  # avoid log(0)
        return np.exp(-np.sum(p * np.log(p)))
    else:
        return np.power(np.sum(p**q), 1 / (1 - q))

def loadDiversity(filePath, ko = False):
    TIME_POINTS = 100
    
    BACTERIA_N = 11
    
    if ko:
        BACTERIA_N = 10
    
    
    with open(filePath) as f:
        data = np.zeros((TIME_POINTS, BACTERIA_N))
        counter = 0
        f.readline()
        for line in f:
            a = line.rstrip('\n').split('\t')
            data[counter] = np.array(a[1:BACTERIA_N + 1]).astype(np.float32)
            counter+=1
    
    #data = data.T 
    

        
    return hill_diversity(data[-1])


def getDiversity(folder, ko= False):
    rootPath = os.path.join(Path(os.getcwd()).parents[0], 'files', folder)

    hill = []

    for i in range(1,16):
        hill.append(loadDiversity(os.path.join(rootPath, str(i) + '.tsv'), ko))        
            
    return hill

all_hill = getDiversity('100_env_all_models')            
akkermansia_hill = getDiversity('akkermansia', 1)
anaerostipes_hill = getDiversity('anaerostipes', 1)
bacteroidesFragilis_hill = getDiversity('BacteroidesFragilis',1)
bacteroidesTheta_hill = getDiversity('BacteroidesTheta', 1)
bifidobacterium_hill = getDiversity('Bifidobacterium', 1)
blautia_hill = getDiversity('Blautia', 1)
clostridum_hill = getDiversity('Clostridum', 1)
faecalibacterium_hill = getDiversity('Faecalibacterium', 1)
lactobacillus_hill = getDiversity('Lactobacillus', 1)
prevotella_hill = getDiversity('Prevotella', 1)
roseburia_hill = getDiversity('Roseburia', 1)


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Create a DataFrame for plotting
df = pd.DataFrame({
    'all': all_hill,
    'Akkermansia': akkermansia_hill,
    'Anaerostipes': anaerostipes_hill,
    'BacteroidesFragilis': bacteroidesFragilis_hill,
    'BacteroidesTheta': bacteroidesTheta_hill,
    'Bifidobacterium' : bifidobacterium_hill,
    'Blautia': blautia_hill,
    'Clostridium': clostridum_hill,
    'Faecalibacterium':faecalibacterium_hill,
    'Lactobacillus': lactobacillus_hill,
    'Prevotella' : prevotella_hill,
    'Roseburia' : roseburia_hill
})

# Melt the DataFrame to long format for seaborn
df_melted = df.melt(var_name='Species', value_name='Bray-Curtis Distance')

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Define custom colors
custom_colors = [
    '#A9A9A9',  # Gray for full
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
plt.savefig(os.path.join(Path(os.getcwd()).parents[0], 'files', 'Figures', 'hill.png'), transparent =True, dpi =600)
plt.show()

