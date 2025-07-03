# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 11:10:05 2025

@author: drgarza
"""

import os
from pathlib import Path

import numpy as np

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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


def getAcid(folder, acid ='acetate', max_v = 16):
    rootPath = os.path.join(Path(os.getcwd()).parents[0], 'files', folder)

    acid_c = []

    for i in range(1,max_v):
        _, kinetics = loadKinetics(os.path.join(rootPath, str(i) + '.tsv'))
        acid_c.append(kinetics[acid][-1])        
            
    return acid_c




def plot_acid_concentrations(acid_name, all_values, knockout_dict, custom_colors, filePath = None):
    """
    Plots final concentration of a given acid across all knockouts vs full community.

    Parameters:
        acid_name (str): e.g., 'acetate'
        all_values (list): values for the full community
        knockout_dict (dict): {species_name: list of values}
        custom_colors (list): list of hex colors for knockout boxes (gray added for full community)
    """
    # Combine into a DataFrame
    df = pd.DataFrame(knockout_dict)
    
    # Add full community as the first column
    df.insert(0, 'Full', all_values + [np.nan] * (df.shape[0] - len(all_values)))

    # Melt for seaborn
    df_melted = df.melt(var_name='Condition', value_name=f'{acid_name} Concentration')

    # Define plot
    plt.figure(figsize=(12, 4.5))
    ax = sns.boxplot(
        data=df_melted, x='Condition', y=f'{acid_name} Concentration',
        showcaps=True,
        boxprops=dict(facecolor='white', edgecolor='black', linewidth=2),
        whiskerprops=dict(linewidth=2),
        capprops=dict(linewidth=2),
        medianprops=dict(color='black', linewidth=2),
        flierprops=dict(marker='o', markersize=5, linestyle='none')
    )

    # Color each box
    for i, patch in enumerate(ax.patches[:len(custom_colors)]):
        patch.set_facecolor(custom_colors[i])
        patch.set_alpha(0.7)
        patch.set_edgecolor('black')
        patch.set_linewidth(2)

    # Swarmplot
    sns.swarmplot(
        data=df_melted, x='Condition', y=f'{acid_name} Concentration',
        color='black', size=5, alpha=0.7
    )

    # Hide x-axis labels for later
    ax.set_xticklabels([''] * len(df.columns))

    # Style
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(axis='y', labelsize=12)
    plt.title(f'{acid_name.capitalize()} Concentration Across Knockouts', fontsize=12)
    plt.grid(axis='y', linestyle='--', linewidth=0.5, alpha=0.6)
    plt.tight_layout()
    if filePath is not None:
        plt.savefig(filePath, transparent =True, dpi =600)
    plt.show()




def getCombinedAcid(folder, acid1='(r)-lactate', acid2='(s)-lactate', max_v=16):
    rootPath = os.path.join(Path(os.getcwd()).parents[0], 'files', folder)
    combined = []

    for i in range(1, max_v):
        _, kinetics = loadKinetics(os.path.join(rootPath, str(i) + '.tsv'))
        combined_value = kinetics[acid1][-1] + kinetics[acid2][-1]
        combined.append(combined_value)

    return combined




#acetate#############

ACID = 'acetate'
all_acetate = getAcid('100_env_all_models', acid = ACID)
akkermansia_acetate = getAcid('akkermansia', acid = ACID)
anaerostipes_acetate = getAcid('anaerostipes', acid = ACID)
bacteroidesFragilis_acetate = getAcid('BacteroidesFragilis', acid = ACID)
bacteroidesTheta_acetate = getAcid('BacteroidesTheta', acid = ACID)
bifidobacterium_acetate = getAcid('Bifidobacterium', acid = ACID)
blautia_acetate = getAcid('Blautia', acid = ACID)
clostridum_acetate = getAcid('Clostridum', acid = ACID)
faecalibacterium_acetate = getAcid('Faecalibacterium', acid = ACID)
lactobacillus_acetate = getAcid('Lactobacillus', acid = ACID)
prevotella_acetate = getAcid('Prevotella', acid = ACID)
roseburia_acetate = getAcid('Roseburia', acid = ACID)

acid = 'acetate'

knockouts = {
    'Akkermansia': akkermansia_acetate,
    'Anaerostipes': anaerostipes_acetate,
    'BacteroidesFragilis': bacteroidesFragilis_acetate,
    'BacteroidesTheta': bacteroidesTheta_acetate,
    'Bifidobacterium': bifidobacterium_acetate,
    'Blautia': blautia_acetate,
    'Clostridium': clostridum_acetate,
    'Faecalibacterium': faecalibacterium_acetate,
    'Lactobacillus': lactobacillus_acetate,
    'Prevotella': prevotella_acetate ,
    'Roseburia': roseburia_acetate ,
}

colors = [
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

plot_acid_concentrations(acid_name=ACID, all_values=all_acetate, knockout_dict=knockouts, custom_colors=colors,
                         filePath = os.path.join(Path(os.getcwd()).parents[0], 'files', 'Figures', ACID + '.png'))


#butyrate#############

ACID = 'butyrate'
all_butyrate = getAcid('100_env_all_models', acid = ACID)
akkermansia_butyrate = getAcid('akkermansia', acid = ACID)
anaerostipes_butyrate = getAcid('anaerostipes', acid = ACID)
bacteroidesFragilis_butyrate = getAcid('BacteroidesFragilis', acid = ACID)
bacteroidesTheta_butyrate = getAcid('BacteroidesTheta', acid = ACID)
bifidobacterium_butyrate = getAcid('Bifidobacterium', acid = ACID)
blautia_butyrate = getAcid('Blautia', acid = ACID)
clostridum_butyrate = getAcid('Clostridum', acid = ACID)
faecalibacterium_butyrate = getAcid('Faecalibacterium', acid = ACID)
lactobacillus_butyrate = getAcid('Lactobacillus', acid = ACID)
prevotella_butyrate = getAcid('Prevotella', acid = ACID)
roseburia_butyrate = getAcid('Roseburia', acid = ACID)

acid = 'butyrate'

knockouts = {
    'Akkermansia': akkermansia_butyrate,
    'Anaerostipes': anaerostipes_butyrate,
    'BacteroidesFragilis': bacteroidesFragilis_butyrate,
    'BacteroidesTheta': bacteroidesTheta_butyrate,
    'Bifidobacterium': bifidobacterium_butyrate,
    'Blautia': blautia_butyrate,
    'Clostridium': clostridum_butyrate,
    'Faecalibacterium': faecalibacterium_butyrate,
    'Lactobacillus': lactobacillus_butyrate,
    'Prevotella': prevotella_butyrate ,
    'Roseburia': roseburia_butyrate ,
}

colors = [
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

plot_acid_concentrations(acid_name=ACID, all_values=all_butyrate, knockout_dict=knockouts, custom_colors=colors,
                         filePath = os.path.join(Path(os.getcwd()).parents[0], 'files', 'Figures', ACID + '.png'))


#formate#############

ACID = 'formate'
all_formate = getAcid('100_env_all_models', acid = ACID)
akkermansia_formate = getAcid('akkermansia', acid = ACID)
anaerostipes_formate = getAcid('anaerostipes', acid = ACID)
bacteroidesFragilis_formate = getAcid('BacteroidesFragilis', acid = ACID)
bacteroidesTheta_formate = getAcid('BacteroidesTheta', acid = ACID)
bifidobacterium_formate = getAcid('Bifidobacterium', acid = ACID)
blautia_formate = getAcid('Blautia', acid = ACID)
clostridum_formate = getAcid('Clostridum', acid = ACID)
faecalibacterium_formate = getAcid('Faecalibacterium', acid = ACID)
lactobacillus_formate = getAcid('Lactobacillus', acid = ACID)
prevotella_formate = getAcid('Prevotella', acid = ACID)
roseburia_formate = getAcid('Roseburia', acid = ACID)

acid = 'formate'

knockouts = {
    'Akkermansia': akkermansia_formate,
    'Anaerostipes': anaerostipes_formate,
    'BacteroidesFragilis': bacteroidesFragilis_formate,
    'BacteroidesTheta': bacteroidesTheta_formate,
    'Bifidobacterium': bifidobacterium_formate,
    'Blautia': blautia_formate,
    'Clostridium': clostridum_formate,
    'Faecalibacterium': faecalibacterium_formate,
    'Lactobacillus': lactobacillus_formate,
    'Prevotella': prevotella_formate ,
    'Roseburia': roseburia_formate ,
}

colors = [
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

plot_acid_concentrations(acid_name=ACID, all_values=all_formate, knockout_dict=knockouts, custom_colors=colors,
                         filePath = os.path.join(Path(os.getcwd()).parents[0], 'files', 'Figures', ACID + '.png'))

#fumarate#############

ACID = 'fumarate'
all_fumarate = getAcid('100_env_all_models', acid = ACID)
akkermansia_fumarate = getAcid('akkermansia', acid = ACID)
anaerostipes_fumarate = getAcid('anaerostipes', acid = ACID)
bacteroidesFragilis_fumarate = getAcid('BacteroidesFragilis', acid = ACID)
bacteroidesTheta_fumarate = getAcid('BacteroidesTheta', acid = ACID)
bifidobacterium_fumarate = getAcid('Bifidobacterium', acid = ACID)
blautia_fumarate = getAcid('Blautia', acid = ACID)
clostridum_fumarate = getAcid('Clostridum', acid = ACID)
faecalibacterium_fumarate = getAcid('Faecalibacterium', acid = ACID)
lactobacillus_fumarate = getAcid('Lactobacillus', acid = ACID)
prevotella_fumarate = getAcid('Prevotella', acid = ACID)
roseburia_fumarate = getAcid('Roseburia', acid = ACID)

acid = 'fumarate'

knockouts = {
    'Akkermansia': akkermansia_fumarate,
    'Anaerostipes': anaerostipes_fumarate,
    'BacteroidesFragilis': bacteroidesFragilis_fumarate,
    'BacteroidesTheta': bacteroidesTheta_fumarate,
    'Bifidobacterium': bifidobacterium_fumarate,
    'Blautia': blautia_fumarate,
    'Clostridium': clostridum_fumarate,
    'Faecalibacterium': faecalibacterium_fumarate,
    'Lactobacillus': lactobacillus_fumarate,
    'Prevotella': prevotella_fumarate ,
    'Roseburia': roseburia_fumarate ,
}

colors = [
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

plot_acid_concentrations(acid_name=ACID, all_values=all_fumarate, knockout_dict=knockouts, custom_colors=colors,
                         filePath = os.path.join(Path(os.getcwd()).parents[0], 'files', 'Figures', ACID + '.png'))

#propionate#############

ACID = 'propionate'
all_propionate = getAcid('100_env_all_models', acid = ACID)
akkermansia_propionate = getAcid('akkermansia', acid = ACID)
anaerostipes_propionate = getAcid('anaerostipes', acid = ACID)
bacteroidesFragilis_propionate = getAcid('BacteroidesFragilis', acid = ACID)
bacteroidesTheta_propionate = getAcid('BacteroidesTheta', acid = ACID)
bifidobacterium_propionate = getAcid('Bifidobacterium', acid = ACID)
blautia_propionate = getAcid('Blautia', acid = ACID)
clostridum_propionate = getAcid('Clostridum', acid = ACID)
faecalibacterium_propionate = getAcid('Faecalibacterium', acid = ACID)
lactobacillus_propionate = getAcid('Lactobacillus', acid = ACID)
prevotella_propionate = getAcid('Prevotella', acid = ACID)
roseburia_propionate = getAcid('Roseburia', acid = ACID)

acid = 'propionate'

knockouts = {
    'Akkermansia': akkermansia_propionate,
    'Anaerostipes': anaerostipes_propionate,
    'BacteroidesFragilis': bacteroidesFragilis_propionate,
    'BacteroidesTheta': bacteroidesTheta_propionate,
    'Bifidobacterium': bifidobacterium_propionate,
    'Blautia': blautia_propionate,
    'Clostridium': clostridum_propionate,
    'Faecalibacterium': faecalibacterium_propionate,
    'Lactobacillus': lactobacillus_propionate,
    'Prevotella': prevotella_propionate ,
    'Roseburia': roseburia_propionate ,
}

colors = [
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

plot_acid_concentrations(acid_name=ACID, all_values=all_propionate, knockout_dict=knockouts, custom_colors=colors,
                         filePath = os.path.join(Path(os.getcwd()).parents[0], 'files', 'Figures', ACID + '.png'))

#######################
acid = 'lactate_total'

# Full community
all_lactate = getCombinedAcid('100_env_all_models')

# Knockouts
knockouts = {
    'Akkermansia': [r + s for r, s in zip(
        getAcid('akkermansia', acid='(r)-lactate'),
        getAcid('akkermansia', acid='(s)-lactate')
    )],
    'Anaerostipes': [r + s for r, s in zip(
        getAcid('anaerostipes', acid='(r)-lactate'),
        getAcid('anaerostipes', acid='(s)-lactate')
    )],
    'BacteroidesFragilis': [r + s for r, s in zip(
        getAcid('BacteroidesFragilis', acid='(r)-lactate'),
        getAcid('BacteroidesFragilis', acid='(s)-lactate')
    )],
    'BacteroidesTheta': [r + s for r, s in zip(
        getAcid('BacteroidesTheta', acid='(r)-lactate'),
        getAcid('BacteroidesTheta', acid='(s)-lactate')
    )],
    'Bifidobacterium': [r + s for r, s in zip(
        getAcid('Bifidobacterium', acid='(r)-lactate'),
        getAcid('Bifidobacterium', acid='(s)-lactate')
    )],
    'Blautia': [r + s for r, s in zip(
        getAcid('Blautia', acid='(r)-lactate'),
        getAcid('Blautia', acid='(s)-lactate')
    )],
    'Clostridium': [r + s for r, s in zip(
        getAcid('Clostridum', acid='(r)-lactate'),
        getAcid('Clostridum', acid='(s)-lactate')
    )],
    'Faecalibacterium': [r + s for r, s in zip(
        getAcid('Faecalibacterium', acid='(r)-lactate'),
        getAcid('Faecalibacterium', acid='(s)-lactate')
    )],
    'Lactobacillus': [r + s for r, s in zip(
        getAcid('Lactobacillus', acid='(r)-lactate'),
        getAcid('Lactobacillus', acid='(s)-lactate')
    )],
    'Prevotella': [r + s for r, s in zip(
        getAcid('Prevotella', acid='(r)-lactate'),
        getAcid('Prevotella', acid='(s)-lactate')
    )] ,
    'Roseburia': [r + s for r, s in zip(
        getAcid('Roseburia', acid='(r)-lactate'),
        getAcid('Roseburia', acid='(s)-lactate')
    )] ,
}

# Same colors as before
colors = [
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

# Plot it
<<<<<<< HEAD
plot_acid_concentrations(acid_name=acid, all_values=all_lactate, knockout_dict=knockouts, custom_colors=colors)


#co2#############

ACID = 'co2+'
all_co2 = getAcid('100_env_all_models', acid = ACID)
akkermansia_co2 = getAcid('akkermansia', acid = ACID)
anaerostipes_co2 = getAcid('anaerostipes', acid = ACID)
bacteroidesFragilis_co2 = getAcid('BacteroidesFragilis', acid = ACID)
bacteroidesTheta_co2 = getAcid('BacteroidesTheta', acid = ACID)
bifidobacterium_co2 = getAcid('Bifidobacterium', acid = ACID)
blautia_co2 = getAcid('Blautia', acid = ACID)
clostridum_co2 = getAcid('Clostridum', acid = ACID)
faecalibacterium_co2 = getAcid('Faecalibacterium', acid = ACID)
lactobacillus_co2 = getAcid('Lactobacillus', acid = ACID)
prevotella_co2 = getAcid('Prevotella', acid = ACID)
roseburia_co2 = getAcid('Roseburia', acid = ACID)

acid = 'co2'

knockouts = {
    'Akkermansia': akkermansia_co2,
    'Anaerostipes': anaerostipes_co2,
    'BacteroidesFragilis': bacteroidesFragilis_co2,
    'BacteroidesTheta': bacteroidesTheta_co2,
    'Bifidobacterium': bifidobacterium_co2,
    'Blautia': blautia_co2,
    'Clostridium': clostridum_co2,
    'Faecalibacterium': faecalibacterium_co2,
    'Lactobacillus': lactobacillus_co2,
    'Prevotella': prevotella_co2 ,
    'Roseburia': roseburia_co2 ,
}

colors = [
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

plot_acid_concentrations(acid_name=ACID, all_values=all_co2, knockout_dict=knockouts, custom_colors=colors)

#carbon_dioxide#############

ACID = 'carbon dioxide'
all_carbon_dioxide = getAcid('100_env_all_models', acid = ACID)
akkermansia_carbon_dioxide = getAcid('akkermansia', acid = ACID)
anaerostipes_carbon_dioxide = getAcid('anaerostipes', acid = ACID)
bacteroidesFragilis_carbon_dioxide = getAcid('BacteroidesFragilis', acid = ACID)
bacteroidesTheta_carbon_dioxide = getAcid('BacteroidesTheta', acid = ACID)
bifidobacterium_carbon_dioxide = getAcid('Bifidobacterium', acid = ACID)
blautia_carbon_dioxide = getAcid('Blautia', acid = ACID)
clostridum_carbon_dioxide = getAcid('Clostridum', acid = ACID)
faecalibacterium_carbon_dioxide = getAcid('Faecalibacterium', acid = ACID)
lactobacillus_carbon_dioxide = getAcid('Lactobacillus', acid = ACID)
prevotella_carbon_dioxide = getAcid('Prevotella', acid = ACID)
roseburia_carbon_dioxide = getAcid('Roseburia', acid = ACID)

acid = 'carbon_dioxide'

knockouts = {
    'Akkermansia': akkermansia_carbon_dioxide,
    'Anaerostipes': anaerostipes_carbon_dioxide,
    'BacteroidesFragilis': bacteroidesFragilis_carbon_dioxide,
    'BacteroidesTheta': bacteroidesTheta_carbon_dioxide,
    'Bifidobacterium': bifidobacterium_carbon_dioxide,
    'Blautia': blautia_carbon_dioxide,
    'Clostridium': clostridum_carbon_dioxide,
    'Faecalibacterium': faecalibacterium_carbon_dioxide,
    'Lactobacillus': lactobacillus_carbon_dioxide,
    'Prevotella': prevotella_carbon_dioxide ,
    'Roseburia': roseburia_carbon_dioxide ,
}

colors = [
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

plot_acid_concentrations(acid_name=ACID, all_values=all_carbon_dioxide, knockout_dict=knockouts, custom_colors=colors)

#proton#############

ACID = 'proton'
all_proton = getAcid('100_env_all_models', acid = ACID)
akkermansia_proton = getAcid('akkermansia', acid = ACID)
anaerostipes_proton = getAcid('anaerostipes', acid = ACID)
bacteroidesFragilis_proton = getAcid('BacteroidesFragilis', acid = ACID)
bacteroidesTheta_proton = getAcid('BacteroidesTheta', acid = ACID)
bifidobacterium_proton = getAcid('Bifidobacterium', acid = ACID)
blautia_proton = getAcid('Blautia', acid = ACID)
clostridum_proton = getAcid('Clostridum', acid = ACID)
faecalibacterium_proton = getAcid('Faecalibacterium', acid = ACID)
lactobacillus_proton = getAcid('Lactobacillus', acid = ACID)
prevotella_proton = getAcid('Prevotella', acid = ACID)
roseburia_proton = getAcid('Roseburia', acid = ACID)

acid = 'proton'

knockouts = {
    'Akkermansia': akkermansia_proton,
    'Anaerostipes': anaerostipes_proton,
    'BacteroidesFragilis': bacteroidesFragilis_proton,
    'BacteroidesTheta': bacteroidesTheta_proton,
    'Bifidobacterium': bifidobacterium_proton,
    'Blautia': blautia_proton,
    'Clostridium': clostridum_proton,
    'Faecalibacterium': faecalibacterium_proton,
    'Lactobacillus': lactobacillus_proton,
    'Prevotella': prevotella_proton ,
    'Roseburia': roseburia_proton ,
}

colors = [
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

plot_acid_concentrations(acid_name=ACID, all_values=all_proton, knockout_dict=knockouts, custom_colors=colors)





#hydrogen#############

ACID = 'hydrogen'
all_hydrogen = getAcid('100_env_all_models', acid = ACID)
akkermansia_hydrogen = getAcid('akkermansia', acid = ACID)
anaerostipes_hydrogen = getAcid('anaerostipes', acid = ACID)
bacteroidesFragilis_hydrogen = getAcid('BacteroidesFragilis', acid = ACID)
bacteroidesTheta_hydrogen = getAcid('BacteroidesTheta', acid = ACID)
bifidobacterium_hydrogen = getAcid('Bifidobacterium', acid = ACID)
blautia_hydrogen = getAcid('Blautia', acid = ACID)
clostridum_hydrogen = getAcid('Clostridum', acid = ACID)
faecalibacterium_hydrogen = getAcid('Faecalibacterium', acid = ACID)
lactobacillus_hydrogen = getAcid('Lactobacillus', acid = ACID)
prevotella_hydrogen = getAcid('Prevotella', acid = ACID)
roseburia_hydrogen = getAcid('Roseburia', acid = ACID)

acid = 'hydrogen'

knockouts = {
    'Akkermansia': akkermansia_hydrogen,
    'Anaerostipes': anaerostipes_hydrogen,
    'BacteroidesFragilis': bacteroidesFragilis_hydrogen,
    'BacteroidesTheta': bacteroidesTheta_hydrogen,
    'Bifidobacterium': bifidobacterium_hydrogen,
    'Blautia': blautia_hydrogen,
    'Clostridium': clostridum_hydrogen,
    'Faecalibacterium': faecalibacterium_hydrogen,
    'Lactobacillus': lactobacillus_hydrogen,
    'Prevotella': prevotella_hydrogen ,
    'Roseburia': roseburia_hydrogen ,
}

colors = [
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

plot_acid_concentrations(acid_name=ACID, all_values=all_hydrogen, knockout_dict=knockouts, custom_colors=colors)

#ammonium#############

ACID = 'ammonium'
all_ammonium = getAcid('100_env_all_models', acid = ACID)
akkermansia_ammonium = getAcid('akkermansia', acid = ACID)
anaerostipes_ammonium = getAcid('anaerostipes', acid = ACID)
bacteroidesFragilis_ammonium = getAcid('BacteroidesFragilis', acid = ACID)
bacteroidesTheta_ammonium = getAcid('BacteroidesTheta', acid = ACID)
bifidobacterium_ammonium = getAcid('Bifidobacterium', acid = ACID)
blautia_ammonium = getAcid('Blautia', acid = ACID)
clostridum_ammonium = getAcid('Clostridum', acid = ACID)
faecalibacterium_ammonium = getAcid('Faecalibacterium', acid = ACID)
lactobacillus_ammonium = getAcid('Lactobacillus', acid = ACID)
prevotella_ammonium = getAcid('Prevotella', acid = ACID)
roseburia_ammonium = getAcid('Roseburia', acid = ACID)

acid = 'ammonium'

knockouts = {
    'Akkermansia': akkermansia_ammonium,
    'Anaerostipes': anaerostipes_ammonium,
    'BacteroidesFragilis': bacteroidesFragilis_ammonium,
    'BacteroidesTheta': bacteroidesTheta_ammonium,
    'Bifidobacterium': bifidobacterium_ammonium,
    'Blautia': blautia_ammonium,
    'Clostridium': clostridum_ammonium,
    'Faecalibacterium': faecalibacterium_ammonium,
    'Lactobacillus': lactobacillus_ammonium,
    'Prevotella': prevotella_ammonium ,
    'Roseburia': roseburia_ammonium ,
}

colors = [
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

plot_acid_concentrations(acid_name=ACID, all_values=all_ammonium, knockout_dict=knockouts, custom_colors=colors)

=======
plot_acid_concentrations(acid_name=ACID, all_values=all_lactate, knockout_dict=knockouts, custom_colors=colors,
                         filePath = os.path.join(Path(os.getcwd()).parents[0], 'files', 'Figures', ACID + '.png'))
>>>>>>> 7c73007c24dba4102e05f7e4c0904e5028c560cd
