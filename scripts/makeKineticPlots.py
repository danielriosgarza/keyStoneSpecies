# -*- coding: utf-8 -*-
"""
Created on Sun Jun 29 13:00:35 2025

@author: drgarza
"""

import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import math

from aquarel import load_theme


theme = load_theme("boxy_light")
theme.apply()

def loadBacterialKinetics(filePath):
    TIME_POINTS = 100
    BACTERIA_N = 11
    time = np.zeros(TIME_POINTS)
    
    with open(filePath) as f:
        bacteria = f.readline().rstrip('\n').split('\t')[1:BACTERIA_N + 1]
        data = np.zeros((TIME_POINTS, BACTERIA_N))
        counter = 0
        for line in f:
            a = line.rstrip('\n').split('\t')
            time[counter] = float(a[0])
            data[counter] = np.array(a[1:BACTERIA_N + 1]).astype(np.float32)
            counter+=1
    return time, data.T, bacteria
        




def plot_kinetics_grid(kinetics, filePath=None, n_cols=4, subplot_size=3.5, title_prefix='Sample'):
    """
    Plot bacterial kinetics clearly, from top to bottom and left to right.

    Parameters:
    - kinetics: dict {sample_id: {'time': ..., 'data': ..., 'bac': ...}}
    - filePath: path to save figure (optional)
    - n_cols: number of columns
    - subplot_size: size of each subplot
    - title_prefix: optional title prefix
    """
    n_samples = len(kinetics)
    n_rows = math.ceil(n_samples / n_cols)

    figsize = (n_cols * subplot_size, n_rows * subplot_size + 2)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, sharex=False, sharey=False)

    # Define 11 clearly distinguishable colors
    distinguishable_colors = [
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
        '#FF1493',  # Deep Pink
    ]

    sorted_samples = sorted(kinetics.keys())
    axes = axes.reshape((n_rows, n_cols))

    idx = 0
    for col in range(n_cols):
        for row in range(n_rows):
            if idx >= n_samples:
                axes[row, col].axis('off')
                continue

            ax = axes[row, col]
            sample_id = sorted_samples[idx]
            time = kinetics[sample_id]['time']
            data = kinetics[sample_id]['data']
            labels = kinetics[sample_id]['bac']

            for i, label in enumerate(labels):
                color = distinguishable_colors[i % len(distinguishable_colors)]
                ax.plot(time, data[i], label=label, lw=4, color=color)

            ax.tick_params(axis='both', labelsize=12)
            idx += 1

    # Global X and Y labels
    fig.text(0.5, 0.05, 'Time (h)', ha='center', fontsize=25)
    fig.text(0.06, 0.5, 'Biomass', va='center', rotation='vertical', fontsize=25)

    # Global legend below subplots
    #handles, labels = axes[0,0].get_legend_handles_labels()
    #fig.legend(handles, labels,
    #           loc='lower center',
    #           fontsize=25,
    #           ncol=min(len(labels), 6),
    #           bbox_to_anchor=(0.53, 0.03))

    # Explicit subplot spacing
    plt.subplots_adjust(left=0.12, right=0.98, top=0.97, bottom=0.1, hspace=0.2, wspace=0.2)

    if filePath:
        plt.savefig(filePath, dpi=450, bbox_inches='tight')

    plt.show()







rootPath = os.path.join(Path(os.getcwd()).parents[0], 'files', '100_env_all_models')

kinetics = {}
counter = 0

for i in range(12):
    

    for i in range(1,16):
        kinetics[counter] = {}
        kinetics[counter]['time'],kinetics[counter]['data'], kinetics[counter]['bac'] = loadBacterialKinetics(os.path.join(rootPath, str(i) + '.tsv'))
        counter+=1


# for i in range(1,25):
#     kinetics[counter] = {}
#     kinetics[counter]['time'],kinetics[counter]['data'], kinetics[counter]['bac'] = loadBacterialKinetics(os.path.join(rootPath, str(i) + '.tsv'))
#     counter+=1

# plot_kinetics_grid(kinetics, n_cols=5, filePath=os.path.join(Path(os.getcwd()).parents[0], 'files', 'Figures', '100_envs.png'))

plot_kinetics_grid(kinetics, n_cols=12, filePath=os.path.join(Path(os.getcwd()).parents[0], 'files', 'Figures', 'mock.png'))