# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 15:48:49 2025

@author: drgarza
"""

import os
from pathlib import Path
os.environ['GRB_LICENSE_FILE'] = os.path.join(Path(os.getcwd()), 'licenses', 'gurobi.lic')

import cobra
from cobra.flux_analysis import pfba
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Constants
MAX_SEC = 10
MONOD_CONSTANT = 10
YIELD_FACTOR = 0.1
TIME_SPAN = (0, 120)
N_POINTS = 100
INIT_BIOMASS = 0.01
q_in = 2.08
VOLUME = 100


def apply_env(model: cobra.Model,
              env: dict[str, float],
              upper_bound: float = MAX_SEC):
    for rid, conc in env.items():
        if model.reactions.has_id(rid):
            rxn = model.reactions.get_by_id(rid)
            rxn.lower_bound = -abs(conc)
            rxn.upper_bound = upper_bound
    model.reactions.get_by_id('EX_h2o(e)').lower_bound=-10
    return pfba(model, objective=model.objective)


def avoid_zero_met(states, changes):
    return np.where(states + changes < 0, -states, changes)


def get_growth(model,
               env_ids,
               env_state,
               monod_c=MONOD_CONSTANT,
               biomassName='biomass345'):
    uptake_bounds = {
        rid: env_state.get(rid, 0.0) / (env_state.get(rid, 0.0) + monod_c[idx])
        for idx,rid in enumerate(env_ids)
        }
    #print(uptake_bounds)
    sol = apply_env(model, uptake_bounds)
    mu = sol.fluxes[biomassName]
    fluxes = {}
    for rid in env_ids:
        if model.reactions.has_id(rid):
            v = sol.fluxes.get(rid, 0.0)
            fluxes[rid] = v*YIELD_FACTOR if v > 0 else v
        else:
            fluxes[rid] = 0.0
    
    return mu, fluxes


def simulate(models_orig,
             env_ids,
             env_names,
             env_initial,
             monod_c,
             *,
             t_span=TIME_SPAN,
             n_points=N_POINTS,
             show_plots=False,
             biomassNames=None,
             feed_medium=None):

    n_models = len(models_orig)
    

    # Default feed is the initial environment
    feed_medium = feed_medium or env_initial

    
    V = VOLUME

    def ode(t, y):
        biomasses = np.maximum(y[:n_models], 0.0)
        mets = np.maximum(y[n_models:], 0.0)
        env_state = dict(zip(env_ids, mets))

        d_bio = np.zeros(n_models)
        d_mets = np.zeros_like(mets)

        for i in range(n_models):
            mu, fluxes = get_growth(models_orig[i], 
                                    env_ids, 
                                    env_state, 
                                    monod_c = monod_c,
                                    biomassName=biomassNames[i])
            d_bio[i] = mu * biomasses[i]
            d_mets += biomasses[i] * np.array([fluxes[r] for r in env_ids])

        # Add dilution terms
        d_bio -= (q_in / V) * biomasses
        d_mets += (q_in / V) * (np.array([feed_medium[r] for r in env_ids]) - mets)

        d_mets = avoid_zero_met(mets, d_mets)

        return np.concatenate((d_bio, d_mets))

    # Initial condition
    y0 = np.concatenate((
        [INIT_BIOMASS] * n_models,
        [env_initial.get(r, 0.0) for r in env_ids]
    ))

    sol = solve_ivp(ode, t_span, y0, t_eval=np.linspace(*t_span, n_points))

    if show_plots:
        # Plot biomass
        plt.figure(figsize=(6, 4))
        for i in range(n_models):
            plt.plot(sol.t, sol.y[i], label=f"Pop {i+1}")
        plt.xlabel("Time (h)")
        plt.ylabel("Biomass")
        plt.title("Group Population Growth")
        plt.legend()
        plt.tight_layout()
        plt.show()

        # Plot metabolites (only those that change)
        plt.figure(figsize=(6, 4))
        threshold = 3.0  # minimum change to be considered "meaningful"
        plotted_any = False
        for j, (rid, label) in enumerate(zip(env_ids, env_names)):
            trace = sol.y[n_models + j]
            print(f"{j:3d} {rid:20s} {label:20s} Δ={trace.max() - trace.min():.3f}")
            if (trace.max() - trace.min()) > threshold:
                plt.plot(sol.t, trace, label=label)
                plotted_any = True

        plt.xlabel("Time (h)")
        plt.ylabel("Concentration")
        plt.title("Environment Metabolite Dynamics")
        if plotted_any:
            plt.legend(ncol=2, fontsize=7)
        else:
            plt.text(0.5, 0.5, "No changing metabolites",
                     ha='center', va='center', transform=plt.gca().transAxes)
        plt.tight_layout()
        plt.show()

    return sol
