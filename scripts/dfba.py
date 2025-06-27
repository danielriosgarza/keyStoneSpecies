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
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

MAX_SEC = 10
MONOD_CONSTANT = 1
TIME_SPAN = 100
N_POINTS = 100
INIT_BIOMASS = 0.01


def apply_env(model: cobra.Model, 
              env: dict[str, float], 
              upper_bound: float = MAX_SEC):
    for rid, conc in env.items():
        if model.reactions.has_id(rid):
            rxn = model.reactions.get_by_id(rid)
            rxn.lower_bound = -abs(conc)
            rxn.upper_bound = upper_bound
    return pfba(model, objective =model.objective)#.optimize()


def avoid_zero_met(states, changes):
    return np.where(states + changes < 0, -states, changes)


def get_growth(model, 
               env_names, 
               env_state, 
               sol0, 
               monod_c = MONOD_CONSTANT, 
               biomassName = 'biomass345'):
    uptake_bounds = {rid: env_state.get(rid, 0.0) / (env_state.get(rid, 0.0) + monod_c) for rid in env_names}
    sol = apply_env(model, uptake_bounds)
    mu = pfba.fluxes[biomassName]#sol.objective_value
    fluxes = {rid: (sol.fluxes.get(rid, 0.0) if model.reactions.has_id(rid) else 0.0) for rid in env_names}
    return mu / sol0, fluxes


def simulate(models_orig, env_names, env_initial,
             *, t_span=TIME_SPAN, n_points=N_POINTS, show_plots=False):
    n_models = len(models_orig)
    base_models = [m.copy() for m in models_orig]

    base_models_applied_env = [apply_env(m, env_initial) for m in base_models]
    sol0s = [max(m.objective_value, 1) for m in base_models_applied_env]
    
    def ode(t, y):
        biomasses = np.maximum(y[:n_models], 0.0)
        mets = np.maximum(y[n_models:], 0.0)
        env_state = dict(zip(env_names, mets))
        d_bio = np.zeros(n_models)
        d_mets = np.zeros_like(mets)
        for i in range(n_models):
            mu_norm, fluxes = get_growth(base_models[i], env_names, env_state, sol0s[i])
            d_bio[i] = mu_norm * biomasses[i]
            d_mets += biomasses[i] * np.array([fluxes[r] for r in env_names])
        d_mets = avoid_zero_met(mets, d_mets)
        return np.concatenate((d_bio, d_mets))

    y0 = np.concatenate(([INIT_BIOMASS] * n_models,
                         [env_initial.get(r, 0.0) for r in env_names]))
    
    sol = solve_ivp(ode, t_span, y0, t_eval=np.linspace(*t_span, n_points))

    if show_plots:
        # Figure 1: Biomass
        plt.figure(figsize=(6, 4))
        for i in range(n_models):
            plt.plot(sol.t, sol.y[i], label=f"Pop {i+1}")
        plt.xlabel("Time (h)"); plt.ylabel("Biomass")
        plt.title("Group Population Growth")
        plt.legend(); plt.tight_layout(); plt.show()

        # Figure 2: Metabolites
        plt.figure(figsize=(6, 4))
        for j, rid in enumerate(env_names):
            plt.plot(sol.t, sol.y[n_models + j], label=rid)
        plt.xlabel("Time (h)"); plt.ylabel("Concentration")
        plt.title("Environment Metabolite Dynamics")
        plt.legend(ncol=2, fontsize=7); plt.tight_layout(); plt.show()

    return sol
