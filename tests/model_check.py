#!/usr/bin/env python
# coding: utf-8

# In[1]:

import click
import mdtraj as md
import pyemma as pm
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Optional, Union, Tuple
from pathlib import Path
import pickle

from msmtools.estimation import transition_matrix as _transition_matrix
from msmtools.analysis import timescales as _timescales


# In[2]:


def featurizer(hp_dict: Dict, traj_paths: List[str], top_path: str) -> List[np.ndarray]:
    if hp_dict['feature__value'] == 'dihedrals':
        assert hp_dict['dihedrals__which'] == 'all'
        def f(traj: md.Trajectory, **kwargs) -> np.ndarray:
            _, phi = md.compute_phi(traj)
            _, psi = md.compute_psi(traj)
            _, chi1 = md.compute_chi1(traj)
            _, chi2 = md.compute_chi2(traj)
            _, chi3 = md.compute_chi3(traj)
            _, chi4 = md.compute_chi4(traj)
            _, chi5 = md.compute_chi5(traj)
            ftraj = np.concatenate([phi, psi, chi1, chi2, chi3, chi4, chi5], axis=1)
            ftraj = np.concatenate([np.cos(ftraj), np.sin(ftraj)], axis=1)
            return ftraj

    elif hp_dict['feature__value'] == 'distances':
        def f(traj: md.Trajectory, **kwargs):
            scheme = kwargs['distances__scheme']
            transform = kwargs['distances__transform']
            centre = kwargs['distances__centre']
            steepness = kwargs['distances__steepness']
            ftraj, _ = md.compute_contacts(traj, scheme=scheme)
            if transform=='logistic':
                ftraj = 1.0/(1+np.exp(-steepness*(ftraj - centre)))
            return ftraj
    else:
        raise ValueError
    ftrajs = []
    for traj_path in traj_paths:
        traj = md.load(traj_path, top=top_path)
        ftrajs.append(f(traj, **hp_dict))
    return ftrajs


def tica(hp_dict: Dict[str, Union[float, int, str]], ftrajs: List[np.ndarray]) -> List[np.ndarray]:
    lag = hp_dict['tica__lag']
    stride = hp_dict['tica__stride']
    dim = hp_dict['tica__dim']
    tica = pm.coordinates.tica(ftrajs, lag=lag, dim=dim, kinetic_map=True)
    ttrajs = tica.get_output()
    return ttrajs, tica

def kmeans(hp_dict: Dict, ttrajs: List[np.ndarray], seed: int) -> List[np.ndarray]:
    k = hp_dict['cluster__k']
    max_iter = hp_dict['cluster__max_iter']
    stride = hp_dict['cluster__stride']
    kmeans = pm.coordinates.cluster_kmeans(ttrajs, k=k, max_iter=max_iter, stride=stride, fixed_seed=seed, n_jobs=1)
    dtrajs = kmeans.dtrajs
    return dtrajs, kmeans


def its(dtrajs: List[np.ndarray], lags: List[int], nits: int) -> np.ndarray:
    its_obj = pm.msm.timescales_msm(dtrajs=dtrajs, lags=lags, nits=nits)
    return its_obj.timescales


def score(dtrajs: List[np.ndarray], lags: List[int], nits: int) -> np.ndarray:
    all_vs = []
    for lag in lags: 
        m = pm.msm.estimate_markov_model(dtrajs, lag=lag)
        vs = np.array([m.score(dtrajs, score_k=k) for k in range(2, nits+2)])
        vs = vs.reshape(1, -1)
        all_vs.append(vs)
    all_vs = np.concatenate(all_vs, axis=0)
    return all_vs
        


def bootstrap(ftrajs: List[np.ndarray], rng: np.random.Generator) -> List[np.ndarray]:
    probs = np.array([x.shape[0] for x in ftrajs])
    probs = probs/np.sum(probs)
    ix = np.arange(len(ftrajs))
    new_ix = rng.choice(ix,size=len(ftrajs), p=probs, replace=True)
    return [ftrajs[i] for i in new_ix]



def summarise(df):
    df_summary = df.groupby(['hp_ix', 'lag', 'process']).agg(median=(0, lambda x: np.quantile(x, 0.5)),
                                                                   lb=(0, lambda x: np.quantile(x, 0.025)),
                                                                   ub=(0, lambda x: np.quantile(x, 0.975)), 
                                                                   count =(0, lambda x: x.shape[0]-x.isna().sum()))
    return df_summary


def samples_to_summary(samples: np.ndarray, lags: List[int],  hp_ix: int)-> pd.DataFrame: 
    """
    samples=np.ndarray[lagtime, process, bs_sample]
    """
    df = pd.concat({(hp_ix, lags[i], j+2): pd.DataFrame(samples[i, j, :]) for i in range(samples.shape[0]) for j in range(samples.shape[1])})
    df.index.rename(('hp_ix', 'lag', 'process', 'bs_ix'), inplace=True)
    df_summary = summarise(df)
    return df_summary
                    

@click.command()
@click.argument("protein", type=str)
@click.argument("hp_ix", type=int)
def run(protein, hp_ix):
    seed = 49587
    n_bootstraps = 100
    nits=25
    
    rng = np.random.default_rng(seed)
    lags = list(range(1, 102, 10))


    hps = pd.read_hdf('../data/msms/hpsample.h5')
    top_path = f'/home/rob/Data/DESRES/DESRES-Trajectory_{protein.upper()}-0-protein/{protein.upper()}-0-protein/protein.pdb'
    traj_paths = list(Path('/home/rob/Data/DESRES/').rglob(f'*{protein.upper()}*/**/*.xtc'))
    traj_paths = [str(x) for x in traj_paths]
    traj_paths.sort()
    assert traj_paths

    source_ts = pd.DataFrame(pd.read_hdf(f'../data/msms/{protein}/summary.h5', key='timescales'))
    source_vs = pd.DataFrame(pd.read_hdf(f'../data/msms/{protein}/summary.h5', key='vamps'))


    ftrajs_all = featurizer(hps.iloc[hp_ix, :].to_dict(), traj_paths, top_path)
    # Bootstrap results
    ts_samples = []
    vs_samples = []
    for i in range(n_bootstraps):
        print(i, end=', ')
        ftrajs = bootstrap(ftrajs_all, rng)

        assert len(ftrajs) == len(ftrajs_all)
        ttrajs, tica_mod = tica(hps.iloc[hp_ix, :].to_dict(), ftrajs)
        dtrajs, kmeans_mod = kmeans(hps.iloc[hp_ix, :].to_dict(), ttrajs, seed)
        ts = its(dtrajs, lags, nits=nits)
        vs = score(dtrajs, lags, nits=nits)

        ts_samples.append(ts[..., np.newaxis])
        vs_samples.append(vs[..., np.newaxis])

    # Summarise values
    ts_samples = np.concatenate(ts_samples, axis=-1)
    vs_samples = np.concatenate(vs_samples, axis=-1)

    target_ts = samples_to_summary(ts_samples, lags, hp_ix)
    target_vs = samples_to_summary(vs_samples, lags,  hp_ix)

    # Compare to msmsense values
    comp_ts = pd.merge(target_ts, source_ts, left_index=True, right_index=True, how='left')
    comp_vs = pd.merge(target_vs, source_vs, left_index=True, right_index=True, how='left')

    # Plot vamps comaprison
    fig, axes = plt.subplots(10, 2, figsize=(8, 20), sharex=True, sharey=True)
    num_plots = axes.flatten().shape[0]
    for k, v in comp_vs.groupby(['hp_ix', 'process']):
        if k[1]-2 < num_plots:
            lags = v.index.get_level_values(level=1)
            ax = axes.flatten()[k[1]-2]
            ax.plot(lags, v['median_y']/v['median_x'], label=f'median {k[1]}', marker='o')
            ax.plot(lags, v['lb_y']/v['lb_x'], label=f'lower {k[1]}', marker='o')
            ax.plot(lags, v['ub_y']/v['ub_x'], label=f'upper {k[1]}', marker='o')
            ax.grid()
            ax.legend()

    plt.tight_layout()
    plt.savefig(f"hp_{hp_ix}_vamps_comparison.pdf", bbox_inches='tight')

    # Plot timescales comparison
    fig, axes = plt.subplots(10, 2, figsize=(8, 20), sharex=True, sharey=True)
    num_plots = axes.flatten().shape[0]
    for k, v in comp_ts.groupby(['hp_ix', 'process']):
        if k[1]-2 < num_plots:
            lags = v.index.get_level_values(level=1)
            ax = axes.flatten()[k[1]-2]
            ax.plot(lags, v['median_y']/v['median_x'], label=f'median {k[1]}', marker='o')
            ax.plot(lags, v['lb_y']/v['lb_x'], label=f'lower {k[1]}', marker='o')
            ax.plot(lags, v['ub_y']/v['ub_x'], label=f'upper {k[1]}', marker='o')
            ax.grid()
            ax.legend()

    plt.tight_layout()
    plt.savefig(f"hp_{hp_ix}_timescales_comparison.pdf", bbox_inches='tight')


# In[3]:

if __name__ == '__main__':
    run()



# In[7]:


# df = target_ts.copy(deep=True)

# df = df.droplevel(level=0)

# for i in range(100):
#     df2 = source_ts.copy(deep=True)
#     df2 = df2.loc[(i, slice(None), slice(None)), :]
#     df2 = df2.droplevel(level=0)
#     tmp = df.merge(df2, left_index=True, right_index=True, how='left')
#     # print(tmp.head())
#     diff = np.abs(tmp['median_x']/tmp['median_y']).mean()
#     print(i, diff)
    


# In[ ]:




