

import pyemma as pm
from typing import *
import numpy as np

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


def fit(ix, ftrajs_all, hp_dict, seed, lag, score_k):
    ftrajs = [ftrajs_all[i] for i in ix]

    ttrajs, tica_mod = tica(hp_dict, ftrajs)
    dtrajs, kmeans_mod = kmeans(hp_dict, ttrajs, seed)    

    msm = pm.msm.estimate_markov_model(dtrajs, lag=lag)

    vamp_eq = np.sum(msm.eigenvalues(score_k)**2) # if score_k = 3 then I want:  1 + lambda_2^2 + lambda_3^2
    evs = msm.eigenvalues(score_k+1) # if score_k = 3 then gap = lambda_3/lambda_4 (so need first 4 evs) 
    gap = evs[-2]/evs[-1]
    return vamp_eq, gap
