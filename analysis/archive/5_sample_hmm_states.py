import pyemma as pm
import mdtraj as md
import numpy as np
from pathlib import Path
import pandas as pd
from typing import Dict, Union, List
from msmsense.featurizers import dihedrals, distances
from msmsense.bootstrap_cmatrices import get_sub_dict, get_trajs
from functools import partial
import functions as funcs
import sys
from pyemma import config as pmconfig
pmconfig.show_progress_bars = False 

import pickle
import time
import seaborn as sns



def get_feature_dict(df, row_num):
    row_dict = df.filter(regex='__', axis=1).iloc[row_num, :].to_dict()
    feature_dict = get_sub_dict(row_dict, 'feature')
    if feature_dict['value'] == 'distances':
        feature_dict.update(get_sub_dict(row_dict, 'distances'))
    if feature_dict['value'] == 'dihedrals':
        feature_dict.update(get_sub_dict(row_dict, 'dihedrals'))
    return feature_dict

def get_kws_dict(df, row_num, kws):
    row_dict = df.filter(regex='__', axis=1).iloc[row_num, :].to_dict()   
    kws_dict = get_sub_dict(row_dict, kws)
    return kws_dict

def set_proper_dtypes(df):
    """
    forgot to save integers as integers. Only the distances feature columns have true floats. 
    """
    float_cols = list(df.filter(regex='distances.*', axis=1)) + list(df.filter(regex='.*vamp.*', axis=1)) 
    float_cols = float_cols + list(df.filter(regex='.*gap.*', axis=1)) 
    potential_integer_cols = df.columns.difference(float_cols)
    for col in potential_integer_cols:
        if str(df[col].dtype) != 'object':
            df[col] = df[col].astype(int)
    return df

def get_trajs_top(traj_dir: Path, protein_dir: str, rng: Union[np.random.Generator, None]=None):
    trajs = list(traj_dir.rglob(f"*{protein_dir.upper()}*/*.xtc"))
    trajs.sort()
    if rng is not None:
        ix = rng.choice(np.arange(len(trajs)), size=len(trajs), replace=True)
        trajs = [trajs[i] for i in ix]
    
    top = list(traj_dir.rglob(f"*{protein_dir.upper()}*/*.pdb"))[0]
    
    return {'trajs': trajs, 'top': top}
    
    
def get_random_traj(trajs: List[md.Trajectory], num_frames: int, rng: np.random.Generator)-> md.Trajectory: 
    traj_ix = np.arange(len(trajs))
    frame_ix = [np.arange(traj.n_frames) for traj in trajs]
    
    rand_ix = [(ix, rng.choice(frame_ix[ix])) for ix in rng.choice(traj_ix, size=num_frames)]
    rand_traj = md.join([trajs[x[0]][x[1]] for x in rand_ix])
    return rand_traj
    


# In[3]:


class MSM(object):
    
    def __init__(self, lag: int, num_evs: int, top: str, feature_kws: Dict[str, Union[str, int, float]], 
                 tica_kws: Dict[str, Union[str, int, float]], cluster_kws: Dict[str, Union[str, int, float]], seed: int):
        """
        Defines the whole MSM pipeline.
        lag: markov lag time 
        num_evs: number of eigenvectors in VAMP score. This includes stationary distribution. note: all projections are done onto the processes 1 - num_evs, i.e., exclude the stationary distribution (process 0)
        
        """
        self.lag = lag
        self.num_evs = num_evs
        self.top = top
        self.feature_kws = feature_kws
        self.tica_kws = tica_kws
        self.cluster_kws = cluster_kws
        self.featurizer = None
        self._set_featurizer()
        self.seed = seed

        self.msm = None
        
    def _set_featurizer(self):
        feature_kws = self.feature_kws.copy()
        feature = feature_kws.pop('value')
        
        if feature == 'distances':
            def f(traj_paths, top):
                dists = []
                for path in traj_paths:
                    traj = md.load(path, top=self.top)
                    dists.append(distances([traj], **feature_kws)[0])
                return dists
            
            self.featurizer = f
            
        elif feature == 'dihedrals':
            def f(traj_paths, top):
                diheds = []
                for path in traj_paths:
                    traj = md.load(path, top=self.top)
                    diheds.append(dihedrals([traj], **feature_kws)[0])
                return diheds
            
            self.featurizer = f
        else:
            raise NotImplementedError('Unrecognized feature')
        

    def fit(self, trajs):
        ftrajs = self.featurizer(trajs, self.top)
        tica = pm.coordinates.tica(data=ftrajs, **self.tica_kws)
        ttrajs = tica.get_output()
        cluster = pm.coordinates.cluster_kmeans(data=ttrajs, **self.cluster_kws, fixed_seed=self.seed)
        dtrajs = cluster.dtrajs
        self.msm = pm.msm.estimate_markov_model(dtrajs=dtrajs, lag=self.lag)

        

def model_kwargs(mod_defs, top, row_num, seed):

    lag = int(mod_defs.chosen_lag.values[row_num])
    num_evs = int(mod_defs.new_num_its.values[row_num])
    feat_kws = get_feature_dict(mod_defs, row_num)
    tica_kws = get_kws_dict(mod_defs, row_num, 'tica')
    cluster_kws = get_kws_dict(mod_defs, row_num, 'cluster')

    kwargs = dict(lag = lag, num_evs=num_evs, top=top, feature_kws=feat_kws, tica_kws=tica_kws, cluster_kws=cluster_kws, seed=seed)
    return kwargs


def fit_model(kwargs, trajs):
    model = MSM(**kwargs)
    model.fit(trajs)
    return model

def cktest(model, mlags=10):
    pass


def get_trajectories(traj_dir, protein_dir, rng = None):
    traj_paths = get_trajs_top(traj_dir, protein_dir, rng)
    traj_paths_str = dict(top=str(traj_paths['top']), trajs=[str(x) for x in traj_paths['trajs']])
    top = md.load(str(traj_paths['top']))
#     trajs = [md.load(str(x), top=top) for x in traj_paths['trajs']]
    return traj_paths_str['trajs'], traj_paths_str['top']
             


def get_model_defs(all_models, protein, feature):
    mod_defs = all_models.loc[all_models.protein==protein, :].copy()
    if feature in mod_defs['feature'].values:
        row_num = np.where(mod_defs['feature'].values==feature)[0][0]
    else: 
        row_num = None
    return mod_defs, row_num


def sample_states(traj_dir, protein_dir, rng, mod_defs, row_num, num_structures=100):
    
    

    trajs, top = get_trajectories(traj_dir, protein_dir)
    kwargs = model_kwargs(mod_defs, top, row_num, seed)
    model = fit_model(kwargs, trajs)
    
    num_procs = model.num_evs
    hmm = model.msm.coarse_grain(num_procs)
            
    sd = hmm.stationary_distribution

    # Sort states and get distributions
    state_ix = np.argsort(sd)[::-1]
    sd = sd[state_ix]
    ms_dist = hmm.metastable_distributions
    ms_dist = ms_dist[state_ix, :]

    # Get samples
    samples = [pm.coordinates.save_traj(trajs, idist, outfile=None, top=top)
              for idist in model.msm.sample_by_distributions(ms_dist, num_structures)] 
    
    
    return samples

    


if __name__ == '__main__': 

    traj_dir = Path('/Volumes/REA/MD/12FF/strided/')
    m1_sel = set_proper_dtypes(pd.read_hdf('./summaries/m1_model_selection.h5'))
    m2_sel = set_proper_dtypes(pd.read_hdf('./summaries/m2_model_selection.h5'))
    m3_sel = set_proper_dtypes(pd.read_hdf('./summaries/m3_model_selection.h5'))
    model_selections={'m1': m1_sel, 'm2': m2_sel, 'm3': m3_sel}
    
    prot_dict = dict(zip(funcs.PROTEIN_LABELS, funcs.PROTEIN_DIRS))
    prot_pdbs = {'Chignolin': '5awl.pdb', 'BBA': '1fme.pdb', 'Trp-cage': '2jof.pdb',  
                 'Villin': '2f4k.pdb', 'WW-domain': '2f21.pdb' , 'BBL': '2wxc.pdb' , 
                 'Homeodomain': '2p6j.pdb', 'Protein-B': '1prb.pdb'}
    
    num_structures = 100 
    seed = 12098345
    
    features = m1_sel.feature.unique()
    
    rng = np.random.default_rng(seed)
    
    protein = sys.argv[1]
    
    # Setup output directory
    protein_dir = prot_dict[protein]
    root_dir = Path(f"compare_structures/{protein}")
    root_dir.mkdir(exist_ok=True, parents=True)
    
    # loop over feature
    for feature in features: 
        print(feature)
        # loop over selection methods
        for method, selection in model_selections.items(): 
            print(method)
            # Get model definition
            mod_defs, row_num = get_model_defs(selection, protein, feature)
            if not row_num is None:
                # bootstrap cktest
                samples = sample_states(traj_dir, protein_dir, rng, mod_defs, row_num, num_structures=num_structures)
                for i, traj in enumerate(samples):
                    traj.save_xtc(str(root_dir.joinpath(f"{method}_model_{feature}_state_{i}.xtc")))

    
