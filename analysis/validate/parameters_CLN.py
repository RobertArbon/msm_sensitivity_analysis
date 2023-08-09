# MSM parameters

import numpy as np
from pathlib import Path
import matplotlib as mpl
import seaborn as sns

protein = 'CLN'
seed = 49587
rng = np.random.default_rng(seed)
lags = list(range(1, 102, 10))
lags_test = [1, 5, 10, 20, 50, 100, 200, 500, 1000]
n_bootstraps = 100
nits= 25
lag = 31
hp_ix = 52
    # CLN 0-130 random, 131-231 optimised
    # CLN before 52, after 218
n_ts = 10
test_labels = '5AWL'
    
# PCCA+ parameters
n_sets = 2
core_membership_cutoff = 0.9

# Input path 
summary_path = r'../CLN/maximize_t2.h5'
top_path = '../../../chignolin/protein.pdb'
ref_path = [r'../../../chignolin/5AWL_H.pdb']

traj_paths = list(Path(r'../../../chignolin').rglob(f'CLN025-0-protein-*.xtc'))
traj_paths = [str(x) for x in traj_paths]
traj_paths.sort()
    
# Drawing parameters
mpl.rcParams['savefig.bbox'] = 'tight'
sns.set_style("white")
sns.set_style({'font.family':'sans-serif', 'font.serif':'Arial'})