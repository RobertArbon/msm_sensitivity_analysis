# MSM parameters

import numpy as np
from pathlib import Path
import matplotlib as mpl
import seaborn as sns

protein = 'BBA'
seed = 49587
rng = np.random.default_rng(seed)
lags = list(range(1, 102, 10))
lags_test = [1, 5, 10, 20, 50, 100, 200, 500, 1000]
n_bootstraps = 100
nits= 25
lag = 41
hp_ix = 227
    # BBA 0-135 random, 136-235 optimised 
    # BBA before 1st 24 (24), 2nd 124 (262), after 227
    # BBA dihedral before 6 (6), after 141 
    # BBA distance before 58 (85), after 58
n_ts = 10
test_labels = '1FME'

# PCCA+ parameters
n_sets = 2
core_membership_cutoff = 0.9

# Input path 
summary_path = r'../BBA/t2.h5'
top_path = '../../../1fme/protein.pdb'
ref_path = [r'../../../1fme/1fme.pdb' ]

traj_paths = list(Path(r'../../../1fme').rglob(f'1FME-*.xtc'))
traj_paths = [str(x) for x in traj_paths]
traj_paths.sort()

# Drawing parameters
mpl.rcParams['savefig.bbox'] = 'tight'
sns.set_style("white")
sns.set_style({'font.family':'sans-serif', 'font.serif':'Arial'})