import numpy as np
from typing import Dict, List, Optional, Union, Tuple
from pathlib import Path

protein = 'CLN'
db_name = protein
storage_name = "sqlite:///../{}/{}.db".format(protein, protein)
study_name = 'vampeq2+vamp23gap'

n_hours = 24
timeout = n_hours*60*60
n_trials = 100
n_boot = 20
seed = 49587
new_gamma = False

DATA_PATH = '/home/rzhu/Desktop/projects/msm_opt/chignolin'
lag = 31
top_path = f'protein.pdb'
traj_path_glob = 'CLN025-*.xtc'
score_k = 2
process = 2

summary_path = f'../CLN/summary.h5'
hp_path = f'../CLN/hpsample_stride1.h5'

fig_path = Path(f"../{protein}/{study_name}")
fig_path.mkdir(parents=True, exist_ok=True)

def hyperopt_gamma(x: int) -> int:
    return min(int(np.ceil(0.25 * np.sqrt(x))), 25)