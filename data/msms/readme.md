This directory contains run scripts for both scoring and validating models with `msmsense` and analysing the data using `sensetools` (i.e., this repo)

Create a conda env: 

```
conda create -n msmsense python==3.9 pyemma pandas numpy click mdtraj -y
conda activate msmsense
```

install this repo and [msmsense](https://github.com/RobertArbon/msm_sensitivity) : 

```
cd ../../
pip install -e .

cd /my/repo/directory/

git clone https://github.com/RobertArbon/msm_sensitivity
cd msm_sensitivity
pip install -e . 
```

`0_convert_trajs.sh` - uses `mdconvert` to export dcds to xtcs with a stride of 5. 

`1_sample_hyperparameters.sh` - samples hyperparameters.  This is done naively and will result in an unbalanced set of hyperaprameters.  Use `truncate_hp_sample.ipynb` in the `analysis` directory in order to truncate the sample to a more balanced set. 

`2_score_hyperparameters.sh` - the main function.  Scores hyperparameters at a set of lags using vamp2 score.  Will score all resolvable processes up to 20. 

`3_summarise_hyperparameter.sh` - summarises the scores - the vamp scores and timescales. Also contains model selection tools but these are a bit useless for the paper now. 

`4_sample_states.sh` - coarse grain msm and sample metastable states.  Note - the ordering of the states is arbitrary. 

`5_sample_evs.sh` - samples eigenvectors. Finds frames that lie exclusively (as far as possible) along a given eigenvector. 

`6_dump_dtrajs.sh` - dumps out the discretized trajectories for use in other analysis.  

`6_project_evs.sh` - projects one set of eigenvectors onto a different model.  Use after `5_sample_evs.sh`. 

`7_cktests.sh` - runs a CK test - not tested yet. 




