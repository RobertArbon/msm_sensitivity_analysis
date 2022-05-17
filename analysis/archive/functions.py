import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
import seaborn as sns
# import plotly.express as px
from typing import List
import re
import patsy as pt
# import pymc3 as pm
from sklearn import preprocessing


cols = sns.color_palette('colorblind')

# PROTEIN_DIRS = ['1fme', '2f4k', '2jof', '2wav', 'cln025', 'gtt', 'prb', 'uvf', 'lambda', 'ntl9', 'nug2', 'a3d']
# PROTEIN_LABELS = ['BBA', 'Villin', 'Trp-cage', 'BBL', 'Chignolin', 'WW-domain', 'Protein-B', 'Homeodomain', '$\lambda$-repressor', 'NTL9', 'Protein-G', r'$\alpha$3D']
PROTEIN_DIRS = ['cln025',
 '2jof',
 '1fme',
 '2f4k',
 'gtt',
 'ntl9',
 '2wav',
 'prb',
 'uvf',
 'nug2',
 'a3d',
 'lambda']


PROTEIN_LABELS = ['Chignolin',
 'Trp-cage',
 'BBA',
 'Villin',
 'WW-domain',
 'NTL9',
 'BBL',
 'Protein-B',
 'Homeodomain',
 'Protein-G',
 '$\\alpha$3D',
 '$\\lambda$-repressor']


LETTERS = list('abcdefghijklmnopqrstuvwxyz')

FIG_DIR =Path(__file__).absolute().parents[1].joinpath('figures')
assert FIG_DIR.exists(), "path to figures directory doesn't exist"



def get_df(hdf_path: Path, results: str) -> pd.DataFrame:
    ts = pd.read_hdf(hdf_path, key=results)
    return ts


def get_hp_index(paths: List[Path], parent: int=0) -> List[int]:
    index = [re.findall('[0-9]+', x.parents[parent].stem)[0] for x in paths]
    index = [int(x) for x in index]
    return index


def get_results_df(results_paths, results: str):
    all_ts = []
    indices = get_hp_index(results_paths)
    for i, path in zip(indices, results_paths):
        ts = get_df(path, results)
        hp = pd.read_hdf(path, key='hp')
        hp['hp_index'] = i #path.parents[0].stem
        df = ts.join(hp).ffill()
        all_ts.append(df)
    all_ts = pd.concat(all_ts, axis=0)
    return all_ts


def timescale_gradient(ts_df: pd.DataFrame, x: str, log: bool = True, denom: str='one') -> pd.DataFrame:
    """
    Takes difference in timescales with respect to x.  Dataframe must be suitably subset before passing!
    """
    t = ts_df.loc[:, ['protein', x, 'value', 'hp_index', 'iteration']]
    
    dupes = t.duplicated(subset=['protein', x, 'hp_index', 'iteration'])
    
    assert not np.any(dupes.values), f'Duplicate values found: protein, {x}, index, & iteration, columns must be unique'
    
    if log:
        t['value'] = np.log(t['value'])
        
    t.sort_values(by=['protein', 'hp_index', 'iteration', x], inplace=True)
    t['delta_t'] = t.groupby(['protein', 'hp_index', 'iteration']).diff()['value']
    t['delta_x'] = t.groupby(['protein', 'hp_index', 'iteration']).diff()[x]
    
    if denom=='one':
        t['grad_t'] = t['delta_t']/1.0
    elif denom=='x':
        t['grad_t'] = t['delta_t']/t[x]
    elif denom=='delta_x':
        t['grad_t'] = t['delta_t']/t['delta_x']
        
    t.dropna(axis = 0, how = 'any', inplace = True)
    return t



def plot_timescales(ts, ax, label, errorbars=True):
#     lags = ts.lag.unique()
    nits = ts.num_its.unique()
    for its in nits:
        ax = plot_timescale(ts, its, ax, label, errorbars)
    return ax


def plot_timescale(ts, its, ax, label, errorbars=True):
    """
    ts = implied timescales df
    its = the implied timescale of interest
    """
    ix = ts.num_its == its
    x = ts.loc[ix, "lag"]
    y = ts.loc[ix, "median"]
    if errorbars:
        yerr = ts.loc[ix, ['del_lower', 'del_upper']].values.T
        ax.errorbar(x, y, yerr, lw=1, marker='o', elinewidth=2, color=cols[0])
    else:
        ax.plot(x, y, lw=1, label=label, color=cols[0])
    ax.set_yscale('log')
    return ax
    
    
def get_tau(tss, lags):
    """
    tss = implied timescales
    lags = associated markov lag time
    returned tau can't be the first lag
    """
    del_ts = tss[1:]/tss[:-1]
    min_grad = np.min(del_ts)
    tau = lags[1:][np.argmin(del_ts)]
    return tau, min_grad


def plot_tau(ts, its, ax, label):

    ix = ts.num_its == its
    lags = ts.loc[ix, "lag"].values
    tss = ts.loc[ix, "median"].values
    tau, grad = get_tau(tss, lags)
    ax.scatter(tau, grad, marker='o', s=50, alpha=1)
    return ax


def plot_taus(ts, ax, label):
    nits = ts.num_its.unique()
    for its in nits:
        ax = plot_tau(ts, its, ax, label)
    return ax



def gamma(alpha, beta):
    def g(x):
        return pm.Gamma(x, alpha=alpha, beta=beta)
    return g

def hcauchy(beta):
    def g(x):
        return pm.HalfCauchy(x, beta=beta)
    return g



# draws=1000, step=None, init='auto', n_init=200000, start=None, trace=None, chain_idx=0, chains=None, cores=None, tune=1000, progressbar=True, model=None, random_seed=None, discard_tuned_samples=True, compute_convergence_checks=True, callback=None, jitter_max_retries=10, *, return_inferencedata=None, idata_kwargs: Optional[dict] = None, mp_ctx=None, pickle_backend: str = 'pickle', **kwargs

def fit_gp(y, X, l_prior, eta_prior, sigma_prior, kernel_type='M52', bayes_kws=dict(draws=1000, tune=1000, chains=2, cores=1), prop_Xu=None):
    """
    function to return a pymc3 model
    y : dependent variable
    X : independent variables
    prop_Xu : number of inducing varibles to use. If None, use full marginal likelihood. If not none, use FTIC. 
    bayes_kw : kws for pm.sample
    X, y are dataframes. We'll use the column names. 
    """
    kernel_type = kernel_type.lower()
    with pm.Model() as model:
        # Covert arrays
        X_a = X.values
        y_a = y.values.flatten()
        X_cols = list(X.columns)
        
        # Globals
# #         prop_Xu = 0.1
#         l_prior = gamma(1, 0.05)
#         eta_prior = hcauchy(2)
#         sigma_prior = hcauchy(2)

        # Kernels
        # 3 way interaction
        eta = eta_prior('eta')
        cov = eta**2
        for i in range(X_a.shape[1]):
            var_lab = 'l_'+X_cols[i]
            if kernel_type=='rbf':
                cov = cov*pm.gp.cov.ExpQuad(X_a.shape[1], ls=l_prior(var_lab), active_dims=[i])
            if kernel_type=='exponential':
                cov = cov*pm.gp.cov.Exponential(X_a.shape[1], ls=l_prior(var_lab), active_dims=[i])
            if kernel_type=='m52':
                cov = cov*pm.gp.cov.Matern52(X_a.shape[1], ls=l_prior(var_lab), active_dims=[i])
            if kernel_type=='m32':
                cov = cov*pm.gp.cov.Matern32(X_a.shape[1], ls=l_prior(var_lab), active_dims=[i])

        # Covariance model
        cov_tot = cov 
        
        # Noise model
        sigma_n =sigma_prior('sigma_n')

        # Model
        if not (prop_Xu is None):
            # Inducing variables
            num_Xu = int(X_a.shape[0]*prop_Xu)
            Xu = pm.gp.util.kmeans_inducing_points(num_Xu, X_a)
            gp = pm.gp.MarginalSparse(cov_func=cov_tot, approx="FITC")
            y_ = gp.marginal_likelihood('y_', X=X_a, y=y_a, Xu=Xu, noise=sigma_n)
        else:
            gp = pm.gp.Marginal(cov_func=cov_tot)
            y_ = gp.marginal_likelihood('y_', X=X_a, y=y_a, noise=sigma_n)
            
        
        if not (bayes_kws is None):
            trace = pm.sample(**bayes_kws)
            result = trace
        else:
            mp = pm.find_MAP()
            result = mp
            
    return gp, result, model



def create_dmatrices(df, formula):
    y, X = pt.dmatrices(formula, data=df, return_type='dataframe')
    X = X.rename(columns=lambda x: re.sub('C','',x))
    return y, X


def scale_dmatrix(X, scaler=None):
    # scales matrices and returns scaler
    if scaler is None: 
        scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))
        scaler.fit(X.values)
    
    Xs = scaler.transform(X.values)
    Xs = pd.DataFrame(Xs, columns=[x+'_s' for x in X.columns], index=X.index)
    return Xs, scaler

def create_grid(search_space):
    # creates prediction grid from search space.
    Xnew = np.meshgrid(*search_space.values())
    Xnew = np.concatenate([x.reshape(-1, 1) for x in Xnew], axis=1)
    Xnew = pd.DataFrame(Xnew, columns=search_space.keys())
    for k, v in search_space.items():
        Xnew.loc[:, k] = Xnew.loc[:, k].astype(v.dtype)
    return Xnew

