from pathlib import Path
from typing import List, Optional

import pandas as pd
import numpy as np
import seaborn as sns
import pyemma as pm
import h5py

import matplotlib.pyplot as plt
import matplotlib as mpl
from pandas.api.types import CategoricalDtype


def feature_labeller(x):
    if x['feature__value'] == 'dihedrals':
        return 'dihed.'
    if x['feature__value'] == 'distances':
        if x['distances__transform'] == 'linear':
            return 'dist.'
        elif x['distances__transform'] == 'logistic':
            return 'logit(dist.)'


def vamps_and_features(summary: Path, hp_definitions: Path) -> pd.DataFrame:
    # Select vamp scores
    vamps = pd.DataFrame(pd.read_hdf(summary, key='vamps'))

    # Add labels for features
    hps = pd.read_hdf(hp_definitions, key='hyperparameters')
    features = hps.apply(feature_labeller, axis=1)
    print(hps.head())
    print(features.rename('feature').head())
    print(vamps.head())
    vamps = vamps.merge(right=features.rename('feature'), left_index=True, right_index=True, how='left')
    return vamps


def plot_vamps_ranked(summary: Path, hp_definitions: Path, output_directory: Path) -> None:
    with h5py.File(summary, 'r') as f:
        grp = f['summary']
        k = int(grp.attrs['chosen_k'])
        lag = int(grp.attrs['chosen_lag'])
        name = grp.attrs['protein_name']

    vamps = vamps_and_features(summary, hp_definitions)

    vamps.reset_index(level=['lag', 'process'], inplace=True)
    vamps = vamps.loc[(vamps.process == k) & (vamps.lag == lag), :]
    vamps.drop(labels=['lag', 'process'], inplace=True, axis=1)
    vamps.sort_index(inplace=True)

    # Cleaning for plotting
    cat_type = CategoricalDtype(categories=np.sort(vamps.feature.unique()), ordered=True)
    vamps['feature'] = vamps['feature'].astype(cat_type)

    # Create variables for plotting
    vamps['lb_err'] = vamps['median'] - vamps['lb']
    vamps['ub_err'] = -vamps['median'] + vamps['ub']
    vamps['rank'] = vamps['median'].rank(ascending=False)

    with sns.plotting_context('paper', font_scale=1.25):
        fig, ax = plt.subplots(1)
        ax.set_ylabel('VAMP-2')
        ax.set_xlabel('Rank')

        cols = sns.color_palette('colorblind')
        for code, cat in enumerate(vamps.feature.cat.categories):
            ix = vamps.feature == cat
            x = vamps.loc[ix, 'rank']
            y = vamps.loc[ix, 'median']
            yerr = vamps.loc[ix, ['lb_err', 'ub_err']].values.T
            ax.errorbar(x, y, yerr, lw=0, marker='o', elinewidth=1, c=cols[code], ms=3, label=cat)
            
        xlim = ax.get_xlim()
        ax.hlines(k, *xlim, color='k', label='Max. VAMP-2')
        ax.set_xlim(xlim)
        ylim = ax.get_ylim()
        ax.set_ylim(0.95, ylim[1])
        ax.legend(bbox_to_anchor=(1, 1))
        plt.tight_layout()

        out_fname = output_directory.joinpath(f'{name}_vamps_ranked.pdf')
        print(f"Saved image at {out_fname}")
        plt.savefig(out_fname, bbox_inches='tight')
        
        
def plot_vamps_vs_gap(summary: Path, hp_definitions: Path, output_directory: Path, n_ks: Optional[int] = None) -> None:
    with h5py.File(summary, 'r') as f:
        grp = f['summary']
        lag = int(grp.attrs['chosen_lag'])
        if n_ks is None:
            n_ks =int(grp.attrs['chosen_k']) + 1
        name = grp.attrs['protein_name']

    # Select vamp scores with gaps
    vamps = pd.read_hdf(summary, key='vamps')
    gaps = pd.read_hdf(summary, key='timescale_ratio')
    suffixes = ['_vamp', '_gap']
    df = vamps.merge(gaps, left_index=True, right_index=True, suffixes=suffixes)

    df.reset_index(level=['lag', 'process'], inplace=True)
    df = df.loc[(df.lag == lag), :]
    df.drop(labels=['lag'], inplace=True, axis=1)
    df.sort_index(inplace=True)

    # Add labels for features
    hps = pd.read_hdf(hp_definitions, key='hyperparameters')
    features = hps.apply(feature_labeller, axis=1)
    df = df.merge(right=features.rename('feature'), left_index=True, right_index=True, how='left')

    # Cleaning for plotting
    cat_type = CategoricalDtype(categories=np.sort(df.feature.unique()), ordered=True)
    df['feature'] = df['feature'].astype(cat_type)

    # Create variables for plotting
    for suffix in suffixes:
        df[f'lb_err{suffix}'] = df[f'median{suffix}'] - df[f'lb{suffix}']
        df[f'ub_err{suffix}'] = -df[f'median{suffix}'] + df[f'ub{suffix}']


    # Plots
    with sns.plotting_context('paper', font_scale=1.25):
        n_cols = 2
        n_rows = n_ks // n_cols + (n_ks % n_cols) - 1

        ks = np.arange(2, n_ks + 1)
        features = np.unique(df.feature.values)
        
        cols = sns.color_palette('colorblind', features.shape[0] + 4)

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6, 1+2*n_rows), sharey=False, sharex=True)
        if axes.ndim==1:
            axes = axes.reshape(1, -1)

        for i in range(ks.shape[0]):
            ax = axes.flatten()[i]

            for j, feat in enumerate(features):

                ix = (df.process == ks[i] ) & (df.feature == feat)
                x = df.loc[ix, 'median_gap'].values
                xerr = df.loc[ix, ['lb_err_gap', 'ub_err_gap']].values.T

                y = df.loc[ix, 'median_vamp'].values
                yerr = df.loc[ix, ['lb_err_vamp', 'ub_err_vamp']].values.T

                ax.errorbar(x, y, xerr=xerr, yerr=yerr, lw=0, elinewidth=0.5, marker='o', color=cols[j], alpha=0.5, label=feat, ms=5)
                ax.set_xscale('log')

                if i % n_cols == 0:
                    ax.set_ylabel(f'VAMP-2')

                if i // n_rows == 1:
                    ax.set_xlabel('Timescale ratio')
                    
                ax.annotate(text=f"$k={{{ks[i]}}}$", xy=(0.9, 0.95), ha='right', va='top', xycoords='axes fraction')

        axes[0, -1].legend(bbox_to_anchor=(1, 1), loc='upper left')
        plt.tight_layout()

        out_fname = output_directory.joinpath(f'{name}_vamps_vs_gap.pdf')
        print(f"Saved image at {out_fname}")
        plt.savefig(out_fname, bbox_inches='tight')