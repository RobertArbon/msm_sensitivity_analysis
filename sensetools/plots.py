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

LETTERS = list('abcdefghiklmnopqrstuvwxy')


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

    selections = pd.read_hdf(summary, key='model_selection')

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
        markers = {'worst': 'P', 'fixed_k': 'o', 'timescale_gap': 'X'}

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
                if i == 0:
                    ax.errorbar(x, y, xerr=xerr, yerr=yerr, lw=0, elinewidth=0.5, marker='o', color=cols[j], alpha=0.5, label=feat, ms=5)
                else:
                    ax.errorbar(x, y, xerr=xerr, yerr=yerr, lw=0, elinewidth=0.5, marker='o', color=cols[j], alpha=0.5,
                                ms=5)
                ax.set_xscale('log')

                if i % n_cols == 0:
                    ax.set_ylabel(f'VAMP-2')

                if i // n_rows == 1:
                    ax.set_xlabel('Timescale ratio')
                    
                ax.annotate(text=f"$k={{{ks[i]}}}$", xy=(0.9, 0.95), ha='right', va='top', xycoords='axes fraction')

        all_methods = selections['method'].unique()
        added_methods = []
        for i, row in selections.iterrows():
            k = row['process']
            x = row['median_gap']
            y = row['median_vamp']
            label = row['method']
            hp_ix = row['hp_ix']
            feature = row['feature']
            feat_ix = np.where(features == feature)[0][0]
            ax_ix = np.where(ks == k)[0][0]
            ax = axes.flatten()[ax_ix]
            if (label not in added_methods) and (feature == features[0]):
                ax.scatter(x, y, marker=markers[label], color=cols[feat_ix],
                           label=label.replace('_', ' ').capitalize(),
                           edgecolor='k', zorder=10, alpha=1, s=50)
                added_methods.append(label)
            else:
                ax.scatter(x, y, marker=markers[label], color=cols[feat_ix],
                           edgecolor='k', zorder=10, alpha=1, s=50)
        handles = []
        labels = []
        for ax in axes.flatten():
            h, l = ax.get_legend_handles_labels()
            handles.extend(h)
            labels.extend(l)

        axes[0, -1].legend(handles=handles, labels=labels, bbox_to_anchor=(1, 1), loc='upper left')
        plt.tight_layout()

        out_fname = output_directory.joinpath(f'{name}_vamps_vs_gap.pdf')
        print(f"Saved image at {out_fname}")
        plt.savefig(out_fname, bbox_inches='tight')


def implied_timescales(summary: Path, output_directory: Path) -> None:

    with h5py.File(summary, 'r') as f:
        grp = f['summary']
        lag = int(grp.attrs['chosen_lag'])
        name = grp.attrs['protein_name']

    selections = pd.read_hdf(summary, key='model_selection')
    timescales = pd.read_hdf(summary, key='timescales')
    gaps = pd.read_hdf(summary, key='timescale_ratio')

    hp_ixs = selections['hp_ix'].unique()

    def plot(hp_ix, timescales, gaps, use_lag, out_dir, protein_name):
        ts = timescales.reset_index()
        ts = ts.loc[ts.hp_ix == hp_ix, :]

        ts_ratio = gaps.reset_index()
        ts_ratio = ts_ratio.loc[(ts_ratio.hp_ix == hp_ix) & (ts_ratio.lag == use_lag), :]

        ts_kwargs = dict(lw=1, elinewidth=2, capsize=3, marker='o', markersize=4)
        ratio_kwargs = dict(lw=0, elinewidth=2, capsize=3, marker='o', markersize=4)
        jitter_sd = 0.5

        with sns.plotting_context('paper', font_scale=1.25):
            sns.set_style('whitegrid')
            # Setup axes
            n_rows, n_cols = 2, 1
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(6, 8), sharey='row', sharex='row')
            ts_axes = axes[0]
            ratio_axes = axes[1]

            # Implied timescales
            for proc, df in ts.groupby('process'):
                x = df['lag'] + np.random.normal(loc=0, scale=jitter_sd, size=df.shape[0])
                y = df['median']
                yerr = (df['median']-df['lb'], df['ub'] - df['median'])
                ts_axes.errorbar(x, y, yerr, **ts_kwargs, label=f"{int(proc)}")
            ts_axes.set_yscale('log')
            ts_axes.set_ylabel('Timescale (ns)')
            ts_axes.set_xlabel('Lag time (ns)')

            # Chosen lag
            ylim = ts_axes.get_ylim()
            ts_axes.vlines(use_lag, *ylim, color='k', label=r'$\tau_{\mathrm{M}}$')
            ts_axes.set_ylim(*ylim)

            # Gap
            x = ts_ratio['process']
            y = ts_ratio['median']
            yerr = (ts_ratio['median']-ts_ratio['lb'], ts_ratio['ub']-ts_ratio['median'])
            ratio_axes.errorbar(x, y, yerr, **ratio_kwargs)
            ratio_axes.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))

            # Gap threshold
            xlim = ratio_axes.get_xlim()
            ratio_axes.hlines(1.5, *xlim, color='k', label='Ratio = 1.5')
            ratio_axes.set_xlim(xlim)
            ratio_axes.set_yscale('log')

            # Adjust ylims
            ylim = ratio_axes.get_ylim()
            upper = int((np.log10(ylim[0])+1))+1
            new_upper = max(10**upper, ylim[1])
            ratio_axes.set_ylim(ylim[0], new_upper)
            ratio_axes.set_ylabel(r'Timescale ratio $t_{i}/t_{i+1}$')
            ratio_axes.set_xlabel(r'Process $i$')

            # Legends
            h, l = ts_axes.get_legend_handles_labels()
            n_leg = min(15, len(h))
            ts_axes.legend(h[:n_leg], l[:n_leg], bbox_to_anchor=(1, 1), loc='upper left')

            axes[1].legend(bbox_to_anchor=(1, 1), loc='upper left')

            for i, ax in enumerate(axes.flatten()):
                ax.annotate(text=f"({LETTERS[i]})", xy=(0.05, 0.9), xycoords='axes fraction')
            plt.tight_layout()

            plt.savefig(out_dir.joinpath(f'{protein_name}_{hp_ix}_its.pdf'), bbox_inches='tight')

    for hp_ix in hp_ixs:
        plot(hp_ix, timescales, gaps, lag, output_directory, name)

