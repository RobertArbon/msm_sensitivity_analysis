from pathlib import Path
from typing import List

import pandas as pd
import numpy as np
import seaborn as sns
import pyemma as pm
import h5py

import matplotlib.pyplot as plt
import matplotlib as mpl
from pandas.api.types import CategoricalDtype

from .plots import vamps_and_features


def select_lag(summary: Path, output_directory: Path, cutoff: float, proc: int = 2, units: str = 'ns') -> None:
    grad_cutoff = np.log(1+cutoff)

    # create dataframe for plotting and analysis
    ts_grad = pd.DataFrame(pd.read_hdf(summary, key='timescale_gradient'))  # read_hdf is typed as returning an 'object' which borks the annotations!
    ts_grad.dropna(inplace=True)
    ts_grad.reset_index(inplace=True)
    ts_grad = ts_grad.loc[ts_grad.process == proc, :]
    ts_grad.sort_values(by=['hp_ix', 'lag'], inplace=True)

    # calculate minimum lag that crosses cutoff
    ix = (ts_grad['median'] < grad_cutoff) & (ts_grad['median'] > 0)
    chosen_lag = ts_grad.loc[ix, :].groupby('hp_ix', as_index=False).first()
    chosen_lag = np.min(chosen_lag.lag)

    # Add to summary
    with h5py.File(summary, 'a') as f:
        g = f['summary']
        g.attrs['chosen_lag'] = chosen_lag
        name = g.attrs['protein_name']

    print(f'Chosen lag time for {name} is {chosen_lag} {units}')
    # plot
    with sns.plotting_context('paper', font_scale=1.25):
        fig, ax = plt.subplots(1, figsize=(6, 4))
        sns.lineplot(data=ts_grad, units='hp_ix', x='lag', y='median', estimator=None, alpha=0.5, ax=ax)
        # Threshold
        xlim = ax.get_xlim()
        ax.hlines(grad_cutoff, *xlim, color='k', label='Threshold', lw=2)
        ax.set_xlim()
        # chose lag
        ylim = ax.get_ylim()
        ax.vlines(chosen_lag, *ylim, color='g', label='Lag time', lw=2)
        ax.set_ylim(ylim)
        ax.legend(bbox_to_anchor=(1, 1), loc='upper left')

        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(20))
        ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(10))
        ax.grid()

        ax.set_ylabel(r'$\Delta \log{{(t_{{{0}}})}}$' + r'/$\Delta \tau$')
        ax.set_xlabel(f'lag time ({units})')
        plt.tight_layout()
        print(f"Saved image at {output_directory.joinpath(f'{name}_timescale_gradient.pdf')}")
        plt.savefig(output_directory.joinpath(f'{name}_timescale_gradient.pdf'), bbox_inches='tight')


def select_dominant(summary: Path, output_directory: Path, cutoff: int) -> None:
    SMALL_COUNT = 5
    MIN_PROCS = 10
    with h5py.File(summary, 'r') as f:
        grp = f['summary']
        try:
            chosen_lag = grp.attrs['chosen_lag']
        except KeyError:
            raise RuntimeError('No chosen lag - run `plot_timescale_gradient`')
        name = grp.attrs['protein_name']
    ts_ratio = pd.DataFrame(pd.read_hdf(summary, key='timescale_ratio'))
    ts_ratio.reset_index(inplace=True)
    ts_ratio.dropna(inplace=True)

    # Select lag
    ts_ratio = ts_ratio.loc[ts_ratio.lag == chosen_lag, :]
    # Rank gaps (we only want rank 1)
    ts_ratio['rank'] = ts_ratio.groupby(['hp_ix'], as_index=False)['median'].rank(ascending=False)

    # Some formatting for pretty plots
    cat_type = CategoricalDtype(categories=np.sort(ts_ratio.process.unique()).astype(int), ordered=True)
    ts_ratio['process'] = ts_ratio['process'].astype(int).astype(cat_type)

    # Select the process which is has the biggest gap
    top_gaps = ts_ratio.loc[ts_ratio['rank'] == 1, :]

    # Find the most common gap above threshold
    gap_counts = top_gaps.groupby('process', as_index=False).count()
    max_count = gap_counts.loc[gap_counts['rank'] > cutoff, 'process'].max()
    print(f'Chosen num dominant processes {name} is {max_count}')

    with h5py.File(summary, 'a') as f:
        grp = f['summary']
        grp.attrs['chosen_k'] = int(max_count)

    # Plot
    with sns.plotting_context('paper', font_scale=1.25):
        fig, ax = plt.subplots(1, figsize=(6, 4))
        sns.countplot(data=top_gaps, x='process', order=np.arange(2, max(max_count, MIN_PROCS) + 1).astype(int),
                      ax=ax, color=sns.color_palette('colorblind')[0])
        ax.set_ylabel('Count')
        ax.set_xlabel('Num. dominant processes')

        heights = [x.get_height() for x in ax.containers[0]]
        labels = [int(x) if x < SMALL_COUNT else '' for x in heights]
        ax.bar_label(ax.containers[0], labels)

        max_count_x = max_count - 2
        ax.annotate(text=f"*", xy=(max_count_x, heights[max_count_x]), xytext=(0, 1),
                    textcoords='offset points', ha='center')

        plt.ylim(0, 115)
        plt.tight_layout()
        print(f"Saved image at {output_directory.joinpath(f'{name}_timescale_gap.pdf')}")
        plt.savefig(output_directory.joinpath(f'{name}_timescale_gap.pdf'), bbox_inches='tight')


def select_model(summary: Path, hp_definitions: Path, cutoff: float) -> None:
    with h5py.File(summary, 'r') as f:
        grp = f['summary']
        lag = int(grp.attrs['chosen_lag'])
        k = int(grp.attrs['chosen_k'])
        name = grp.attrs['protein_name']

    results = dict()

    # Get VAMP scores
    vamps = vamps_and_features(summary, hp_definitions)
    vamps.reset_index(level=['lag'], inplace=True)
    vamps = vamps.loc[vamps.lag == lag, :]

    ####################################################################################################################
    # FIXED K
    ####################################################################################################################
    # Fix k and select
    vamps_k = vamps.reset_index(level=['process'])
    vamps_k = vamps_k.loc[vamps_k.process == k, :]
    vamps_k = vamps_k.drop(labels=['lag', 'process'], axis=1)

    vamps_k['rank'] = vamps_k['median'].rank(ascending=False)
    vamps_k.sort_values(by='rank', inplace=True, ascending=True)
    vamps_k['rank_by_feature'] = vamps_k.groupby(['feature'], as_index=False)['median'].rank(ascending=False)

    # Record best hps
    hps = list(vamps_k.loc[vamps_k.rank_by_feature == 1, :].index)
    results['fixed_k'] = list(zip([k]*len(hps), hps))

    # Record worst hps
    hps = list(vamps_k.loc[vamps_k['rank'] == vamps_k['rank'].max(), :].index)
    results['worst'] = list(zip([k]*len(hps), hps))


    ####################################################################################################################
    # TIMESCALE GAP
    ####################################################################################################################
    # Add gaps in
    gaps = pd.DataFrame(pd.read_hdf(summary, key='timescale_ratio'))
    gaps.reset_index(level=['lag'], inplace=True)
    gaps = gaps.loc[gaps.lag == lag, :]
    gaps = gaps.drop(labels=['lag'], axis=1)
    gaps.dropna(inplace=True)

    suffixes = ['_vamp', '_gap']
    vamps = vamps.merge(gaps, left_index=True, right_index=True, how='inner', suffixes=suffixes)
    vamps.reset_index(inplace=True)

    # Select vamps within threshold
    vamps['score_loss'] = 1-vamps['median_vamp']/vamps['process']
    top = vamps.loc[vamps['score_loss'] < cutoff, :].copy()

    # Select top gaps from remaining
    top['gap_per_feature_rank'] = top.groupby('feature')['median_gap'].rank(ascending=False)
    top.sort_values(by='score_loss', ascending=True, inplace=True)

    hps = top.loc[top['gap_per_feature_rank'] == 1.0, 'hp_ix'].values
    ks = top.loc[top['gap_per_feature_rank'] == 1.0, 'process'].values
    results['timescale_gap'] = list(zip(ks, hps))

    # Results dataframe
    vamps['vamp_rank'] = vamps.groupby('process')['median_vamp'].rank(ascending=False)
    results_df = []
    for method, idxs in results.items():
        for k, hp_ix in idxs:
            df = vamps.loc[(vamps.hp_ix == hp_ix) & (vamps.process == k), :].copy()
            df['method'] = method
            results_df.append(df)
    results_df = pd.concat(results_df)
    results_df.to_hdf(summary, key='model_selection')

    print(f'Selecting models based on fixed k ({k}) at lag {lag}:')
    print(results_df.loc[results_df['method'] == 'fixed_k', :])
    print()
    print(f"Selecting the worst model overall based on fixed k ({k}) and lag {lag}:")
    print(results_df.loc[results_df['method'] == 'worst', :])
    print()
    print(f"Selecting models based on the timescale gap:")
    print(results_df.loc[results_df['method'] == 'timescale_gap', :])
