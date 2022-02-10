from pathlib import Path

import pandas as pd
import numpy as np
import seaborn as sns
import pyemma as pm
import h5py

import matplotlib.pyplot as plt
import matplotlib as mpl


def timescale_gradient(summary: Path, output_directory: Path, cutoff: float, proc: int = 2, units: str = 'ns') -> None:
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


    