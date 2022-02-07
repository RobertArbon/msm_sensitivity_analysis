from typing import List, Dict
import pickle
from pathlib import Path

import pandas as pd
import numpy as np


def extract_result(results: Dict, key: str) -> pd.Series:
    df = pd.concat({(int(res['hp_ix']), int(bs_ix), k): pd.DataFrame(v, index=[0]).T for bs_ix, res in results.items()
                       for k, v in res[key].items()})
    df.index.rename(('hp_ix', 'bs_ix', 'lag', 'process'), inplace=True)
    s = df[0]
    return s


def gradient(s: pd.Series, denominator: str) -> pd.Series:
    controls = list(s.index.names)
    controls.remove(denominator)
    index_names = controls + [denominator]
    s = s.sort_index(level=index_names, inplace=False)
    dy = s.groupby(controls).diff()
    dx = s.index.get_level_values(denominator).to_series().diff().values
    grad = pd.Series(dy/dx)
    return grad


def collate_results(hp_directory: Path) -> Dict[str, pd.Series]:
    bs_paths = hp_directory.glob('*.pkl')
    results = {x.stem: pickle.load(x.open('rb')) for x in bs_paths}
    vamps = extract_result(results, key='vamp')
    ts = extract_result(results, key='ts')
    ts_proc_grad = np.exp(gradient(np.log(ts), denominator='process'))
    ts_lag_grad = gradient(np.log(ts), denominator='lag')
    results = dict(
        vamps=vamps,
        timescales=ts,
        timescale_ratio=ts_proc_grad,
        timescale_gradient=ts_lag_grad
    )
    return results


def series_summary(s: pd.Series) -> pd.DataFrame:
    controls = list(s.index.names)
    controls.remove('bs_ix')
    df = s.groupby(controls).agg(
            median = lambda x: np.quantile(x, 0.5),
            lb = lambda x: np.quantile(x, 0.025),
            ub = lambda x: np.quantile(x, 0.975)
    )
    return df


def summarise_results(results: Dict[str, pd.Series]) -> Dict[str, pd.DataFrame]:
    summary = {k: series_summary(v) for k, v in results.items()}
    return summary


def write_summaries(out_dir: Path, summaries: Dict[str, pd.DataFrame]) -> None:
    for k, v in summaries.items():
        v.to_hdf(str(out_dir.joinpath("summary.h5")), key=k)


def main(directory: Path) -> None:
    hp_dirs = [x for x in directory.glob('*') if x.is_dir()]
    all_summaries = dict()

    for hp_dir in hp_dirs:
        results = collate_results(hp_dir)
        summaries = summarise_results(results)
        for k, v in summaries.items():
            try:
                all_summaries[k] = pd.concat([all_summaries[k], v])
            except KeyError:
                all_summaries[k] = v

    write_summaries(directory, all_summaries)






