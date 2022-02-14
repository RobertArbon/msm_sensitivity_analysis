from typing import List, Dict, Union
import pickle
from pathlib import Path

import pandas as pd
import numpy as np
import h5py


def extract_result(results: Dict, key: str) -> pd.Series:
    df = pd.concat({(int(res['hp_ix']), int(bs_ix), k): pd.DataFrame(v, index=[0]).T for bs_ix, res in results.items()
                       for k, v in res[key].items()})
    df.index.rename(('hp_ix', 'bs_ix', 'lag', 'process'), inplace=True)
    s = df[0]
    return s


def gradient(s: pd.Series, denominator: str, periods: int = 1) -> pd.Series:
    controls = list(s.index.names)
    controls.remove(denominator)
    index_names = controls + [denominator]
    s = s.sort_index(level=index_names, inplace=False)
    dy = s.groupby(controls).diff(periods=periods)
    dx = s.index.get_level_values(denominator).to_series().diff(periods=periods).values
    grad = pd.Series(dy/dx)
    return grad


def collate_results(hp_directory: Path) -> Dict[str, pd.Series]:
    bs_paths = hp_directory.glob('*.pkl')
    results = {x.stem: pickle.load(x.open('rb')) for x in bs_paths}
    vamps = extract_result(results, key='vamp')
    ts = extract_result(results, key='ts')
    # The periods variable aligns the ratio with the correct process.
    # The (-1) gets the ratio the correct way round (ts_2/ts_3) as the periods affects both the demon and numer.
    ts_proc_grad = np.exp((-1)*gradient(np.log(ts), denominator='process', periods=-1))
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
            median=lambda x: np.quantile(x, 0.5),
            lb=lambda x: np.quantile(x, 0.025),
            ub=lambda x: np.quantile(x, 0.975),
            count=lambda x: x.shape[0]-x.isna().sum()

    )
    return df


def summarise_results(results: Dict[str, pd.Series]) -> Dict[str, pd.DataFrame]:
    summary = {k: series_summary(v) for k, v in results.items()}
    return summary


def write(out_dir: Path, summaries: Dict[str, Union[pd.DataFrame, pd.Series]], protein_name: str, protein_label: str,
          file_name: str = 'summary') -> None:
    for k, v in summaries.items():
        v.to_hdf(str(out_dir.joinpath(f"{file_name}.h5")), key=k)
    with h5py.File(out_dir.joinpath(f"{file_name}.h5"), 'a') as f:
        grp = f.create_group(f"{file_name}")
        grp.attrs['protein_name'] = protein_name
        grp.attrs['protein_label'] = protein_label


def main(directory: Path, dump_raw: bool) -> None:
    label = directory.stem
    names_by_labels = {'1fme': 'BBA'}

    hp_dirs = [x for x in directory.glob('*') if x.is_dir()]
    all_summaries = dict()
    all_raw = dict()
    for hp_dir in hp_dirs:
        results = collate_results(hp_dir)
        if dump_raw:
            for k, v in results.items():
                try:
                    all_raw[k] = pd.concat([all_raw[k], v])
                except KeyError:
                    all_raw[k] = v

        summaries = summarise_results(results)
        for k, v in summaries.items():
            try:
                all_summaries[k] = pd.concat([all_summaries[k], v])
            except KeyError:
                all_summaries[k] = v
    if dump_raw:
        print("Dumping raw results")
        write(directory, all_raw, protein_name=names_by_labels[label], protein_label=label, file_name='raw')

    write(directory, all_summaries, protein_name=names_by_labels[label], protein_label=label)






