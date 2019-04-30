from .plugin_setup import plugin
from ._semantics import TSVMolecules
import pandas as pd


def _read_dataframe(fh):
    # Using `dtype=object` and `set_index` to avoid type casting/inference
    # of any columns or the index.
    df = pd.read_csv(fh, sep='\t', header=0, dtype='str')
    df.set_index(df.columns[0], drop=True, append=False, inplace=True)
    df.index.name = 'id'
    return df

# define a transformer from pd.DataFrame to -> TSVMolecules
@plugin.register_transformer
def _1(data: pd.DataFrame) -> TSVMolecules:
    ff = TSVMolecules()
    with ff.open() as fh:
        data.to_csv(fh, sep='\t', header=True)
    return ff

# define a transformer from TSVMolecules -> pd.DataFrame
@plugin.register_transformer
def _2(ff: TSVMolecules) -> pd.DataFrame:
    with ff.open() as fh:
        df = _read_dataframe(fh)
        return df

from ._semantics import FingerprintNetworkEdgesFile
import pandas as pd

@plugin.register_transformer
def dataframe_to_edges(data: pd.DataFrame) -> FingerprintNetworkEdgesFile:
    ff = FingerprintNetworkEdgesFile()
    with ff.open() as fh:
        data.to_csv(path_or_buf=fh, sep="\t", index=False)
    return ff