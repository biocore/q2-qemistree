from .plugin_setup import plugin
from ._semantics import FingerprintNetworkEdgesFile
import pandas as pd

@plugin.register_transformer
def dataframe_to_edges(data: pd.DataFrame) -> FingerprintNetworkEdgesFile:
    ff = FingerprintNetworkEdgesFile()
    with ff.open() as fh:
        data.to_csv(path_or_buf=fh, sep="\t", index=False)
    return ff
