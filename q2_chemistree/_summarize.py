# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage

from ._collate_fingerprint import collate_fingerprint
from ._semantics import CSIDirFmt

def output_fingerprints(csi_results: CSIDirFmt) -> pd.DataFrame:
    all_fingerprints = collate_fingerprint(csi_results, False)
    all_fingerprints["#FeatureID"] = all_fingerprints.index
    print(all_fingerprints.head())
    return all_fingerprints