# ----------------------------------------------------------------------------
# Copyright (c) 2024, Stephanie Hereira-Pacheco.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd

def alpha_taxa(table: pd.DataFrame, q: float) -> pd.Series:
    # check if is a df
    if isinstance(table, np.ndarray):
        comm = pd.DataFrame(table)

    # abundances to props
    prop = table.div(table.sum(axis=1), axis=0)

    # Calculation of Hill numbers
    if q == 0:
        hill_numbers = (table > 0).sum(axis=1)  # Richness
    elif q == 1:
        hill_numbers = np.exp(-np.sum(prop * np.log(prop.replace(0, np.nan)), axis=1))  # Shannon
    else:
        hill_numbers = (prop**q).sum(axis=1) ** (1 / (1 - q))  # General formula

    # Return as serie
    return pd.Series(hill_numbers, index=table.index, name=f"q={q}")

