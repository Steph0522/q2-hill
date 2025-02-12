# ----------------------------------------------------------------------------
# Copyright (c) 2024, Stephanie Hereira-Pacheco.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import pandas.testing as pdt

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform
from q2_types.feature_table import BIOMV100Format

from q2_hill._methods import hill_numbers



