# ----------------------------------------------------------------------------
# Copyright (c) 2024, Stephanie Hereira-Pacheco.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------# ----------------------------------------------------------------------------
# Copyright (c) 2024, Stephanie Hereira-Pacheco.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Plugin, Float
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.sample_data import SampleData, AlphaDiversity
import pandas as pd
from q2_hill._methods import alpha_taxa

plugin = Plugin(
    name="hill",
    version="0.0.1",
    website="https://github.com/tu-repositorio",
    package="q2_hill",
    description="Calculate Hill numbers from a Feature Table with Frequencies",
)

plugin.methods.register_function(
    function=alpha_taxa,
    inputs={"table": FeatureTable[Frequency]},
    parameters={"q": Float},
    outputs=[("alpha_diversity", SampleData[AlphaDiversity])],
    name="alpha_taxa",
    description="Calculate Hill numbers for a q value (order)",
    input_descriptions={"table": "The feature table containing the samples for which hill numbers should be computed."},
    parameter_descriptions={"q": "Order q of diversity (float between 0 and inf)"},
    output_descriptions={"alpha_diversity": "Hill numbers calculated for the order indicated"},
)

