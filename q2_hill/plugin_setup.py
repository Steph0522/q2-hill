from qiime2.plugin import (Plugin, Str, Properties, Choices, Int, Bool, Range,
                           Float, Set, Visualization, Metadata, MetadataColumn,
                           Categorical, Numeric, Citations, Threads)
from typing import Union
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.tree import Phylogeny, Rooted
from q2_types.sample_data import SampleData, AlphaDiversity
from q2_hill._methods import alpha_taxa, alpha_phylo, alpha_functional, compute_distance
from q2_types.distance_matrix import DistanceMatrix


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

plugin.methods.register_function(
    function=alpha_phylo,
    inputs={
        "table": FeatureTable[Frequency],
        "phylogeny": Phylogeny[Rooted],
    },
    parameters={
        "q": Float,
        "metric": Str % Choices(["PD", "qDT"])
    },
    outputs=[("alpha_diversity", SampleData[AlphaDiversity])],
    input_descriptions={
        "table": "Feature table with species abundances.",
        "phylogeny": "Rooted phylogenetic tree corresponding to species in the table."
    },
    parameter_descriptions={
        "q": "Order of Hill number (q ≥ 0).",
        'metric': "Metric to calculate: 'PD' for hillR o 'qDT' para hilldiv2. Default PD."
    },
    output_descriptions={
        "alpha_diversity": "Alpha diversity values per sample."
    },
    name="Phylogenetic Hill Diversity",
    description="Computes phylogenetic diversity using Hill numbers with a given order q."
)

plugin.methods.register_function(
    function=alpha_functional,
    inputs={
        "table": FeatureTable[Frequency]
    },
    parameters={
        "traits": Metadata,
        "q": Float,
        "dist": Str % Choices(["euclidean", "manhattan", "gower"]),
        "metric": Str % Choices(["FD", "FD_q", "D_q", "Q", "MD_q", "FDis"]),
        "tau": Float % Range(0, 1e10) | Str % Choices(["min", "max", "mean"])
    },
    outputs=[("diversity", SampleData[AlphaDiversity])],
    input_descriptions={
        "table": "Feature table with species abundances."
    },
    parameter_descriptions={
        "traits": "Metadata containing trait data with species as rows and trait variables as columns.",
        "q": "Order of Hill number (q ≥ 0).",
        "dist": "Distance metric to compute: 'euclidean', 'manhattan' or 'gower'.",
        "metric": "Metric to calculate: 'FD' for hilldiv2, or 'FD_q', 'D_q', 'Q', 'MD_q' or 'FDis' for hillR.",
        "tau": "Threshold for distances. Can be a positive number or None if no threshold is applied."
    },
    output_descriptions={
        "diversity": "Alpha diversity values per sample."
    },
    name="Functional Hill Diversity",
    description="Computes functional diversity (Hill numbers) using a feature table and traits metadata."
)

