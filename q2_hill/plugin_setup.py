import qiime2.plugin
from qiime2.plugin import Plugin, Str, Choices, Float, Range, Metadata
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency, PresenceAbsence
from q2_types.tree import Phylogeny, Rooted
from q2_types.sample_data import SampleData, AlphaDiversity
from q2_types.distance_matrix import DistanceMatrix
from q2_hill._methods import (
    alpha_taxa,
    alpha_phylo,
    alpha_functional,
    beta_taxa
)

plugin = Plugin(
    name="hill",
    version="0.0.1",
    website="https://github.com/tu-repositorio",
    package="q2_hill",
    description="Compute Hill numbers of different types - taxonomic, phylogenetic, and functional - for alpha and beta diversity metrics at different orders of q",
)

plugin.methods.register_function(
    function=alpha_taxa,
    inputs={'table':
            FeatureTable[Frequency | RelativeFrequency | PresenceAbsence]},
    parameters={"q": Float},
    outputs=[("alpha_diversity", SampleData[AlphaDiversity])],
    name="alpha_taxa",
    description="Calculate Hill taxonomic numbers for a q value (order)",
    input_descriptions={
        "table": "The feature table containing the samples for which Hill "
        "numbers should be computed."
    },
    parameter_descriptions={
        "q": "Order q of diversity (float between 0 and inf)"
    },
    output_descriptions={
        "alpha_diversity": "Taxonomic Hill numbers calculated for the order indicated"
    },
)

plugin.methods.register_function(
    function=alpha_phylo,
    inputs={'table':
            FeatureTable[Frequency | RelativeFrequency | PresenceAbsence],
        "phylogeny": Phylogeny[Rooted],
    },
    parameters={
        "q": Float,
        "metric": Str % Choices(["PD", "qDT"]),
    },
    outputs=[("alpha_diversity", SampleData[AlphaDiversity])],
    input_descriptions={
        "table": "The feature table containing the samples for which Hill "
        "numbers should be computed.",
        "phylogeny": (
            "Rooted phylogenetic tree corresponding to species in the table."
        ),
    },
    parameter_descriptions={
        "q": "Order q of diversity (float between 0 and inf)",
        "metric": (
            "Metric to calculate: 'PD' as hillR (effective total branch length) or 'qDT' as hilldiv2 (effective number of lineages). "
            "Default PD."
        ),
    },
    output_descriptions={
        "alpha_diversity": "Phylogenetic Hill numbers calculated for the order indicated."
    },
    name="Phylogenetic Hill Diversity",
    description=(
        "Computes phylogenetic diversity using Hill numbers with a "
        "given order q."
    ),
)

plugin.methods.register_function(
    function=alpha_functional,
    inputs={"table": FeatureTable[Frequency]},
    parameters={
        "traits": Metadata,
        "q": Float,
        "dist": Str % Choices(["euclidean", "manhattan", "gower"]),
        "metric": Str % Choices(["FD", "FD_q", "D_q", "Q", "MD_q", "FDis"]),
        "tau": (Float % Range(0, 1e10))
        | (Str % Choices(["min", "max", "mean"])),
    },
    outputs=[("diversity", SampleData[AlphaDiversity])],
    input_descriptions={"table": "The feature table containing the samples for which Hill "
        "numbers should be computed."},
    parameter_descriptions={
        "traits": (
            "Metadata containing trait data with species as rows and "
            "trait variables as columns."
        ),
        "q": "Order of Hill number (q ≥ 0).",
        "dist": (
            "Distance metric to compute: 'euclidean', 'manhattan' or "
            "'gower'."
        ),
        "metric": (
            "Metric to calculate: 'FD' for hilldiv2, or 'FD_q', 'D_q', "
            "'Q', 'MD_q' or 'FDis' for hillR."
        ),
        "tau": (
            "Threshold for distances. Can be a positive number or one of "
            "'min', 'max', or 'mean'."
        ),
    },
    output_descriptions={"diversity": "Alpha diversity values per sample."},
    name="Functional Hill Diversity",
    description=(
        "Computes functional diversity (Hill numbers) using a feature "
        "table and traits metadata."
    ),
)
plugin.methods.register_function(
    function=beta_taxa,
    inputs={"table": FeatureTable[Frequency | RelativeFrequency | PresenceAbsence]},
    parameters={"q": Float,
                "metric": Str % Choices(["C", "S", "V", "U"]),
                "similarity": qiime2.plugin.Bool},
    outputs=[("distance_matrix", DistanceMatrix)],
    input_descriptions={"table": "The feature table containing the samples for which Hill "
        "numbers should be computed. "},
    parameter_descriptions={
        "q": "Order of diversity (float between 0 and inf).",
        "metric": "Metric to calculate: 'C':Sørensen-type overlap-complement, 'S':Jaccard-type turnover, 'V':Sørensen-type turnover or 'U':Jaccard-type overlap-complement.",
    },
    output_descriptions={
        "distance_matrix": "Pairwise distance matrix based on Hill numbers for the metric chosen"
    },
    name="hillpair_taxa",
    description="Computes overall dissimilarity metrics based on the Hill taxonomic numbers beta diversity following Chiu et al. (2014).  "
)
