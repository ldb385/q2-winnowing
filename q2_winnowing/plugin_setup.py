
import importlib

import qiime2.plugin
from qiime2.plugin import MetadataColumn, Categorical
from q2_types.feature_table import FeatureTable, RelativeFrequency

import q2_winnowing
from ._type import Winnowed
from ._format import ( WinnowedDirectoryFormat, WinnowedFeatureOrderingFormat,
                       WinnowedAucOrderingFormat, WinnowedPermanovaOrderingFormat )
from q2_winnowing.winnow import processing
from q2_winnowing._summarize._visualizer import summarize

# cites = qiime2.plugin.Citations.load("citations.bib", package="q2_winnowing")
# Note: this can be replaced with a bibliography when the thesis is completed.

# Setup Choice lists for easy accessibility
# <><><> NOTE: THESE ARE FOR STEPS 1-3 <><><>
_METRIC_TYPES_ = ["graph_centrality", "pca_importance", "log_transform", "abundance"]
_CONDITIONING_TYPES_ = ["add_one", "hellinger"]
_EVALUATION_TYPES_ = ["kl_divergence"]
_CENTRALITY_TYPES_ = ["betweenness", "closeness", "degree", "eigenvector"]
_SMOOTHING_TYPES_ = ["sliding_window"]
_CORRELATION_TYPES_ = ["spearman", "pearson", "kendall", "MIC"]
_CORRELATION_PROPERTIES_ = ["negative", "positive", "both"]
_ALL_OR_INT_ = ["all", "0,1,2,3,..."]
_BOOLEAN_ = ["True", "False"]


plugin = qiime2.plugin.Plugin(
    name="winnowing",
    version=q2_winnowing.__version__,
    website="https://github.com/ldb385/q2-winnowing",
    package="q2_winnowing",
    description="A Qiime2 plugin that will infer the interaction type of microbial communities through statistical analysis. " +
                "This will allow for a better understanding of taxa interaction at a micro scale.",
    short_description="Plugin for inferring the interaction type of microbial communities"
)


# <><><> Register semantic types and formats <><><>
plugin.register_semantic_types(Winnowed)

# Register the formats defined
plugin.register_formats( WinnowedFeatureOrderingFormat, WinnowedAucOrderingFormat,
                         WinnowedPermanovaOrderingFormat, WinnowedDirectoryFormat )

# Register directory format with semantic type
plugin.register_semantic_type_to_format(
    Winnowed,
    artifact_format=WinnowedDirectoryFormat
)


# <><><> Register functions <><><>

plugin.methods.register_function(
    name='winnowing processing',
    description=("Infer the interaction type of microbial communities through statistical analysis. "
                 "This will allow for a better understanding of taxa interaction at a micro scale."),
    function=processing,
    inputs={
        "infile1": FeatureTable[RelativeFrequency],
        "infile2": FeatureTable[RelativeFrequency]
    },
    outputs=[
        ("result", Winnowed )
    ],
    input_descriptions={
        "infile1": ("This is the biom file which will have OTU info extracted from and analyzed to generate "
                    "an interaction table of taxom."),
        "infile2": ("This is only used in the case of an A/B analysis and will not be used if ab_comp is False.")
    },
    parameters={
        "name": qiime2.plugin.Str,
        "sample_types": MetadataColumn[Categorical],
        "ab_comp": qiime2.plugin.Bool,
        "metric": qiime2.plugin.Str % qiime2.plugin.Choices(_METRIC_TYPES_),
        "evaluation": qiime2.plugin.Str % qiime2.plugin.Choices(_EVALUATION_TYPES_),
        "min_count": qiime2.plugin.Int,
        "conditioning": qiime2.plugin.Str % qiime2.plugin.Choices(_CONDITIONING_TYPES_),
        "total_select": qiime2.plugin.Str,
        "iteration_select": qiime2.plugin.Set[ qiime2.plugin.Int ],
        "pca_components": qiime2.plugin.Int,
        "smooth_type": qiime2.plugin.Str % qiime2.plugin.Choices(_SMOOTHING_TYPES_),
        "window_size": qiime2.plugin.Int,
        "centrality": qiime2.plugin.Str % qiime2.plugin.Choices(_CENTRALITY_TYPES_),
        "keep_threshold": qiime2.plugin.Float,
        "correlation": qiime2.plugin.Str % qiime2.plugin.Choices(_CORRELATION_TYPES_),
        "weighted": qiime2.plugin.Bool,
        "correlation_prop": qiime2.plugin.Str % qiime2.plugin.Choices(_CORRELATION_PROPERTIES_),
        "min_connected": qiime2.plugin.Float,
        "detailed": qiime2.plugin.Bool,
        "verbose": qiime2.plugin.Bool
    },
    parameter_descriptions={
        "name": ("This is the string that will be attached to output files. This is used especially in "
                 "the case of detailed. "),
        "sample_types": ("This data provides a legend for which samples are of different types. This allows for "
                         "permutations in PERMANOVA calculations to be descriptive. Data is labelled as "
                         "invaded and natural "),
        "ab_comp": ("Boolean representing whether to perform AB comparison on the data. "),
        "metric": ("This is the metric to use. "),
        "evaluation": ("This is the evaluation type to use. "),
        "min_count": ("Features with counts below this number will be removed. "),
        "conditioning": ("Conditioning type to use on the data. "),
        "total_select": ("Number of features to select in total. "
                         "Possible selections are: "
                         f"{_ALL_OR_INT_}"),
        "iteration_select": ("Number of features to select for each time the metric is called. "
                             "Note: a set of values must be given to be later used in a kappa calcultaion "
                             "An example of an input could be [1, 4, 16, 64, 128] "),
        "pca_components": ("Number of pca components to find. "),
        "smooth_type": ("Type of Smoothing to be used to remove noise. "),
        "window_size": ("If Smoothing type is a sliding window, this is the size of the window. "),
        "centrality": ("If graph_centrality is the metric type, this is the type of Centrality to use. "),
        "keep_threshold": ( "If graph_centrality is the metric type, this is the threshold to use to remove weak edges."),
        "correlation":
            ("If graph_centrality is the metric type, this is the type of correlation to use to build the graph. "),
        "weighted": (
            "If graph_centrality is the metric type, this specifies if weighted edges should be used to create the graph. "),
        "correlation_prop": (
            "If graph centrality is the metric, this specifies if positive, negative, or both types of correlation should be used. "),
        "min_connected": (
            "The minimum percentage of connectedness of the graph that should be considered before the winnowing process is aborted. "),
        "detailed": ("Notifies plugin to output diagrams and csv files to each steps respective output folder throughout "
                     "computation. If not enabled files will not be generated. "),
        "verbose": ("Notifies plugin to generate dump files for every step. These will contain all data that previously "
                    "may have been output through print statements during execution. Each dump.txt file is stored in "
                    "output foleders that correspond with each step. ")
    },
    output_descriptions={
        "result": ("This is a directory containing an feature ordering based off influential taxom. "
                                       "Feature ordering also includes kappa and agreement values of iteration selection. "
                                       "AUC and Permanova values from the highest iteration selection. ")
    }
)


# <><><> Register visualizers <><><>
plugin.visualizers.register_function(
    function=summarize,
    inputs={ "data": Winnowed },
    parameters={},
    input_descriptions={
        "data": ( "this is the measure OTU's along with their interactions and the meausure of interaction" )
    },
    parameter_descriptions={
    },
    name="Winnow Interaction Visualization",
    description=("Visualize the output from processing function through graphs.")
)




importlib.import_module('q2_winnowing._transformer') # Avoid circular dependencies