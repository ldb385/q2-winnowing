

import qiime2.plugin
from q2_types.feature_table import FeatureTable, RelativeFrequency, BIOMV210DirFmt
from q2_types.tree import Phylogeny, Rooted, Unrooted

import q2_winnowing
from q2_winnowing.winnow import winnow_processing

# cites = qiime2.plugin.Citations.load("citations.bib", package="q2_winnowing")
# Note: this can be replaced with a bibliography when the thesis is completed.

# Setup Choice lists for easy accesability
# <><><> NOTE: THESE ARE FOR STEPS 1-3 <><><>
_METRIC_TYPES_ = ["graph_centrality", "pca_importance"]
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

# <><><> Register functions <><><>

# pipeline: step 1-3
plugin.methods.register_function(
    name='processing',
    description=("Infer the interaction type of microbial communities through statistical analysis. "
                 "This will allow for a better understanding of taxa interaction at a micro scale."),
    function=winnow_processing,
    inputs={
        "infile1": FeatureTable[RelativeFrequency],
        "infile2": FeatureTable[RelativeFrequency]
    },
    outputs=[
        # TODO: Verify this is the proper output
        ("interaction_of_taxa_result", FeatureTable[RelativeFrequency])
    ],
    input_descriptions={
        "infile1": ("This is the biom file which will have OTU info extracted from and analyzed to generate "
                    "an interaction table of taxom."),
        "infile2": ("This is only used in the case of an A/B analysis and will not be used if ab_comp is False.")
    },
    parameters={
        "name": qiime2.plugin.Str,
        "ab_comp": qiime2.plugin.Bool,
        "metric_name": qiime2.plugin.Str % qiime2.plugin.Choices(_METRIC_TYPES_),
        "evaluation_type": qiime2.plugin.Str % qiime2.plugin.Choices(_EVALUATION_TYPES_),
        "min_count": qiime2.plugin.Int,
        "c_type": qiime2.plugin.Str % qiime2.plugin.Choices(_CONDITIONING_TYPES_),
        "total_select": qiime2.plugin.Str,
        "iteration_select": qiime2.plugin.Str,
        "pca_components": qiime2.plugin.Int,
        "smooth_type": qiime2.plugin.Str % qiime2.plugin.Choices(_SMOOTHING_TYPES_),
        "window_size": qiime2.plugin.Int,
        "centrality_type": qiime2.plugin.Str % qiime2.plugin.Choices(_CENTRALITY_TYPES_),
        "keep_threshold": qiime2.plugin.Float,
        "correlation": qiime2.plugin.Str % qiime2.plugin.Choices(_CORRELATION_TYPES_),
        "weighted": qiime2.plugin.Bool,
        "corr_prop": qiime2.plugin.Str % qiime2.plugin.Choices(_CORRELATION_PROPERTIES_),
        "plot_metric": qiime2.plugin.Bool,
        "create_graph": qiime2.plugin.Bool,
        "plot_pca": qiime2.plugin.Bool,
        "naming_file": qiime2.plugin.Str,
        "proc_id": qiime2.plugin.Int,
        "min_connected": qiime2.plugin.Float,
        "detailed": qiime2.plugin.Bool,
        "verbose": qiime2.plugin.Bool
    },
    parameter_descriptions={
        "name": ("This is the string that will be attached to output files. This is used especially in "
                 "the case of detailed."),
        "ab_comp": ("Boolean representing whether to perform AB comparison on the data."),
        "metric_name": ("This is the metric to use."),
        "evaluation_type": ("This is the evaluation type to use."),
        "min_count": ("Features with counts below this number will be removed."),
        "c_type": ("Conditioning type to use on the data."),
        "total_select": ("Number of features to select in total."
                         "Possible selections are:"
                         f"\t {_ALL_OR_INT_}"),
        "iteration_select": ("Number of features to select for each time the metric is called."
                             "Possible selections are:"
                             f"\t {_ALL_OR_INT_}"),
        "pca_components": ("Number of pca components to find"),
        "smooth_type": ("Type of Smoothing to be used to remove noise."),
        "window_size": ("If Smoothing type is a sliding window, this is the size of the window."),
        "centrality_type": ("If graph_centrality is the metric type, this is the type of Centrality to use."),
        "keep_threshold": ( "If graph_centrality is the metric type, this is the threshold to use to remove weak edges."),
        "correlation":
            ("If graph_centrality is the metric type, this is the type of correlation to use to build the graph."),
        "weighted": (
            "If graph_centrality is the metric type, this specifies if weighted edges should be used to create the graph."),
        "corr_prop": (
            "If graph centrality is the metric, this specifies if positive, negative, or both types of correlation should be used."),
        "plot_metric": (
            "Including this parameter will create a line plot of the metric values for the selected features."),
        "create_graph": ("If graph centrality is the metric, including this parameter will create a graph image of "
                         "the selected features (using the same correlation type used to select the features)."),
        "plot_pca": ("If PCA is the metric, including this parameter will create a scatter plot image of "
                     "the first two principle components."),
        "naming_file": (
            "The file to be used to name the features. If not used, the features will be outputted with the names the input file."),
        "proc_id": ("The identifying number to use in the output file names."),
        "min_connected": (
            "The minimum percentage of connectedness of the graph that should be considered before the winnowing process is aborted."),
        "detailed": ("Notifies plugin to output diagrams and csv files to each steps respective output folder throughout"
                     "computation. If not enabled files will not be generated"),
        "verbose": ("Notifies plugin to generate dump files for every step. These will contain all data that previously "
                    "may have been output through print statements during execution. Each dump.txt file is stored in "
                    "output foleders that correspond with each step.")
    },
    output_descriptions={
        "interaction_of_taxa_result": ("This is the completed table of different taxa with their corresponding interactions"
                                       "to other taxa in their enviroment. This was computed through the use of methods"
                                       "such as: Abundance analysis, AUC analysis, F-Score ordering,"
                                       "PERMANOVA calculation, Jacobian matrices, and SEM analysis ")
    }
)


