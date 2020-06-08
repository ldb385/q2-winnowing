import qiime2.plugin
from q2_types.feature_table import FeatureTable, RelativeFrequency
from q2_types.tree import Phylogeny, Rooted, Unrooted

import q2_winnowing
from q2_winnowing.winnow_data import winnow_data

# cites = qiime2.plugin.Citations.load("citations.bib", package="q2_winnowing")
# Note: this can be replaced with a bibliography when the thesis is completed.

# Setup Choice lists for easy accesability
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
# Currently all code is stored in a main pipeline, therefore only one function will be registered with many parameters.

# main
plugin.methods.register_function(
    function=winnow_data,
    inputs={
        "biomFile1": FeatureTable[RelativeFrequency]
    },
    outputs=[
        # TODO: Verify this is the proper output
        ("metric_original_values_result", FeatureTable[RelativeFrequency]),
        ("abundance_values_result", FeatureTable[RelativeFrequency]),
        ("graph_network_visual_result", Phylogeny[Unrooted]),
        ("metric_network_values_result", FeatureTable[RelativeFrequency]),
        ("metric_values_result", Phylogeny[Unrooted]),
        ("metric_network_visual_result", Phylogeny[Unrooted]),
        ("parameter_list_result", FeatureTable[RelativeFrequency])
    ],
    input_descriptions={
        "biomFile1": ("This is the path to the csv file of OTU's that the data is read from."
                      "First file read in.")
    },
    parameters={
        "ab_comp": qiime2.plugin.Bool,
        "biomFile2": qiime2.plugin.Str,
        "metric_name": qiime2.plugin.Str % qiime2.plugin.Choices(_METRIC_TYPES_),
        "evaluation_type": qiime2.plugin.Str % qiime2.plugin.Choices(_EVALUATION_TYPES_),
        "min_count": qiime2.plugin.Int,
        "c_type": qiime2.plugin.Str % qiime2.plugin.Choices(_CONDITIONING_TYPES_),
        "total_select": qiime2.plugin.Metadata,
        "iteration_select": qiime2.plugin.Metadata,
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
        "min_connected": qiime2.plugin.Float
    },
    parameter_descriptions={
        "ab_comp": ("Boolean representing whether to perform AB comparison on the data."
                    "Possible inputs are:"
                    f"\t {_BOOLEAN_}"),
        "biomFile2": ("Second file to user for A/B comparison."),
        "metric_name": ("This is the metric to use."
                        "Possible metrics include: "
                        f"\t {_METRIC_TYPES_}"),
        "evaluation_type": ("This is the evaluation type to use."
                            "Possible evaluation types are:"
                            f"\t {_EVALUATION_TYPES_}"),
        "min_count": ("Features with counts below this number will be removed."),
        "c_type": ("Conditioning type to use on the data."
                   "Possible conditioning types are:"
                   f"\t {_CONDITIONING_TYPES_}"),
        "total_select": ("Number of features to select in total."
                         "Possible selections are:"
                         f"\t {_ALL_OR_INT_}"),
        "iteration_select": ("Number of features to select for each time the metric is called."
                             "Possible selections are:"
                             f"\t {_ALL_OR_INT_}"),
        "pca_components": ("Number of pca components to find"),
        "smooth_type": ("Type of Smoothing to be used to remove noise."
                        "Possible smoothing:"
                        f"\t {_SMOOTHING_TYPES_}"),
        "window_size": ("If Smoothing type is a sliding window, this is the size of the window."),
        "centrality_type": ("If graph_centrality is the metric type, this is the type of Centrality to use."
                            "Possible centrality types include:"
                            f"\t {_CENTRALITY_TYPES_}"),
        "keep_threshold": ( "If graph_centrality is the metric type, this is the threshold to use to remove weak edges."),
        "correlation":
            ("If graph_centrality is the metric type, this is the type of correlation to use to build the graph."
             "Possible correlations are:"
             f"\t {_CORRELATION_TYPES_}"),
        "weighted": (
            "If graph_centrality is the metric type, this specifies if weighted edges should be used to create the graph."
            "Possible inputs are:"
            f"\t {_BOOLEAN_}"),
        "corr_prop": (
            "If graph centrality is the metric, this specifies if positive, negative, or both types of correlation should be used."
            "Possible correlation properties are:"
            f"\t {_CORRELATION_PROPERTIES_}"),
        "plot_metric": (
            "Including this parameter will create a line plot of the metric values for the selected features."),
        "create_graph": ("If graph centrality is the metric, including this parameter will create a graph image of "
                         "the selected features (using the same correlation type used to select the features)."
                         "Possible inputs are:"
                         f"\t {_BOOLEAN_}"),
        "plot_pca": ("If PCA is the metric, including this parameter will create a scatter plot image of "
                     "the first two principle components."
                     "Possible inputs are:"
                     f"\t {_BOOLEAN_}"),
        "naming_file": (
            "The file to be used to name the features. If not used, the features will be outputted with the names the input file."),
        "proc_id": ("The identifying number to use in the output file names."),
        "min_connected": (
            "The minimum percentage of connectedness of the graph that should be considered before the winnowing process is aborted.")
    },
    # TODO: fill in actual return descriptions
    output_descriptions={
        "metric_original_values_result": ("A"),
        "abundance_values_result": ("B"),
        "graph_network_visual_result": ("C"),
        "metric_network_values_result": ("D"),
        "metric_values_result": ("E"),
        "metric_network_visual_result": ("F"),
        "parameter_list_result": ("G"),
    },
    name='Winnowing',
    description=("Infer the interaction type of microbial communities through statistical analysis. "
                 "This will allow for a better understanding of taxa interaction at a micro scale.")
)
