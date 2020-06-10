

import biom

from qiime2.plugin import Bool, Str, Int, Float

from q2_winnowing._file_conversions import *
from q2_winnowing.step1_3.pipeline import main



def winnow_pipeline( biomFile1: biom.Table , ab_comp: Bool=False, biomFile2: biom.Table=None, metric_name: Str=None,
                 c_type: Str=None, min_count: Int=3, total_select: Str=None, iteration_select: Str=None,
                 pca_components: Int=4, smooth_type: Str="sliding_window", window_size: Int=3, centrality_type: Str=None,
                 keep_threshold: Float=0.5, correlation: Str=None, weighted: Bool=False, corr_prop: Str="both",
                 evaluation_type: Str=None, plot_metric: Bool=False, create_graph: Bool=False, plot_pca: Bool=False,
                 naming_file: Str=None, proc_id: Int=0, min_connected: Int=0
                 ) -> ( biom.Table, biom.Table, None, biom.Table, None, None, biom.Table ):
    # TODO: Implement proper return types
    """
    Note this function executes the main functionality of steps 1-3 in the pipeline of
    winnowing data.

    :param biomFile1:
    :param ab_comp:
    :param biomFile2:
    :param metric_name:
    :param c_type:
    :param min_count:
    :param total_select:
    :param iteration_select:
    :param pca_components:
    :param smooth_type:
    :param window_size:
    :param centrality_type:
    :param keep_threshold:
    :param correlation:
    :param weighted:
    :param corr_prop:
    :param evaluation_type:
    :param plot_metric:
    :param create_graph:
    :param plot_pca:
    :param naming_file:
    :param proc_id:
    :param min_connected:
    :return:
    """

    csvFile1 = biom_to_csv( biomFile1 )
    csvFile2 = None

    if( ab_comp ):
        csvFile2 = biom_to_csv( biomFile2 )

    output = main( csvFile1, csvFile2, ab_comp, metric_name, c_type, min_count,
                 total_select, iteration_select, pca_components, smooth_type,
                 window_size, centrality_type, keep_threshold, correlation,
                 weighted, corr_prop, evaluation_type, plot_metric,
                 create_graph, plot_pca, naming_file, proc_id, min_connected
                 )

    # TODO: return values must coorelate with their respective documents
    # these are what is supposed to be output
    metric_original_values_result = None
    abundance_values_result = None
    graph_network_visual_result = None
    metric_network_values_result = None
    metric_values_result = None
    metric_network_visual_result = None
    parameter_list_result = None

    # TODO: This can be put in an --output-dir
    output = (
        metric_original_values_result, abundance_values_result, graph_network_visual_result,
        metric_network_values_result, metric_values_result, metric_network_visual_result, parameter_list_result
    )


    return output



def winnow_ordering():


    return


def winnow_permanova():


    return



def winnow_sensativity():


    return



def winnow_network_connectivity():


    return