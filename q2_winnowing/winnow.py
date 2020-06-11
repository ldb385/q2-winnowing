

import biom

from qiime2.plugin import Bool, Str, Int, Float

from q2_winnowing._file_conversions import *
from q2_winnowing.step1_3.pipeline import main as step1_3_main
from q2_winnowing.step4_5.Step4and5_DecayCurve import main as step4_5_main



def winnow_pipeline( inFile1: biom.Table, inFile2: biom.Table=None, ab_comp: Bool=False, metric_name: Str=None,
                 c_type: Str=None, min_count: Int=3, total_select: Str=None, iteration_select: Str=None,
                 pca_components: Int=4, smooth_type: Str="sliding_window", window_size: Int=3, centrality_type: Str=None,
                 keep_threshold: Float=0.5, correlation: Str=None, weighted: Bool=False, corr_prop: Str="both",
                 evaluation_type: Str=None, plot_metric: Bool=False, create_graph: Bool=False, plot_pca: Bool=False,
                 naming_file: Str=None, proc_id: Int=0, min_connected: Int=0, verbose: Bool=False
                 ) -> biom.Table:
    # TODO: Implement proper return types
    """
    Note this function executes the main functionality of steps 1-3 in the pipeline of
    winnowing data.

    :param inFile1:
    :param ab_comp:
    :param inFile2:
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

    csvFile1 = biom_to_csv( inFile1 )
    csvFile2 = None

    if( ab_comp ):
        csvFile2 = biom_to_csv( inFile2 )

    output = step1_3_main( csvFile1, csvFile2, ab_comp, metric_name, c_type, min_count,
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



def winnow_ordering( inFile: biom.Table, paramFile: biom.Table, verbose: Bool=False ) -> biom.Table:

    output = step4_5_main( inFile, paramFile, verbose )

    output = csv_to_biom( output )

    return output


def winnow_permanova():


    return



def winnow_sensativity():


    return



def winnow_network_connectivity():


    return