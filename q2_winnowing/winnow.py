

import biom
import qiime2
import numpy as np

from qiime2.plugin import Bool, Str, Int, Float

from q2_winnowing.step1_3.Step1_3_Pipeline import main as step1_3_main
from q2_winnowing.step4_5.Step4and5_DecayCurve import main as step4_5_main

def _dummy_biom_table():
    data = np.arange(40).reshape(10, 4)
    sample_ids = ['S%d' % i for i in range(4)]
    observ_ids = ['O%d' % i for i in range(10)]
    sample_metadata = [{'environment': 'A'}, {'environment': 'B'},
                       {'environment': 'A'}, {'environment': 'B'}]
    observ_metadata = [{'taxonomy': ['Bacteria', 'Firmicutes']},
                       {'taxonomy': ['Bacteria', 'Firmicutes']},
                       {'taxonomy': ['Bacteria', 'Proteobacteria']},
                       {'taxonomy': ['Bacteria', 'Proteobacteria']},
                       {'taxonomy': ['Bacteria', 'Proteobacteria']},
                       {'taxonomy': ['Bacteria', 'Bacteroidetes']},
                       {'taxonomy': ['Bacteria', 'Bacteroidetes']},
                       {'taxonomy': ['Bacteria', 'Firmicutes']},
                       {'taxonomy': ['Bacteria', 'Firmicutes']},
                       {'taxonomy': ['Bacteria', 'Firmicutes']}]
    table = biom.Table(data, observ_ids, sample_ids, observ_metadata,
        sample_metadata, table_id = 'Example Table')

    return table

def winnow_processing( infile1: biom.Table, infile2: biom.Table=None, name: Str="NoNameGiven", ab_comp: Bool=False, metric_name: Str=None,
                 c_type: Str=None, min_count: Int=3, total_select: Str=None, iteration_select: Str=None,
                 pca_components: Int=4, smooth_type: Str="sliding_window", window_size: Int=3, centrality_type: Str=None,
                 keep_threshold: Float=0.5, correlation: Str=None, weighted: Bool=False, corr_prop: Str="both",
                 evaluation_type: Str=None, plot_metric: Bool=False, create_graph: Bool=False, plot_pca: Bool=False,
                 naming_file: Str=None, proc_id: Int=0, min_connected: Int=0, detailed: Bool=False, verbose: Bool=False
                 ) -> biom.Table:
    """

    :param infile1:
    :param infile2:
    :param name:
    :param ab_comp:
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
    :param detailed:
    :param verbose:
    :return:
    """

    # Convert input to dataframes
    dataFrame1 = infile1.to_dataframe()
    dataFrame1.name = f"{name}_1_"
    dataFrame2 = None

    if( ab_comp ):
        dataFrame2 = infile2.to_dataframe()
        dataFrame2.name = f"{name}_2_"

    # Pass data to steps 1 to 3
    metric_result, important_features, abundances = \
        _winnow_pipeline( dataFrame1=dataFrame1, dataFrame2=dataFrame2, ab_comp=ab_comp, metric_name=metric_name,
                          c_type=c_type, min_count=min_count, total_select=total_select, iteration_select=iteration_select,
                          pca_components=pca_components, smooth_type=smooth_type, window_size=window_size,
                          centrality_type=centrality_type, keep_threshold=keep_threshold, correlation=correlation,
                          weighted=weighted, corr_prop=corr_prop, evaluation_type=evaluation_type, plot_metric=plot_metric,
                          create_graph=create_graph, plot_pca=plot_pca, naming_file=naming_file, proc_id=proc_id,
                          min_connected=min_connected, detailed=detailed, verbose=verbose)
    # these are used in: Step7_9, Step4_5, Step6

    # Pass data to steps 4 to 5
    AUC_results, AUC_parameters = \
        _winnow_ordering( dataframe=important_features, name=name, detailed=detailed, verbose=verbose)

    # Pass data to step 6





    return _dummy_biom_table()


def _winnow_pipeline( dataFrame1, dataFrame2, ab_comp: Bool=False, metric_name: Str=None,
                 c_type: Str=None, min_count: Int=3, total_select: Str=None, iteration_select: Str=None,
                 pca_components: Int=4, smooth_type: Str="sliding_window", window_size: Int=3, centrality_type: Str=None,
                 keep_threshold: Float=0.5, correlation: Str=None, weighted: Bool=False, corr_prop: Str="both",
                 evaluation_type: Str=None, plot_metric: Bool=False, create_graph: Bool=False, plot_pca: Bool=False,
                 naming_file: Str=None, proc_id: Int=0, min_connected: Int=0, detailed: Bool=False, verbose: Bool=False
                 ):
    """
    Note this function executes the main functionality of steps 1-3 in the pipeline of
    winnowing data.

    :param dataFrame1:
    :param dataFrame2:
    :param ab_comp:
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
    :param detailed:
    :param verbose:
    :return:
    """

    if( ab_comp ):
        metric_result, important_features, abundances = \
            step1_3_main( dataframe1=dataFrame1, dataframe2=dataFrame2, ab_comp=ab_comp, metric_name=metric_name,
                          c_type=c_type, min_count=min_count, total_select=total_select, iteration_select=iteration_select,
                          pca_components=pca_components, smooth_type=smooth_type, window_size=window_size,
                          centrality_type=centrality_type, keep_threshold=keep_threshold, correlation=correlation,
                          weighted=weighted, corr_prop=corr_prop, evaluation_type=evaluation_type, plot_metric=plot_metric,
                          create_graph=create_graph, plot_pca=plot_pca, naming_file=naming_file, proc_id=proc_id,
                          min_connected=min_connected, detailed=detailed,verbose=verbose )
    else:
        metric_result, important_features, abundances = \
            step1_3_main( dataframe1=dataFrame1, dataframe2=None, ab_comp=False, metric_name=metric_name, c_type=c_type,
                          min_count=min_count, total_select=total_select, iteration_select=iteration_select,
                          pca_components=pca_components, smooth_type=smooth_type, window_size=window_size,
                          centrality_type=centrality_type, keep_threshold=keep_threshold, correlation=correlation,
                          weighted=weighted, corr_prop=corr_prop, evaluation_type=evaluation_type, plot_metric=plot_metric,
                          create_graph=create_graph, plot_pca=plot_pca, naming_file=naming_file, proc_id=proc_id,
                          min_connected=min_connected, detailed=detailed,verbose=verbose )


    return ( metric_result, important_features, abundances )


def _winnow_ordering( dataframe, name, detailed: Bool=False, verbose: Bool=False ):

    # Output files and Parameter files are both generated from this function
    output_result, output_param = step4_5_main( dataframe , name=name, detailed=detailed, verbose=verbose )

    return ( output_result, output_param )


def _winnow_permanova():


    return



def _winnow_sensativity():


    return



def _winnow_network_connectivity():


    return

