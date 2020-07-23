

import biom
import qiime2
import numpy as np
import pandas as pd
import os

from qiime2.plugin import Bool, Str, Int, Float, MetadataColumn, Set

from q2_winnowing.step1_3.Step1_3_Pipeline import main as step1_3_main
from q2_winnowing.step4_5.Step4and5_DecayCurve import main as step4_5_main
from q2_winnowing.step6.Step6_Permanova import main as step6_main
from q2_winnowing.step7_9.Step7_9_Jaccard import main as step7_9_main


def _assemble_artifact_output( combined_metric_df, auc_df, permanova_df, jaccard_df ):

    # Precautionary reset index of dataframes to be joined
    combined_metric_df.reset_index( drop=True, inplace=True )
    jaccard_df.reset_index( drop=True, inplace=True )

    jaccard_kappa = jaccard_df.loc[:,"kappa"] # only new info
    jaccard_agreement = jaccard_df.loc[:,"agreement"] # only new info
    try: # cautionary int index
        instertion_spot = combined_metric_df.columns.get_loc("1")
    except:
        # check for int instead of string value
        instertion_spot = combined_metric_df.columns.get_loc(1)

    # insert jaccard values into dataframe
    combined_metric_df.insert( instertion_spot, "agreement", jaccard_agreement )
    combined_metric_df.insert( instertion_spot, "kappa", jaccard_kappa )

    # Combine output in directory format
    artifact_winnowed = [ ( combined_metric_df, auc_df, permanova_df ) ]
        # this needs to be done since qiime2 has a check for whether len( output views ) == len( semantic types )
        # sadly this fails with tuples so as defined by qiime this is the workaround may be changed in updates
        # len( (1, 2, 3) ) = 3, len( ( (1, 2, 3) ) ) = 3, while len( [ (1, 2, 3) ] ) = 1

    return artifact_winnowed


def _write_to_dump( verbose, dump_path, step ):
    """
    this is simply to write to a dump file if verbose is selected.
    by using this in function it improves readability and removes clutterness of code
    with multiple verbose checks
    :param verbose: if selected will write to dump file
    :param dump_path: the file that is being written to
    :param step: this is which step the program is on, each step corresponds to script ran
    :return: nothing is returned program functions as simple printing function
    """

    if( verbose ):
        with open( dump_path, "a" ) as dump: # allows for closeing if exception is throwns
            if( step == 0 ):
                dump.write("Beginning to convert input to dataframes.\n")
            elif( step == 0.5 ):
                dump.write("Finished converting input to dataframes.\n")
            elif( step == 1 ):
                dump.write("Starting steps 1 to 3\n")
            elif( step == 4 ): # assume starting of step 4 is end of step 1 - 3
                dump.write("Finished steps 1 to 3.\nStarting steps 4 to 5\n")
            elif( step == 6 ): # assume starting of step 6 is end of step 4 - 5
                dump.write("Finished steps 4 to 5.\nStarting step 6\n")
            elif( step == 6.5 ):
                dump.write("Finished step 6.\n")
            elif( step == 7 ):
                dump.write("\nStarting steps 7 to 9\n")
            elif( step == 10 ):
                dump.write(
                    "Output for each winnowing step is written to the respective output folder within each step folder \n"
                    "Example is results form PERMANOVA calculation is written to 'q2_winnowing/step6/output'.\n"
                    "\tThis applies for each step.")
                dump.write("Winnow processing finished.")

    return # Nothing just signifies termination of function


def process(infile1: biom.Table, sample_types: MetadataColumn, metric: Str, conditioning: Str,
               infile2: biom.Table=None, name: Str="-name-", ab_comp: Bool=False, min_count: Int=3,
               total_select: Str="all", iteration_select: Set[Int]=None, pca_components: Int=4,
               smooth_type: Str="sliding_window", window_size: Int=3, centrality: Str=None,
               keep_threshold: Float=0.5, correlation: Str=None, weighted: Bool=False, correlation_prop: Str="both",
               evaluation: Str="kl_divergence", min_connected: Int=0,
               detailed: Bool=False, verbose: Bool=False ) -> list:

    if iteration_select is None: # Since default parameter can't be mutable
        iteration_select = {1, 4, 16, 64, 128}

    outDir = f"{os.path.dirname(os.path.realpath(__file__))}/output"
    # allows for cleaner execution and use of relative paths
    dump = open(f"{outDir}/processing_dump.txt", "w+", encoding="utf-8") # Overwrite file to new empty
    dump.close()
    dump = f"{outDir}/processing_dump.txt"

    _write_to_dump( verbose, dump, step=0)

    # This will be used as part of the PERMANOVA calculation
    sample_types = sample_types.to_dataframe()
    # Make sure input is valid
    num_samples = len( infile1.ids( axis='observation' ) ) # this accounts for abundances being same size as well in later steps
    try:
        if( "type" in sample_types.columns ):
            num_sample_types = len( sample_types.loc[:,"type"] )
        else:
            num_sample_types = len( sample_types.loc[:,"Type"] )
    except:
        raise Exception( "Error: sample metadata must include a column titled Type.")
    if( num_samples != num_sample_types ):
        raise Exception( "Error: each provided sample must have a corresponding type. ( natural/invaded ) ")

    metricOutput = pd.DataFrame() # dataframe to write metrics new
    aucOutput = pd.DataFrame() # Keep most accurate AUC
    permanovaOutput = pd.DataFrame() # Keep most accurate PERMANOVA value
    _write_to_dump( verbose, dump, step=0.5 )

    for iteration_selected in sorted( iteration_select ):

        # Convert input to dataframes
        dataFrame1 = infile1.to_dataframe().to_dense()
        dataFrame1.name = f"{name}_1_{iteration_selected}_"
        dataFrame2 = None
        if (ab_comp):
            dataFrame2 = infile2.to_dataframe().to_dense()
            dataFrame2.name = f"{name}_2_{iteration_selected}_"

        newName = f"{name}_{iteration_selected}_" # will allow for easier iteration selection

        # <><><> Pass data to steps 1 to 3 <><><>
        _write_to_dump(verbose, dump, step=1)
        metric_result, important_features, abundances = \
            _winnow_pipeline( dataFrame1=dataFrame1, dataFrame2=dataFrame2, ab_comp=ab_comp, metric_name=metric,
                              c_type=conditioning, min_count=min_count, total_select=total_select, iteration_select=iteration_selected,
                              pca_components=pca_components, smooth_type=smooth_type, window_size=window_size,
                              centrality_type=centrality, keep_threshold=keep_threshold, correlation=correlation,
                              weighted=weighted, corr_prop=correlation_prop, evaluation_type=evaluation,
                              min_connected=min_connected, detailed=detailed, verbose=verbose)
        # these are used in: Step7_9, Step4_5, Step6

        if( metricOutput.empty ): # create a dataframe of import OTU's for jaccard step
            metricOutput = metric_result
        else:
            if( len(metricOutput.columns) < len(metric_result.columns) ):
                metricOutput.columns = metric_result.columns # this accounts for differing # of OTUs
            metricOutput = metricOutput.append(metric_result, ignore_index=True ) # assign back since does not perform in place

        # <><><> Pass data to steps 4 to 5 <><><>
        _write_to_dump( verbose, dump, step=4 )
        AUC_results, AUC_parameters = \
            _winnow_ordering( dataframe=important_features, name=newName, detailed=detailed, verbose=verbose)
        # these are used in: Step6, None
        aucOutput = AUC_results

        # Note: sample types correspond with abundances being passed
        # print( abundances, AUC_results, sample_types )

        # <><><> Pass data to step 6 <><><>
        _write_to_dump( verbose, dump, step=6 )
        PERMANOVA_results = \
            _winnow_permanova( df_AUC_ordering=AUC_results, df_abundances=abundances, df_samples=sample_types,
                               centralityType=centrality, name=newName, detailed=detailed, verbose=verbose )
        permanovaOutput = PERMANOVA_results
        _write_to_dump( verbose, dump, step=6.5 )


    # <><><>  Pass data to steps 7 to 9 <><><>
    _write_to_dump( verbose, dump, step=7 )
    Jaccard_results = _winnow_sensativity(
        metricOutput, name=f"{metric}_{correlation}_{str(keep_threshold)}_{centrality}_{name}",
        detailed=detailed, verbose=verbose )

    # Notify user of output path
    _write_to_dump( verbose, dump, step=10 )

    # assemble output and return as artifact
    artifact_directory = _assemble_artifact_output( metricOutput, aucOutput, permanovaOutput, Jaccard_results )
    return artifact_directory


def _winnow_pipeline( dataFrame1, dataFrame2, ab_comp: Bool=False, metric_name: Str=None,
                 c_type: Str=None, min_count: Int=3, total_select: Str=None, iteration_select: Str=None,
                 pca_components: Int=4, smooth_type: Str="sliding_window", window_size: Int=3, centrality_type: Str=None,
                 keep_threshold: Float=0.5, correlation: Str=None, weighted: Bool=False, corr_prop: Str="both",
                 evaluation_type: Str=None, min_connected: Int=0, detailed: Bool=False, verbose: Bool=False
                 ):
    """
    Note this function executes the main functionality of steps 1-3 in the pipeline of
    winnowing data.

    :param dataFrame1: This is the biom file which will have OTU info extracted from and analyzed to generate
        an interaction table of taxom.
    :param dataFrame2: This is only used in the case of an A/B analysis and will not be used if ab_comp is False.
    :param ab_comp: Boolean representing whether to perform AB comparison on the data.
    :param metric_name: This is the metric to use.
    :param c_type: Conditioning type to use on the data.
    :param min_count: Features with counts below this number will be removed.
    :param total_select: Number of features to select in total. ie: 1,2,3,... or 'all'
    :param iteration_select: Number of features to select for each time the metric is called. ie: 1,2,3,...
    :param pca_components: Number of pca components to find
    :param smooth_type: Type of Smoothing to be used to remove noise.
    :param window_size: If Smoothing type is a sliding window, this is the size of the window.
    :param centrality_type: If graph_centrality is the metric type, this is the type of Centrality to use.
    :param keep_threshold: If graph_centrality is the metric type, this is the threshold to use to remove weak edges.
    :param correlation: If graph_centrality is the metric type, this is the type of correlation to
        use to build the graph.
    :param weighted: If graph_centrality is the metric type, this specifies if weighted edges should be used
        to create the graph.
    :param corr_prop: If graph centrality is the metric, this specifies if positive, negative, or both types
        of correlation should be used.
    :param evaluation_type: This is the evaluation type to use.
    :param plot_metric: Including this parameter will create a line plot of the metric values for the selected features.
    :param create_graph: If graph centrality is the metric, including this parameter will create a graph image of the
        selected features (using the same correlation type used to select the features).
    :param plot_pca: If PCA is the metric, including this parameter will create a scatter plot image of the first two
        principle components.
    :param naming_file: The file to be used to name the features. If not used, the features will be outputted with
        the names the input file.
    :param proc_id: The identifying number to use in the output file names.
    :param min_connected: The minimum percentage of connectedness of the graph that should be considered before the winnowing process is aborted.
    :param detailed: Notifies plugin to output diagrams and csv files to each steps respective output folder throughout
        computation. If not enabled files will not be generated
    :param verbose: Notifies plugin to generate dump files for every step. These will contain all data that previously
        may have been output through print statements during execution. Each dump.txt file is stored in output folders
        that correspond with each step.
    :return: This is the completed table of different taxa with their corresponding interactions to other taxa
         in their enviroment. This was computed through the use of methods such as: Abundance analysis, AUC analysis,
         F-Score ordering, PERMANOVA calculation, Jacobian matrices, and SEM analysis
    """

    if( total_select and total_select != "all" ):
        try:
            total_select = int( total_select )
        except:
            total_select = "all"

    if( ab_comp ):
        metric_result, important_features, abundances = \
            step1_3_main( dataframe1=dataFrame1, dataframe2=dataFrame2, ab_comp=ab_comp, metric_name=metric_name,
                          c_type=c_type, min_count=min_count, total_select=total_select, iteration_select=iteration_select,
                          pca_components=pca_components, smooth_type=smooth_type, window_size=window_size,
                          centrality_type=centrality_type, keep_threshold=keep_threshold, correlation=correlation,
                          weighted=weighted, corr_prop=corr_prop, evaluation_type=evaluation_type,
                          min_connected=min_connected, detailed=detailed,verbose_p=verbose )
    else:
        metric_result, important_features, abundances = \
            step1_3_main( dataframe1=dataFrame1, dataframe2=None, ab_comp=False, metric_name=metric_name, c_type=c_type,
                          min_count=min_count, total_select=total_select, iteration_select=iteration_select,
                          pca_components=pca_components, smooth_type=smooth_type, window_size=window_size,
                          centrality_type=centrality_type, keep_threshold=keep_threshold, correlation=correlation,
                          weighted=weighted, corr_prop=corr_prop, evaluation_type=evaluation_type,
                          min_connected=min_connected, detailed=detailed,verbose_p=verbose )


    return ( metric_result, important_features, abundances )


def _winnow_ordering( dataframe, name, detailed: Bool=False, verbose: Bool=False ):
    """
    Each OTU is now ordered by centrality and the AUC of each is calculated.
    :param dataframe: Important features and OTUs that have been passed on to be ordered
    :param name: Each file generated by this step will be identified by this name
    :param detailed: Notifies plugin to output diagrams and csv files to each steps respective output folder throughout
        computation. If not enabled files will not be generated
    :param verbose: Notifies plugin to generate dump files for every step. These will contain all data that previously
        may have been output through print statements during execution. Each dump.txt file is stored in output folders
        that correspond with each step.
    :return: A dataframe representing each OTU ordered by centrality and calculated AUC
    """

    # Output files and Parameter files are both generated from this function
    auc_result, auc_param = step4_5_main( dataframe , name=name, detailed=detailed, verbose=verbose )

    return ( auc_result, auc_param )


def _winnow_permanova( df_AUC_ordering, df_abundances, df_samples, centralityType, name, detailed=False, verbose=False ):
    """
    100 Permanovas are ran, Each time adding in more OTUs at a 1% interval according to their order in step 5
    Essentially a sweet spot of additions will be reached and the most influential OTUs will be identified based off
    when additions of OTUs provides little amounts of change in output.
    :param df_AUC_ordering: AUC ordering that was completed earlier
    :param df_abundances: abundances that were generated from initial step
    :param df_samples: Sample that provide data on whether taxom were invade or uninvaded
    :param name: Each file generated by this step will be identified by this name
    :param detailed: Notifies plugin to output diagrams and csv files to each steps respective output folder throughout
        computation. If not enabled files will not be generated
    :param verbose: Notifies plugin to generate dump files for every step. These will contain all data that previously
        may have been output through print statements during execution. Each dump.txt file is stored in output folders
        that correspond with each step.
    :return: a dataframe organised to represent the most influential OTUs
    """

    # Permanova results are generated here
    permanova_result = step6_main( df_AUC_ordering, df_abundances, df_samples, centralityType, name=name,
                                   detailed=detailed, verbose=verbose)

    return permanova_result



def _winnow_sensativity( df_metric_results, name, detailed=False, verbose=False ):
    """
    tested if the number of OTUs selected at each iteration influenced the results. If the results
    were consistent. if so use  leave-one-out analysis to distribute the centrality measure among samples.
    This is because the SEM (step10) requires a centrality measure for each sample that could link to
    environmental variables, whereas the data until this point is centrality for each OTU.
    :param df_metric_results: result of leave one out dataframe
    :param name: Each file generated by this step will be identified by this name
    :param detailed: Notifies plugin to output diagrams and csv files to each steps respective output folder throughout
        computation. If not enabled files will not be generated
    :param verbose:  Notifies plugin to generate dump files for every step. These will contain all data that previously
        may have been output through print statements during execution. Each dump.txt file is stored in output folders
        that correspond with each step.
    :return: A dataframe with output stating measures used as well as the resultant kappa and agreement values.
    """

    jaccard_result = step7_9_main( df_metric_results, name=name, detailed=detailed, verbose=verbose )

    return jaccard_result # Dataframe with Kappa and Agreement values

