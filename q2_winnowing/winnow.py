

import biom
import numpy as np
import pandas as pd
import os
from pathlib import Path

import qiime2
from qiime2.plugin import Bool, Str, Int, Float, MetadataColumn, Set

# Import different functions from each file
from q2_winnowing.step1_3.pipeline import main as step1_3_main
from q2_winnowing.step4_5.decay_curve import main as step4_5_main
from q2_winnowing.step6.permanova import main as step6_main
from q2_winnowing.step7_9.jaccard import main as step7_9_main


def _assemble_artifact_output( combined_metric_df, auc_df, permanova_df, jaccard_df ):
    """
    This takes care of assembing a Winnowed artifact recognized by qiime2 from the different data
    generated throughout execution.
    :param combined_metric_df: this holds each metric output from the different iteration selections
    :param auc_df: highest iteration auc selection
    :param permanova_df: highest iteration permanova selection
    :param jaccard_df: kappa and agreement values which are measures of consistency of feature ordering
    :return: list( tuple( feature_df, AUC_df, Permanova_df ) ) Note the reasoning for this is explained in return
    """

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


def _verify_output_folders():
    """
    Detailed methods rely on folders being available that can be written to
    with execution steps as well as additive output if the user specifies.

    :return: Nothing simply terminate after all folders are verified
    """

    file_path = f"{os.path.dirname(os.path.realpath(__file__))}" # use relative path never absolute.
    # please see README.md for details on folder structure
    pipeline_output = f"{file_path}/step1_3/output"
    decay_output = f"{file_path}/step4_5/output"
    permanova_output = f"{file_path}/step6/output"
    jaccard_output = f"{file_path}/step7_9/output"

    # Create output folder if one does not already exists
    Path(pipeline_output).mkdir( parents=True, exist_ok=True )
    Path(decay_output).mkdir( parents=True, exist_ok=True )
    Path(permanova_output).mkdir( parents=True, exist_ok=True )
    Path(jaccard_output).mkdir( parents=True, exist_ok=True )

    return # Nothing function is terminated



def _verbose( step ):
    """
    this is simply to print if verbose is selected.
    by using this in function it improves readability and removes clutterness of code
    with multiple verbose checks
    :param step: this is which step the program is on, each step corresponds to script ran
    :return: nothing is returned program functions as simple printing function
    """

    if( step == 0 ):
        print("\nOutput folders verified.\n")
        print("\nBeginning to convert input to dataframes.\n")
    elif( step == 0.5 ):
        print("\nFinished converting input to dataframes.\n")
    elif( step == 1 ):
        print("\nStarting steps 1 to 3 *pipeline*\n")
    elif( step == 1.5 ):
        print("\nUnable to generate feature ordering, moving to next iteration.\n")
    elif( step == 4 ): # assume starting of step 4 is end of step 1 - 3
        print("\nFinished steps 1 to 3 *pipeline*.\nStarting steps 4 to 5 *AUC calculation*\n")
    elif( step == 6 ): # assume starting of step 6 is end of step 4 - 5
        print("\nFinished steps 4 to 5 *AUC calculation*.\nStarting step 6 *PERMANOVA calculation*.\n")
    elif( step == 6.5 ):
        print("\nFinished step 6 *PERMANOVA calculation*.\n")
    elif( step == 7 ):
        print("\nStarting steps 7 to 9 *Jaccard calculation*.\n")
    elif( step == 10 ):
        print(
            "\nOutput for each winnowing step is written to the respective output folder within each step folder \n"
            "Example is results form PERMANOVA calculation is written to 'q2_winnowing/step6/output'.\n"
            "\tThis applies for each step.")
        print("\nWinnow processing finished.\n")

    return # Nothing just signifies termination of function


def _verify_input_is_provided( metric, conditioning, ab_comp, infile2, centrality, correlation ):
    """
    This function provides preventative error checking.
    If a metric is used make sure the parameters for it have been provided.
    Otherwise the process should be terminated immediately in order to prevent wasted time and confusion.
    :param metric: graph_centrality, log_transform, pca_importance, abundance
    :param conditioning: add_one, hellinger
    :param ab_comp: True, False
    :param infile2: Dataframe used only in ab_comp
    :param centrality: betweenness, closeness, degree, eigenvector
    :param correlation: spearman, pearson, kendall, MIC
    :return: Nothing: Success, Otherwise will terminate process with exception
    """
    if( ab_comp ):
        if( infile2 == None ):
            raise Exception( f"Error: Must provide infile2 in order to perform ab_comp. \n"
                             f"Given parameters were ab_comp: {ab_comp}, infile2: {infile2}.")

    if( metric == "graph_centrality" or metric == "log_transform" ):
        if( conditioning == None or centrality == None or correlation == None ):
            raise Exception( f"Error: {metric} requires conditioning, centrality, and correlation parameters. \n"
                             f"Given parameters were {conditioning}, {centrality}, and {correlation}.")
    elif( metric == "pca_importance" ):
        if( conditioning == None ):
            raise Exception( f"Error: A valid conditioning type must be given. Provided conditioning was {conditioning}.")
    else:
        # This is only for when definition is directly called in python, qiime2 will not allow this to be reached.
        raise Exception( f"Error: {metric} is not a valid metric.")
    return # Nothing, just signifies termination of function





def process( infile1: biom.Table, sample_types: MetadataColumn, metric: Str, conditioning: Str,
               infile2: biom.Table=None, name: Str="-name-", ab_comp: Bool=False, min_count: Int=3,
               total_select: Str="all", iteration_select: Set[Int]=None, pca_components: Int=4,
               smooth_type: Str="sliding_window", window_size: Int=3, centrality: Str=None,
               keep_threshold: Float=0.5, correlation: Str=None, weighted: Bool=False, correlation_prop: Str="both",
               evaluation: Str="kl_divergence", min_connected: Int=0,
               detailed: Bool=False ) -> list:
    """
    This is function corresponds with qiime2 function and takes care of file passing between all parts of plugin.

    :param infile1: This is the biom file (qza) which will have OTU info extracted from and analyzed to generate
        an interaction table of taxom.
    :param sample_types: the is metadata representive of samples taken and whether there are invaded/natural
    :param metric: This is the metric to use
    :param conditioning: Conditioning type to use on the data.
    :param infile2: This is only used in the case of an A/B analysis and will not be used if ab_comp is False.
    :param name: This is attached to all detailed output as a means of identification
    :param ab_comp: Boolean representing whether to perform AB comparison on the data.
    :param min_count: Features with counts below this number will be removed.
    :param total_select: Number of features to select in total. ie: 1,2,3,... or 'all'
    :param iteration_select: Number of features to select for each time the metric is called. ie: 1,2,3,...
    :param pca_components: Number of pca components to find
    :param smooth_type:  Type of Smoothing to be used to remove noise.
    :param window_size:  If Smoothing type is a sliding window, this is the size of the window.
    :param centrality: If graph_centrality is the metric type, this is the type of Centrality to use.
    :param keep_threshold: If graph_centrality is the metric type, this is the threshold to use to remove weak edges.
    :param correlation: If graph centrality is the metric, this specifies if positive, negative, or both types
        of correlation should be used.
    :param weighted: If graph_centrality is the metric type, this specifies if weighted edges should be used
        to create the graph.
    :param correlation_prop:
    :param evaluation: This is the evaluation type to use.
    :param min_connected: The minimum percentage of connectedness of the graph that should be considered before the winnowing process is aborted.
    :param detailed: Notifies plugin to output diagrams and csv files to each steps respective output folder throughout
        computation. If not enabled files will not be generated
    :return: return a list of single item with artifact see artifact generation for details on why this is done
    """

    if iteration_select is None: # Since default parameter can't be mutable
        iteration_select = {1, 4, 16, 64, 128}

    # make sure proper file structure is present
    _verify_output_folders()

    out_dir = f"{os.path.dirname(os.path.realpath(__file__))}/output"
    # allows for cleaner execution and use of relative paths

    _verbose( step=0 )

    # This will be used as part of the PERMANOVA calculation
    if( not isinstance( sample_types, pd.DataFrame ) ): # allows for easier testing and input directly to python
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
        raise Exception( "Error: each provided sample must have a corresponding type. ( natural/invaded ).\n"
                         f"Was given {num_samples} samples and {num_sample_types} types. ")

    # Verify parameters are all given
    _verify_input_is_provided( metric, conditioning, ab_comp, infile2, centrality, correlation )

    # if ab_comp is used we will assume that each sample type corresponds with the 1 - n sample of each dataframe
    if( ab_comp ):
        sample_types = pd.concat([sample_types, sample_types], ignore_index=True)

    metric_output = pd.DataFrame() # dataframe to write metrics new
    auc_output = pd.DataFrame() # Keep most accurate AUC
    permanova_output = pd.DataFrame() # Keep most accurate PERMANOVA value
    _verbose( step=0.5 )

    for iteration_selected in sorted( iteration_select ):

        # Convert input to dataframes
        dataframe_1 = infile1.to_dataframe().to_dense()
        dataframe_1.name = f"{name}_1_{iteration_selected}_"
        dataframe_2 = None
        if (ab_comp):
            dataframe_2 = infile2.to_dataframe().to_dense()
            dataframe_2.name = f"{name}_2_{iteration_selected}_"
            if( len(dataframe_1) != len(dataframe_2) ):
                raise Exception(f"Error: Dataframes must be the same size in order to correlate with sample metadata. "
                                f"dataframe1: {len(dataframe_1)} != dataframe2: {len(dataframe_2)}")

        name_new = f"{name}_{iteration_selected}_" # will allow for easier iteration selection

        # <><><> Pass data to steps 1 to 3 <><><>
        _verbose( step=1)
        metric_result, important_features, abundances = \
            _winnow_pipeline( dataframe_1=dataframe_1, dataframe_2=dataframe_2, ab_comp=ab_comp, metric_name=metric,
                              c_type=conditioning, min_count=min_count, total_select=total_select, iteration_select=iteration_selected,
                              pca_components=pca_components, smooth_type=smooth_type, window_size=window_size,
                              centrality_type=centrality, keep_threshold=keep_threshold, correlation=correlation,
                              weighted=weighted, corr_prop=correlation_prop, evaluation_type=evaluation,
                              min_connected=min_connected, detailed=detailed )
        # these are used in: Step7_9, Step4_5, Step6

        if (metric_output.empty):  # create a dataframe of import OTU's for jaccard step
            metric_output = metric_result
        else:
            if (len(metric_output.columns) < len(metric_result.columns)):
                # the dataframe must be extended to be able to hold new data
                new_columns = [col for col in metric_result.columns if not col in metric_output.columns]
                for col in new_columns:
                    metric_output[col] = ""  # Default as empty
            metric_output = pd.concat([metric_output, metric_result], sort=False,
                                      ignore_index=True)  # assign back since does not perform in place

        # check if a metric result was generated before attempting other steps, must be atleast 2 OTUs
        if( 1 in metric_result.columns and 2 in metric_result.columns ):

            # <><><> Pass data to steps 4 to 5 <><><>
            _verbose( step=4 )
            auc_results, auc_parameters = \
                _winnow_ordering( dataframe=important_features, name=name_new, detailed=detailed )
            # these are used in: Step6, None
            auc_output = auc_results

            # Note: sample types correspond with abundances being passed
            # print( abundances, auc_results, sample_types )

            # <><><> Pass data to step 6 <><><>
            _verbose( step=6 )
            permanova_results = \
                _winnow_permanova( auc_ordering_df=auc_results, abundances_df=abundances, samples_df=sample_types,
                                   centrality_type=centrality, name=name_new, detailed=detailed )
            permanova_output = permanova_results
            _verbose( step=6.5 )

        else:
            _verbose( step=1.5)


    # <><><>  Pass data to steps 7 to 9 <><><>
    _verbose( step=7 )
    jaccard_results = _winnow_sensativity(
        metric_output, name=f"{metric}_{correlation}_{str(keep_threshold)}_{centrality}_{name}",
        detailed=detailed )

    # Notify user of output path
    _verbose( step=10 )
    print("############################# DONE #############################")

    # assemble output and return as artifact
    metric_output.replace(r'^\s*$', np.nan, regex=True, inplace=True) # Replace blank with Nan
    artifact_directory = _assemble_artifact_output( metric_output, auc_output, permanova_output, jaccard_results )

    return artifact_directory


def _winnow_pipeline( dataframe_1, dataframe_2, ab_comp: Bool=False, metric_name: Str=None,
                 c_type: Str=None, min_count: Int=3, total_select: Str=None, iteration_select: Str=None,
                 pca_components: Int=4, smooth_type: Str="sliding_window", window_size: Int=3, centrality_type: Str=None,
                 keep_threshold: Float=0.5, correlation: Str=None, weighted: Bool=False, corr_prop: Str="both",
                 evaluation_type: Str=None, min_connected: Int=0, detailed: Bool=False
                 ):
    """
    Note this function executes the main functionality of steps 1-3 in the pipeline of
    winnowing data.

    :param dataframe_1: This is the dataframe which will have OTU info extracted from and analyzed to generate
        an interaction table of taxom.
    :param dataframe_2: This is only used in the case of an A/B analysis and will not be used if ab_comp is False.
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
    :param min_connected: The minimum percentage of connectedness of the graph that should be considered before the winnowing process is aborted.
    :param detailed: Notifies plugin to output diagrams and csv files to each steps respective output folder throughout
        computation. If not enabled files will not be generated
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
            step1_3_main( dataframe1=dataframe_1, dataframe2=dataframe_2, ab_comp=ab_comp, metric_name=metric_name,
                          c_type=c_type, min_count=min_count, total_select=total_select, iteration_select=iteration_select,
                          pca_components=pca_components, smooth_type=smooth_type, window_size=window_size,
                          centrality_type=centrality_type, keep_threshold=keep_threshold, correlation=correlation,
                          weighted=weighted, corr_prop=corr_prop, evaluation_type=evaluation_type,
                          min_connected=min_connected, detailed=detailed )
    else:
        metric_result, important_features, abundances = \
            step1_3_main( dataframe1=dataframe_1, dataframe2=None, ab_comp=False, metric_name=metric_name, c_type=c_type,
                          min_count=min_count, total_select=total_select, iteration_select=iteration_select,
                          pca_components=pca_components, smooth_type=smooth_type, window_size=window_size,
                          centrality_type=centrality_type, keep_threshold=keep_threshold, correlation=correlation,
                          weighted=weighted, corr_prop=corr_prop, evaluation_type=evaluation_type,
                          min_connected=min_connected, detailed=detailed )


    return ( metric_result, important_features, abundances )


def _winnow_ordering( dataframe, name, detailed: Bool=False ):
    """
    Each OTU is now ordered by centrality and the AUC of each is calculated.
    :param dataframe: Important features and OTUs that have been passed on to be ordered
    :param name: Each file generated by this step will be identified by this name
    :param detailed: Notifies plugin to output diagrams and csv files to each steps respective output folder throughout
        computation. If not enabled files will not be generated
    :return: A dataframe representing each OTU ordered by centrality and calculated AUC
    """

    # Output files and Parameter files are both generated from this function
    auc_result, auc_param = step4_5_main( dataframe , name=name, detailed=detailed )

    return ( auc_result, auc_param )


def _winnow_permanova( auc_ordering_df, abundances_df, samples_df, centrality_type, name, detailed=False ):
    """
    100 Permanovas are ran, Each time adding in more OTUs at a 1% interval according to their order in step 5
    Essentially a sweet spot of additions will be reached and the most influential OTUs will be identified based off
    when additions of OTUs provides little amounts of change in output.
    :param auc_ordering_df: AUC ordering that was completed earlier
    :param abundances_df: abundances that were generated from initial step
    :param samples_df: Sample that provide data on whether taxom were invaded or natural
    :param name: Each file generated by this step will be identified by this name
    :param detailed: Notifies plugin to output diagrams and csv files to each steps respective output folder throughout
        computation. If not enabled files will not be generated
    :return: a dataframe organised to represent the most influential OTUs
    """

    # Permanova results are generated here
    permanova_result = step6_main( auc_ordering_df, abundances_df, samples_df, centrality_type, name=name,
                                   detailed=detailed )

    return permanova_result



def _winnow_sensativity( metric_results_df, name, detailed=False ):
    """
    tested if the number of OTUs selected at each iteration influenced the results. If the results
    were consistent. if so use  leave-one-out analysis to distribute the centrality measure among samples.
    This is because the SEM (step10) requires a centrality measure for each sample that could link to
    environmental variables, whereas the data until this point is centrality for each OTU.
    :param metric_results_df: result of leave one out dataframe
    :param name: Each file generated by this step will be identified by this name
    :param detailed: Notifies plugin to output diagrams and csv files to each steps respective output folder throughout
        computation. If not enabled files will not be generated
    :return: A dataframe with output stating measures used as well as the resultant kappa and agreement values.
    """

    jaccard_result = step7_9_main( metric_results_df, name=name, detailed=detailed )

    return jaccard_result # Dataframe with Kappa and Agreement values

