

import csv
import subprocess
import pandas as pd
import shutil
import biom

from qiime2.plugin import Bool, Str, Int, Float


from q2_winnowing.step1_3.pipeline import main


def _csv_to_tsv( csvFile ):
    name = csvFile.split("/")[-1].split(".")[0]
    # create a new tmp tsv and convert csv to tsv
    with open( csvFile, "r", encoding="utf-8" ) as csvInFile, open( f"./tmp/tmp_{name}.txt", "w+", newline="", encoding="utf-8") as tsvOutFile:
        csvInFile = csv.reader( csvInFile )
        tsvOutFile = csv.writer( tsvOutFile, delimiter="\t" )

        for row in csvInFile:
            tsvOutFile.writerow( row )
    # new tsv file is stored in a temporary file
    return f"./tmp/tmp_{name}.txt"

# _csv_to_tsv("./bromeA_all.csv")


def _tsv_to_csv( tsvFile ):
    name = tsvFile.split("/")[-1].split(".")[0]
    # create a new tmp tsv and convert csv to tsv
    csvFile = pd.read_table( tsvFile, sep='\t')
    csvFile.to_csv( f"./tmp/tmp_{name}.csv" , index=False)
    # new tsv file is stored in a temporary file
    return f"./tmp/tmp_{name}.txt"

# _tsv_to_csv( "./test_data/test_otu_table.txt")


def _verify_tsv_file( tsvFile ):
    tsvFile = open( tsvFile, "r", encoding="utf-8" )
    observationIDs = tsvFile.readline().strip("\n").split("\t")
    # remove all duplicates by converting to set
    uniqueObservationIDs = set( observationIDs )
    if( len(observationIDs) == len(uniqueObservationIDs) ):
        return True
    else:
        return False


def _tsv_to_biom( tsvFile ):
    # first verify it can be converted
    if( _verify_tsv_file(tsvFile) ):
        name = tsvFile.split("/")[-1].split(".")[0]
        # biom is currently not available as a python package and must be called through bash
        biomConversionCommand = ["biom", "convert", "-i", tsvFile, "-o", f"./tmp/tmp_{name}.biom", "--table-type=OTU table", "--to-hdf5"]
        process = subprocess.Popen( biomConversionCommand, stdout=subprocess.PIPE )
        output, error = process.communicate()

        if( error ):
            raise Exception( "Error: unexpected action when converting tsv to biom" )

        return f"./tmp/tmp_{name}.biom"
    else:
        raise Exception( "Error: tsv file is not in correct format to be converted" )
    # Function is complete

# _tsv_to_biom( "./test_data/test_otu_table.txt")


def _biom_to_tsv( biomFile ):
    # getting name
    name = biomFile.split("/")[-1].split(".")[0]
    # biom is currently not available as a python package and must be called through bash
    biomConversionCommand = ["biom", "convert", "-i", biomFile, "-o", f"./tmp/tmpB_{name}.txt", "--to-tsv"]
    process = subprocess.Popen(biomConversionCommand, stdout=subprocess.PIPE)
    output, error = process.communicate()

    if (error):
        raise Exception("Error: unexpected action when converting biom to tsv")

    removeHash = open(f"./tmp/tmpB_{name}.txt", "r")
    removeHash.readline()
    # this will truncate the file, so need to use a different file name:
    hashRemoved = open(f"./tmp/tmp_{name}.txt", 'w')

    shutil.copyfileobj(removeHash, hashRemoved)

    # finished function
    return f"./tmp/tmp_{name}.txt"


def _csv_to_biom( csvFile ):
    # Convert csv to tsv
    path = _csv_to_tsv( csvFile )
    # Convert tsv to biom
    path = _tsv_to_biom( path )

    # Function completed, data stored in tmp
    return path

def _biom_to_csv( biomFile ):
    # Convert biom to tsv
    path = _biom_to_tsv( biomFile )
    # Convert tsv to csv
    path = _tsv_to_csv( path )

    # Function completed, data is stored in tmp
    return path

# _biom_to_csv( "./tmp/tmp_test_otu_table.biom")



def _winnow_csv( csvFile1, csvFile2=None, ab_comp=False, metric_name=None, c_type=None, min_count=3,
                 total_select=None, iteration_select=None, pca_components=4, smooth_type="sliding_window",
                 window_size=3, centrality_type=None, keep_threshold=0.5, correlation=None,
                 weighted=False, corr_prop="both", evaluation_type=None, plot_metric=False,
                 create_graph=False, plot_pca=False, naming_file=None, proc_id=0, min_connected=0
                 ):

    # Since pipeline takes csv we can just pass values to pipeline and return output values
    main( ab_comp, csvFile1, csvFile2, metric_name, c_type, min_count,
                 total_select, iteration_select, pca_components, smooth_type,
                 window_size, centrality_type, keep_threshold, correlation,
                 weighted, corr_prop, evaluation_type, plot_metric,
                 create_graph, plot_pca, naming_file, proc_id, min_connected
                 )

    return



def _pack_results( metric_name, correlation, keep_threshold, centrality_type ):

    # this is the path of the directory where output is stored
    outDir = f"Results/{metric_name}_{correlation}_{str(keep_threshold)}_{centrality_type}"

    return outDir




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

    csvFile1 = _biom_to_csv( biomFile1 )
    csvFile2 = None

    if( ab_comp ):
        csvFile2 = _biom_to_csv( biomFile2 )

    output = _winnow_csv( csvFile1, csvFile2, ab_comp, metric_name, c_type, min_count,
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