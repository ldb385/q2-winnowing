

import os
import csv
import subprocess
import pandas as pd
import shutil
tmpPath = f"{os.path.dirname(os.path.realpath(__file__))}"


def _csv_to_tsv( csvFile ):
    name = csvFile.split("/")[-1].split(".")[0]
    # create a new tmp tsv and convert csv to tsv
    with open( csvFile, "r", encoding="utf-8" ) as csvInFile, open( f"{tmpPath}/tmp/tmp_{name}.txt", "w+", newline="", encoding="utf-8") as tsvOutFile:
        csvInFile = csv.reader( csvInFile )
        tsvOutFile = csv.writer( tsvOutFile, delimiter="\t" )

        for row in csvInFile:
            tsvOutFile.writerow( row )
    # new tsv file is stored in a temporary file
    return f"{tmpPath}/tmp/tmp_{name}.txt"

# _csv_to_tsv("./bromeA_all.csv")


def _tsv_to_csv( tsvFile ):
    name = tsvFile.split("/")[-1].split(".")[0]
    # create a new tmp tsv and convert csv to tsv
    csvFile = pd.read_table( tsvFile, sep='\t')
    csvFile.to_csv( f"{tmpPath}/tmp/tmp_{name}.csv" , index=False)
    # new tsv file is stored in a temporary file
    return f"{tmpPath}/tmp/tmp_{name}.txt"

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
        biomConversionCommand = ["biom", "convert", "-i", tsvFile, "-o", f"{tmpPath}/tmp/tmp_{name}.biom", "--table-type=OTU table", "--to-hdf5"]
        process = subprocess.Popen( biomConversionCommand, stdout=subprocess.PIPE )
        output, error = process.communicate()

        if( error ):
            raise Exception( "Error: unexpected action when converting tsv to biom" )

        return f"{tmpPath}/tmp/tmp_{name}.biom"
    else:
        raise Exception( "Error: tsv file is not in correct format to be converted" )
    # Function is complete

# _tsv_to_biom( "./test_data/test_otu_table.txt")


def _biom_to_tsv( biomFile ):
    # getting name
    name = biomFile.split("/")[-1].split(".")[0]
    # biom is currently not available as a python package and must be called through bash
    biomConversionCommand = ["biom", "convert", "-i", biomFile, "-o", f"{tmpPath}/tmp/tmpB_{name}.txt", "--to-tsv"]
    process = subprocess.Popen(biomConversionCommand, stdout=subprocess.PIPE)
    output, error = process.communicate()

    if (error):
        raise Exception("Error: unexpected action when converting biom to tsv")

    removeHash = open(f"{tmpPath}/tmp/tmpB_{name}.txt", "r")
    removeHash.readline()
    # this will truncate the file, so need to use a different file name:
    hashRemoved = open(f"{tmpPath}/tmp/tmp_{name}.txt", 'w')

    shutil.copyfileobj(removeHash, hashRemoved)

    # finished function
    return f"{tmpPath}/tmp/tmp_{name}.txt"


def csv_to_biom( csvFile ):
    # Convert csv to tsv
    path = _csv_to_tsv( csvFile )
    # Convert tsv to biom
    path = _tsv_to_biom( path )

    # Function completed, data stored in tmp
    return path

# csv_to_biom( "C:/Users/liamb/VirtualSharedFolder/q2-winnowing/TEST/bromeA_all.csv")

def biom_to_csv( biomFile ):
    # Convert biom to tsv
    path = _biom_to_tsv( biomFile )
    # Convert tsv to csv
    path = _tsv_to_csv( path )

    # Function completed, data is stored in tmp
    return path

# _biom_to_csv( "./tmp/tmp_test_otu_table.biom")

def _pack_results( metric_name, correlation, keep_threshold, centrality_type ):

    # this is the path of the directory where output is stored
    outDir = f"{tmpPath}/Results/{metric_name}_{correlation}_{str(keep_threshold)}_{centrality_type}"

    return outDir